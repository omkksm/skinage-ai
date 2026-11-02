/* SkinAge-Omics (브라우저판): DEG · 피부노화 패널 · GO/KEGG ORA · 기존 분석 비교 · 기전 역전 */
"use strict";
const Omics = {};
Omics.PANEL = {"ECM 분해 (MMP)": {"MMP1": "up", "MMP2": "up", "MMP3": "up", "MMP9": "up", "MMP12": "up", "MMP13": "up"}, "콜라겐·ECM 합성": {"COL1A1": "down", "COL1A2": "down", "COL3A1": "down", "COL4A1": "down", "COL7A1": "down", "COL17A1": "down", "ELN": "var", "FBN1": "down", "DCN": "down", "LUM": "down", "HAS2": "down", "HAS3": "down"}, "MMP 억제 (TIMP)": {"TIMP1": "var", "TIMP2": "var", "TIMP3": "down"}, "세포노화·세포주기": {"CDKN2A": "up", "CDKN1A": "up", "TP53": "up", "LMNB1": "down", "GLB1": "up", "SERPINE1": "up"}, "SASP·염증": {"IL6": "up", "CXCL8": "up", "IL1A": "up", "IL1B": "up", "TNF": "up", "CCL2": "up", "NFKB1": "up", "PTGS2": "up"}, "산화스트레스·항산화": {"SOD1": "down", "SOD2": "var", "CAT": "down", "GPX1": "down", "NFE2L2": "down", "HMOX1": "var", "NQO1": "var"}, "표피 장벽·분화": {"FLG": "down", "LOR": "down", "IVL": "var", "KRT10": "down", "KRT14": "down", "AQP3": "down", "CLDN1": "down", "TGM1": "var"}, "장수·대사 신호": {"SIRT1": "down", "SIRT3": "down", "SIRT6": "down", "FOXO3": "down", "MTOR": "up", "IGF1": "down", "IGF1R": "var", "PPARGC1A": "down"}, "TGF-β·AP-1 신호": {"TGFB1": "var", "TGFBR2": "down", "SMAD3": "down", "SMAD7": "up", "JUN": "up", "FOS": "up", "CCN2": "down", "CCN1": "up"}, "색소·광노화": {"MITF": "var", "TYR": "var", "DCT": "var", "PMEL": "var"}};
Omics.ALIAS = { IL8: "CXCL8", CTGF: "CCN2", CYR61: "CCN1", P16: "CDKN2A", P21: "CDKN1A", NRF2: "NFE2L2", PGC1A: "PPARGC1A" };
Omics.panelRows = () => Object.entries(Omics.PANEL).flatMap(([cat, g]) => Object.entries(g).map(([gene, dir]) => ({ gene, cat, dir })));
Omics.cache = {};
Omics.fetchText = path => Omics.cache[path] ? Promise.resolve(Omics.cache[path]) : fetch(path).then(r => { if (!r.ok) throw new Error(path + " 불러오기 실패"); return r.text(); }).then(t => (Omics.cache[path] = t));
Omics.GS = { "KEGG (MSigDB v7.5.1)": "data/kegg.gmt.txt", "GO Biological Process 2021": "data/go_bp.gmt.txt", "MSigDB Hallmark v7.0": "data/hallmark.gmt.txt" };
Omics.readGMT = text => { const sets = {}; text.split(/\r?\n/).forEach(l => { const p = l.split("\t"); if (p.length > 2) sets[p[0]] = new Set(p.slice(2).filter(Boolean).map(g => g.split(",")[0].trim().toUpperCase())); }); return sets; };

/* ---------- 결과표 표준화 */
Omics.COL = { gene: ["gene", "genes", "symbol", "gene_symbol", "genesymbol", "gene_name", "id", "external_gene_name", ""], lfc: ["log2foldchange", "logfc", "log2fc", "lfc"], p: ["pvalue", "p.value", "pval", "p"], padj: ["padj", "fdr", "adj.p.val", "qvalue"], base: ["basemean", "aveexpr", "logcpm"] };
Omics.standardizeDEG = text => {
  const { header, rows } = S.parseCSV(text), h = header.map(x => x.toLowerCase()), f = {};
  for (const [k, c] of Object.entries(Omics.COL)) { const i = h.findIndex(x => c.includes(x)); if (i >= 0) f[k] = i; }
  if (f.gene == null) f.gene = 0;
  if (f.lfc == null) throw new Error("log2 fold change 컬럼(log2FoldChange, logFC)을 찾지 못했습니다.");
  const seen = new Set(), out = [];
  rows.forEach(r => { let g = (r[f.gene] || "").trim().toUpperCase(); g = Omics.ALIAS[g] || g; const lfc = parseFloat(r[f.lfc]); if (!g || seen.has(g) || !Number.isFinite(lfc)) return; seen.add(g); out.push({ gene: g, lfc, p: f.p != null ? parseFloat(r[f.p]) : NaN, padj: f.padj != null ? parseFloat(r[f.padj]) : NaN, base: f.base != null ? parseFloat(r[f.base]) : NaN }); });
  if (f.padj == null && f.p != null) { const q = S.bh(out.map(d => Number.isFinite(d.p) ? d.p : 1)); out.forEach((d, i) => d.padj = q[i]); }
  return out;
};
Omics.classify = (deg, padj, lfc) => deg.forEach(d => d.status = (d.padj < padj && Math.abs(d.lfc) >= lfc) ? (d.lfc > 0 ? "Up" : "Down") : "NS");

/* ---------- count matrix → DEG (median-of-ratios + NB-Wald / 무반복 시 고정 dispersion) */
Omics.parseCounts = text => { const { header, rows } = S.parseCSV(text), samples = header.slice(1), genes = [], M = []; const agg = {};
  rows.forEach(r => { const g = (r[0] || "").trim(); if (!g) return; const v = samples.map((_, j) => Math.max(0, Math.round(parseFloat(r[j + 1]) || 0))); if (agg[g]) agg[g].forEach((x, j) => agg[g][j] = x + v[j]); else agg[g] = v; });
  for (const [g, v] of Object.entries(agg)) { genes.push(g); M.push(v); } return { samples, genes, M }; };
Omics.parseMeta = text => { const { header, rows } = S.parseCSV(text); return { cols: header.slice(1), rows: rows.map(r => ({ sample: r[0].trim(), v: Object.fromEntries(header.slice(1).map((h, i) => [h, (r[i + 1] || "").trim()])) })) }; };
Omics.deg = (cnt, groupOf, ref, test, minCount = 10) => {
  const idx = cnt.samples.map((s, j) => [s, j]).filter(([s]) => groupOf[s] === ref || groupOf[s] === test), cols = idx.map(x => x[1]), isT = idx.map(([s]) => groupOf[s] === test);
  if (cols.length < 2) throw new Error("비교할 샘플이 2개 미만입니다.");
  let G = [], Y = [];
  cnt.M.forEach((row, i) => { const y = cols.map(j => row[j]); if (y.reduce((a, b) => a + b, 0) >= minCount) { G.push(cnt.genes[i]); Y.push(y); } });
  const n = cols.length, logs = Y.filter(y => y.every(v => v > 0)).map(y => y.map(Math.log));
  if (logs.length < 10) throw new Error("모든 샘플에서 0이 아닌 유전자가 너무 적습니다.");
  const sf = [...Array(n).keys()].map(j => Math.exp(S.median(logs.map(l => l[j] - S.mean(l)))));
  const nT = isT.filter(Boolean).length, nR = n - nT, noRep = Math.min(nT, nR) < 2;
  const norm = Y.map(y => y.map((v, j) => v / sf[j]));
  let alpha;
  if (noRep) alpha = norm.map(() => 0.01);
  else { /* 그룹 내 적률추정 → 추세 수축 */
    const disp = norm.map(x => { let acc = 0, w = 0; for (const g of [true, false]) { const xs = x.filter((_, j) => isT[j] === g), ss = sf.filter((_, j) => isT[j] === g); if (xs.length < 2) continue; const m = S.mean(xs), v = S.sd(xs) ** 2, inv = S.mean(ss.map(s => 1 / s)); acc += (m > 0 ? (v - m * inv) / (m * m) : 0) * (xs.length - 1); w += xs.length - 1; } return Math.min(Math.max(acc / Math.max(w, 1), 1e-8), 10); });
    const mean = norm.map(x => S.mean(x)), ok = mean.map(m => m > 1); let sxx = 0, sx = 0, sy = 0, sxy = 0, k = 0; mean.forEach((m, i) => { if (!ok[i]) return; const xx = 1 / m; sx += xx; sy += disp[i]; sxx += xx * xx; sxy += xx * disp[i]; k++; });
    const a1 = Math.max((k * sxy - sx * sy) / (k * sxx - sx * sx), 1e-6), a0 = Math.max((sy - a1 * sx) / k, 1e-6), dof = n - 2, sh = dof / (dof + 10);
    alpha = disp.map((d, i) => { const tr = a0 + a1 / Math.max(mean[i], 1e-3); return Math.min(Math.max(Math.exp(sh * Math.log(Math.max(d, 1e-8)) + (1 - sh) * Math.log(tr)), tr), 10); });
  }
  const res = Y.map((y, i) => { const mu = {}, se = {};
    for (const [lab, g] of [["r", false], ["t", true]]) { const yy = y.filter((_, j) => isT[j] === g), ss = sf.filter((_, j) => isT[j] === g); let m = Math.max(yy.reduce((a, b) => a + b, 0) / ss.reduce((a, b) => a + b, 0), 1e-8);
      for (let it = 0; it < 25; it++) { let sc = 0, inf = 0; yy.forEach((v, k) => { const smu = ss[k] * m, den = 1 + alpha[i] * smu; sc += (v - smu) / den; inf += smu / den; }); const st = Math.max(-5, Math.min(5, sc / Math.max(inf, 1e-12))); m *= Math.exp(st); if (Math.abs(st) < 1e-8) break; }
      m = Math.max(m, 1e-8); let inf = 0; ss.forEach(s => { const smu = s * m; inf += smu / (1 + alpha[i] * smu); }); mu[lab] = m; se[lab] = 1 / Math.sqrt(Math.max(inf, 1e-12)); }
    const l = Math.log(mu.t) - Math.log(mu.r), z = l / Math.sqrt(se.t ** 2 + se.r ** 2);
    return { gene: (Omics.ALIAS[G[i].toUpperCase()] || G[i].toUpperCase()), lfc: Math.max(-12, Math.min(12, l / Math.LN2)), p: S.p2(z), base: S.mean(norm[i]) }; });
  const q = S.bh(res.map(r => r.p)); res.forEach((r, i) => r.padj = q[i]);
  return { deg: res, noRep, n, sf };
};

/* ---------- 패널 · ORA · 비교 · 역전 */
Omics.panel = deg => { const m = Object.fromEntries(deg.map(d => [d.gene, d])); return Omics.panelRows().map(r => { const d = m[r.gene]; const st = d ? d.status : "미검출"; const conc = (!d || st === "NS" || r.dir === "var") ? "-" : ((st === "Up") === (r.dir === "up") ? "일치" : "불일치"); return { ...r, lfc: d ? d.lfc : NaN, padj: d ? d.padj : NaN, st, conc }; }); };
Omics.ora = (genes, universe, sets, minS = 10, maxS = 500) => { const U = new Set(universe), G = new Set(genes.filter(g => U.has(g))), N = U.size, n = G.size, out = [];
  for (const [term, gs] of Object.entries(sets)) { let K = 0, hits = []; gs.forEach(g => { if (U.has(g)) { K++; if (G.has(g)) hits.push(g); } }); if (K < minS || K > maxS || !hits.length) continue; out.push({ term: term.replace(/^KEGG_/, "").replace(/_/g, " "), k: hits.length, K, fe: (hits.length / n) / (K / N), p: S.hyperSf(hits.length, N, K, n), genes: hits.sort().join(", ") }); }
  const q = S.bh(out.map(o => o.p)); out.forEach((o, i) => o.padj = q[i]); return out.sort((a, b) => a.p - b.p); };
Omics.compare = (a, b) => { const mb = Object.fromEntries(b.map(d => [d.gene, d])), common = a.filter(d => mb[d.gene]), x = common.map(d => d.lfc), y = common.map(d => mb[d.gene].lfc);
  const sa = new Set(a.filter(d => d.status !== "NS").map(d => d.gene)), sb = new Set(b.filter(d => d.status !== "NS").map(d => d.gene)), both = [...sa].filter(g => sb.has(g));
  const same = both.filter(g => a.find(d => d.gene === g).status === mb[g].status).length, union = new Set([...sa, ...sb]).size;
  const ct = {}; common.forEach(d => { const k = d.status + "|" + mb[d.gene].status; ct[k] = (ct[k] || 0) + 1; });
  return { n: common.length, r: S.pearson(x, y), rho: S.spearman(x, y), da: sa.size, db: sb.size, both: both.length, jac: both.length / union, dir: same / Math.max(both.length, 1), rec: both.length / Math.max(sb.size, 1), prec: both.length / Math.max(sa.size, 1), ct, pts: common.map((d, i) => [y[i], x[i], d.gene]) }; };
Omics.reversal = deg => { const m = Object.fromEntries(deg.map(d => [d.gene, d])), used = []; Omics.panelRows().forEach(r => { if (r.dir === "var" || !m[r.gene]) return; const sg = r.dir === "up" ? 1 : -1, v = -sg * Math.tanh(m[r.gene].lfc); used.push({ gene: r.gene, cat: r.cat, dir: r.dir, lfc: m[r.gene].lfc, rev: v }); }); return { score: used.length ? S.mean(used.map(u => u.rev)) : NaN, used }; };
/* 엑소좀 실데이터: 원자료 보호 — 미리 계산한 요약(유의 DEG·패널 유전자·그래프 표본·통계)만 사용 */
Omics.exoSummary = () => Omics.fetchText("data/exosome_summary.json").then(t => JSON.parse(t));
Omics.exoDEG = sm => sm.genes.map(([gene, lfc, padj, status]) => ({ gene, lfc, padj, status })).concat(sm.ns_sample.map(([lfc, padj]) => ({ gene: "", lfc, padj, status: "NS" })));
Omics.loadExoDEG = () => Omics.exoSummary().then(sm => sm.genes.map(([gene, lfc, padj, status]) => ({ gene, lfc, padj, status })));

/* ---------- 화면 */
Omics.st = { src: "demo", padj: 0.05, lfc: 1, tab: "deg", gs: "KEGG (MSigDB v7.5.1)" };
Omics.init = root => { Omics.root = root; Omics.render(); };
Omics.load = async () => {
  const st = Omics.st;
  if (st.src === "exo") { const sm = await Omics.exoSummary(); return { deg: Omics.exoDEG(sm), noRep: sm.noRep, n: sm.n, sm, protected: true, label: "연구실 실데이터: 마우스 색소세포 α-MSH vs EV17 엑소좀 5 ppm (각 1샘플) · 원자료 보호 모드(결과 열람만 가능, 기준 FDR 0.05·1.5배 고정)" }; }
  if (st.src === "demo") {
    const [c, m] = ["data/demo_counts.csv.txt", "data/demo_metadata.csv.txt"];
    const cnt = Omics.parseCounts(await Omics.fetchText(c)), meta = Omics.parseMeta(await Omics.fetchText(m)), gOf = Object.fromEntries(meta.rows.map(r => [r.sample, r.v.group]));
    const r = Omics.deg(cnt, gOf, "young", "aged"); return { deg: r.deg, noRep: r.noRep, label: "합성 데모 데이터: 젊은 피부 6 vs 노화 피부 6 (정답이 알려진 시뮬레이션)", ref: "young", test: "aged" };
  }
  if (st.src === "upcounts") { if (!st.cnt || !st.meta) return null; const gcol = st.meta.cols.find(c => c.toLowerCase() === "group") || st.meta.cols[0], gOf = Object.fromEntries(st.meta.rows.map(r => [r.sample, r.v[gcol]])), lv = [...new Set(Object.values(gOf))].sort(); const ref = st.ref && lv.includes(st.ref) ? st.ref : lv[0], test = st.test && lv.includes(st.test) && st.test !== ref ? st.test : lv.find(l => l !== ref); const r = Omics.deg(st.cnt, gOf, ref, test); return { deg: r.deg, noRep: r.noRep, label: `업로드 데이터: ${ref} vs ${test}`, ref, test, levels: lv }; }
  if (st.src === "updeg") return st.updeg ? { deg: st.updeg, label: "업로드한 DEG 결과표" } : null;
};
Omics.render = async () => {
  const st = Omics.st, R = Omics.root;
  R.innerHTML = `<div class="split"><aside class="card ctrl"><h3>입력 데이터</h3>
    ${[["demo", "데모 데이터 (합성)"], ["exo", "실데이터 예시: 엑소좀 RNA-seq"], ["upcounts", "Count matrix 업로드"], ["updeg", "DEG 결과표 업로드"]].map(([k, l]) => `<label class="radio"><input type="radio" name="om-src" value="${k}" ${k === st.src ? "checked" : ""}> ${l}</label>`).join("")}
    <div id="om-up"></div>
    <h3>DEG 기준</h3><label for="om-padj">FDR 기준</label><input id="om-padj" type="number" step="0.01" min="0.0001" max="0.5" value="${st.padj}"><label for="om-lfc">|log2FC| 기준 (1 = 2배, 0.585 = 1.5배)</label><input id="om-lfc" type="number" step="0.05" min="0" max="5" value="${st.lfc}">
    <p class="hint">브라우저판은 핵심 분석을 제공합니다. AI 노화점수(LOOCV)와 전체 기능은 Python(Streamlit)판에서 사용할 수 있습니다.</p></aside>
    <section><nav class="subtabs">${[["deg", "① DEG 시각화"], ["panel", "② 피부노화 유전자"], ["ora", "③ GO·KEGG 경로"], ["cmp", "④ 기존 분석과 비교"]].map(([k, l]) => `<button data-tab="${k}" aria-selected="${k === st.tab}">${l}</button>`).join("")}</nav><div id="om-body"><p class="loading card">분석 중…</p></div></section></div>`;
  R.querySelectorAll("input[name=om-src]").forEach(r => r.onchange = () => { st.src = r.value; if (r.value === "exo") st.lfc = 0.585; if (r.value === "demo") st.lfc = 1; Omics.render(); });
  R.querySelectorAll("[data-tab]").forEach(b => b.onclick = () => { st.tab = b.dataset.tab; Omics.render(); });
  R.querySelector("#om-padj").onchange = e => { st.padj = +e.target.value; Omics.render(); };
  R.querySelector("#om-lfc").onchange = e => { st.lfc = +e.target.value; Omics.render(); };
  const up = R.querySelector("#om-up");
  if (st.src === "upcounts") { up.innerHTML = `<label class="btn" for="om-c">Count matrix (행=유전자, 열=샘플)</label><input id="om-c" type="file" accept=".csv,.tsv,.txt" hidden><label class="btn" for="om-m">샘플 정보 (첫 열=샘플명, group 열)</label><input id="om-m" type="file" accept=".csv,.tsv,.txt" hidden><p class="hint">${st.cnt ? "count ✓ " : ""}${st.meta ? "샘플정보 ✓" : ""}</p>`;
    up.querySelector("#om-c").onchange = e => e.target.files[0].text().then(t => { st.cnt = Omics.parseCounts(t); Omics.render(); });
    up.querySelector("#om-m").onchange = e => e.target.files[0].text().then(t => { st.meta = Omics.parseMeta(t); Omics.render(); }); }
  if (st.src === "updeg") { up.innerHTML = `<label class="btn" for="om-d">DESeq2/edgeR/limma 결과표</label><input id="om-d" type="file" accept=".csv,.tsv,.txt" hidden>`; up.querySelector("#om-d").onchange = e => e.target.files[0].text().then(t => { try { st.updeg = Omics.standardizeDEG(t); Omics.render(); } catch (er) { up.insertAdjacentHTML("beforeend", `<p class="err">${esc(er.message)}</p>`); } }); }
  const body = R.querySelector("#om-body");
  let D; try { D = await Omics.load(); } catch (e) { body.innerHTML = `<p class="err card">${esc(e.message)}</p>`; return; }
  if (!D) { body.innerHTML = `<p class="card">왼쪽에서 파일을 올려 주세요.</p>`; return; }
  const deg = D.deg; if (!D.protected) Omics.classify(deg, st.padj, st.lfc);
  if (D.protected) { R.querySelector("#om-padj").disabled = true; R.querySelector("#om-lfc").disabled = true; R.querySelector("#om-padj").value = 0.05; R.querySelector("#om-lfc").value = 0.585; }
  const named = deg.filter(d => d.gene);
  const nU = deg.filter(d => d.status === "Up").length, nD = deg.filter(d => d.status === "Down").length, pn = Omics.panel(named), pnSig = pn.filter(p => p.st === "Up" || p.st === "Down").length;
  let h = `<p class="hint">${esc(D.label)}${D.noRep ? " · 그룹당 반복 1개 → 고정 dispersion(BCV 0.1, edgeR 무반복 권장값) 사용, 탐색적 결과" : ""}</p>
    <div class="kpis"><div class="kpi"><span>검정 유전자</span><b>${(D.n || deg.length).toLocaleString()}</b></div><div class="kpi"><span>증가 DEG</span><b style="color:var(--up)">${nU}</b></div><div class="kpi"><span>감소 DEG</span><b style="color:var(--down)">${nD}</b></div><div class="kpi"><span>피부노화 패널 유의 변화</span><b>${pnSig} / ${pn.length}</b></div></div>`;
  if (st.tab === "deg") {
    const top = deg.filter(d => d.status !== "NS").sort((a, b) => a.padj - b.padj);
    h += `<div class="card">${C.volcano(deg, { padj: st.padj, lfc: st.lfc, highlight: ["MMP1", "COL1A1", "CDKN2A", "IL6"] })}</div>
      <div class="card tblwrap"><div class="row"><h3 style="margin:0">유의 DEG (${top.length})</h3>${D.protected ? '<span class="hint">🔒 원자료 보호: 다운로드 불가</span>' : '<button class="btn" id="om-dl">전체 결과 CSV 저장</button>'}</div><table><thead><tr><th>유전자</th><th class="n">log2FC</th><th class="n">FDR</th><th>판정</th></tr></thead><tbody>${top.slice(0, 150).map(d => `<tr><td>${d.gene}</td><td class="n">${C.fmt(d.lfc, 2)}</td><td class="n">${d.padj.toExponential(2)}</td><td>${d.status === "Up" ? "증가" : "감소"}</td></tr>`).join("")}</tbody></table>${top.length > 150 ? '<p class="hint">상위 150개만 표시. 전체는 CSV로 저장하세요.</p>' : ""}</div>`;
    body.innerHTML = h; if (!D.protected) body.querySelector("#om-dl").onclick = () => App.download("deg_results.csv", S.toCSV(["gene", "log2fc", "pvalue", "padj", "status"], deg.map(d => [d.gene, d.lfc, d.p, d.padj, d.status])), "text/csv");
  } else if (st.tab === "panel") {
    const cats = [...new Set(pn.map(p => p.cat))];
    h += `<div class="card tblwrap"><table><thead><tr><th>범주</th><th class="n">패널</th><th class="n">검출</th><th class="n">유의 변화</th><th class="n">방향 일치</th><th class="n">불일치</th><th class="n">평균 log2FC</th></tr></thead><tbody>${cats.map(c => { const g = pn.filter(p => p.cat === c), l = g.filter(p => Number.isFinite(p.lfc)).map(p => p.lfc); return `<tr><td>${c}</td><td class="n">${g.length}</td><td class="n">${g.filter(p => p.st !== "미검출").length}</td><td class="n">${g.filter(p => p.st === "Up" || p.st === "Down").length}</td><td class="n">${g.filter(p => p.conc === "일치").length}</td><td class="n">${g.filter(p => p.conc === "불일치").length}</td><td class="n">${l.length ? C.fmt(S.mean(l), 2) : "–"}</td></tr>`; }).join("")}</tbody></table></div>
      <div class="card">${C.hbars(pn.filter(p => Number.isFinite(p.lfc)).map(p => ({ k: `${p.gene} ${p.dir === "up" ? "▲" : p.dir === "down" ? "▼" : "◆"}`, v: p.lfc, c: p.st === "Up" ? "var(--up)" : p.st === "Down" ? "var(--down)" : "var(--ns)" })), { dec: 2, xlab: "log2 Fold Change (▲증가·▼감소 문헌 경향, ◆상황의존, 색 = 유의 판정)", left: 130 })}</div>`;
    body.innerHTML = h;
  } else if (st.tab === "ora") {
    if (D.protected) {
      const o = D.sm.ora[st.gs], tbl = (r, dir, col) => `<div><h3>${dir} DEG 과대표현 경로 (FDR<0.05: ${r.filter(x => x.padj < 0.05).length})</h3>${r.length ? C.hbars(r.slice(0, 12).map(x => ({ k: x.term.length > 46 ? x.term.slice(0, 44) + "…" : x.term, v: -Math.log10(x.padj), c: col })), { dec: 1, xlab: "−log10(FDR)", left: 300, W: 640 }) : "<p>과대표현 경로 없음</p>"}<div class="tblwrap"><table><thead><tr><th>경로</th><th class="n">겹침</th><th class="n">크기</th><th class="n">FE</th><th class="n">FDR</th></tr></thead><tbody>${r.map(x => `<tr><td>${esc(x.term)}</td><td class="n">${x.k}</td><td class="n">${x.K}</td><td class="n">${C.fmt(x.fe, 1)}</td><td class="n">${x.padj.toExponential(1)}</td></tr>`).join("")}</tbody></table></div></div>`;
      body.innerHTML = h + `<div class="card"><label for="om-gs">유전자세트</label> <select id="om-gs">${Object.keys(Omics.GS).map(k => `<option ${k === st.gs ? "selected" : ""}>${k}</option>`).join("")}</select><p class="hint">배경 = 검정 유전자 ${D.n.toLocaleString()}개 · 초기하분포 + BH-FDR (원자료 보호 모드: 미리 계산된 결과)</p><div class="grid2">${tbl(o.Up, "증가", "var(--up)")}${tbl(o.Down, "감소", "var(--down)")}</div></div>`;
      body.querySelector("#om-gs").onchange = e => { st.gs = e.target.value; Omics.render(); };
      return;
    }
    body.innerHTML = h + `<div class="card"><label for="om-gs">유전자세트</label> <select id="om-gs">${Object.keys(Omics.GS).map(k => `<option ${k === st.gs ? "selected" : ""}>${k}</option>`).join("")}</select><div id="om-ora"><p class="loading">유전자세트 불러오는 중…</p></div></div>`;
    body.querySelector("#om-gs").onchange = e => { st.gs = e.target.value; Omics.render(); };
    const sets = Omics.readGMT(await Omics.fetchText(Omics.GS[st.gs])), U = deg.map(d => d.gene);
    const part = (dir, col) => { const r = Omics.ora(deg.filter(d => d.status === dir).map(d => d.gene), U, sets); return `<div><h3>${dir === "Up" ? "증가" : "감소"} DEG 과대표현 경로 (FDR<0.05: ${r.filter(x => x.padj < 0.05).length})</h3>${r.length ? C.hbars(r.slice(0, 12).map(x => ({ k: x.term.length > 46 ? x.term.slice(0, 44) + "…" : x.term, v: -Math.log10(x.padj), c: col })), { dec: 1, xlab: "−log10(FDR)", left: 300, W: 640 }) : "<p>과대표현 경로 없음</p>"}<div class="tblwrap"><table><thead><tr><th>경로</th><th class="n">겹침</th><th class="n">크기</th><th class="n">FE</th><th class="n">FDR</th></tr></thead><tbody>${r.slice(0, 25).map(x => `<tr><td title="${esc(x.genes)}">${esc(x.term)}</td><td class="n">${x.k}</td><td class="n">${x.K}</td><td class="n">${C.fmt(x.fe, 1)}</td><td class="n">${x.padj.toExponential(1)}</td></tr>`).join("")}</tbody></table></div></div>`; };
    body.querySelector("#om-ora").innerHTML = `<p class="hint">배경 = 검정된 유전자 ${U.length.toLocaleString()}개 · 초기하분포 검정 + BH-FDR</p><div class="grid2">${part("Up", "var(--up)")}${part("Down", "var(--down)")}</div>`;
  } else if (st.tab === "cmp" && D.protected) {
    const c = D.sm.cmp, k = ["Up", "NS", "Down"];
    body.innerHTML = h + `<div class="card"><h3>기존 edgeR 결과(2024.11)와 비교</h3><div class="kpis"><div class="kpi"><span>공통 유전자</span><b>${c.n.toLocaleString()}</b></div><div class="kpi"><span>log2FC Spearman ρ</span><b>${C.fmt(c.rho, 3)}</b></div><div class="kpi"><span>기존 DEG 재현율</span><b class="good">${C.fmt(c.rec * 100, 1)}%</b></div><div class="kpi"><span>공통 DEG 방향 일치</span><b class="good">${C.fmt(c.dir * 100, 0)}%</b></div></div>
      <div class="grid2"><div>${C.xy({ series: [{ pts: c.pts, color: "var(--accent)", r: 1.6, op: .5 }], diag: true, xlab: "log2FC — 기존 edgeR", ylab: "log2FC — SkinAge-Omics" })}<p class="hint">Pearson r ${C.fmt(c.r, 3)} · 그래프는 1,200개 유전자 무작위 표본</p></div>
      <div class="tblwrap"><table><thead><tr><th>SkinAge-Omics \\ 기존 edgeR</th>${k.map(x => `<th class="n">${x}</th>`).join("")}</tr></thead><tbody>${k.map(a => `<tr><td>${a}</td>${k.map(b => `<td class="n">${(c.ct[a + "|" + b] || 0).toLocaleString()}</td>`).join("")}</tr>`).join("")}</tbody></table><p>앱 DEG ${c.da} · 기존 DEG ${c.db} · 공통 ${c.both} · Jaccard ${C.fmt(c.jac, 3)} · 앱 DEG 중 기존과 일치 ${C.fmt(c.prec * 100, 1)}%</p></div></div></div>`;
  } else if (st.tab === "cmp") {
    body.innerHTML = h + `<div class="card"><div class="row"><label class="btn" for="om-cmp">기존 DEG 결과표 올리기</label><input id="om-cmp" type="file" accept=".csv,.tsv,.txt" hidden></div><div id="om-cmpout"><p class="hint">기존에 수행한 DESeq2/edgeR/limma 결과와 log2FC 상관, DEG 중복도, 방향 일치율을 비교합니다.</p></div></div>`;
    const show = (other, lab) => { Omics.classify(other, st.padj, st.lfc); const c = Omics.compare(deg, other), k = ["Up", "NS", "Down"];
      body.querySelector("#om-cmpout").innerHTML = `<div class="kpis"><div class="kpi"><span>공통 유전자</span><b>${c.n.toLocaleString()}</b></div><div class="kpi"><span>log2FC Spearman ρ</span><b>${C.fmt(c.rho, 3)}</b></div><div class="kpi"><span>${lab} DEG 재현율</span><b class="good">${C.fmt(c.rec * 100, 1)}%</b></div><div class="kpi"><span>공통 DEG 방향 일치</span><b class="good">${C.fmt(c.dir * 100, 0)}%</b></div></div>
        <div class="grid2"><div>${C.xy({ series: [{ pts: c.pts.filter((_, i) => i % 2 === 0), color: "var(--accent)", r: 1.6, op: .5 }], diag: true, xlab: `log2FC — ${lab}`, ylab: "log2FC — SkinAge-Omics" })}<p class="hint">Pearson r ${C.fmt(c.r, 3)} (한쪽 count 0 유전자의 처리 차이 포함)</p></div>
        <div class="tblwrap"><table><thead><tr><th>SkinAge-Omics \\ ${esc(lab)}</th>${k.map(x => `<th class="n">${x}</th>`).join("")}</tr></thead><tbody>${k.map(a => `<tr><td>${a}</td>${k.map(b => `<td class="n">${(c.ct[a + "|" + b] || 0).toLocaleString()}</td>`).join("")}</tr>`).join("")}</tbody></table><p>앱 DEG ${c.da} · ${esc(lab)} DEG ${c.db} · 공통 ${c.both} · Jaccard ${C.fmt(c.jac, 3)}</p></div></div>`; };
    body.querySelector("#om-cmp").onchange = e => e.target.files[0].text().then(t => { try { show(Omics.standardizeDEG(t), "기존 분석"); } catch (er) { body.querySelector("#om-cmpout").innerHTML = `<p class="err">${esc(er.message)}</p>`; } });
  }
};
