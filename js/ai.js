/* SkinAge-AI 예측 엔진 + 화면 (core/predict.py 와 같은 모델의 브라우저 구현) */
"use strict";
const AI = {};

AI.ORG = {
  cell: { label: "세포 (in vitro)", sens: 1.0, deliver: 1.0, dur: [1, 7], unit: "µM", cmax: 100, ecs: 1.0 },
  mouse: { label: "마우스 (in vivo)", sens: 0.75, deliver: 0.75, dur: [14, 84], unit: "% (도포)", cmax: 2, ecs: 0.02 },
  human: { label: "인체 (적용시험)", sens: 0.45, deliver: 0.55, dur: [28, 168], unit: "% (도포)", cmax: 2, ecs: 0.025 }
};
AI.DIS = {
  aging: { label: "피부노화", short: "노화 지수", axes: { ECM: 1.0, COL: 1.2, SEN: 1.0, INF: 0.6, OX: 0.8, PIG: 0.6, BAR: 0.5, MOR: 1.2 }, stresses: ["uv", "ox", "gly"] },
  atopic: { label: "아토피 피부염", short: "아토피 중증도 지수", axes: { TH2: 1.2, IGE: 0.8, ITCH: 1.0, BAR: 1.2, INF: 0.9, MOR: 0.5, OX: 0.4 }, stresses: ["allergen", "barrier", "ox"] }
};
AI.STR = {
  uv: { label: "광노화 (UVB/UVA)", unit: { cell: "UVB 10–50 mJ/cm²", mouse: "UVB 60–180 mJ/cm², 주 3회", human: "일광 누적노출 지수" },
    w: { ECM: 1.0, COL: 0.8, SEN: 0.5, INF: 0.9, OX: 0.8, PIG: 1.0, BAR: 0.7, MOR: 0.9, TH2: 0.15, IGE: 0.05, ITCH: 0.2 } },
  ox: { label: "산화스트레스 (H₂O₂ 등)", unit: { cell: "H₂O₂ 50–400 µM", mouse: "산화 유도 모델 강도", human: "산화스트레스 노출 지수" },
    w: { ECM: 0.6, COL: 0.5, SEN: 1.0, INF: 0.6, OX: 1.0, PIG: 0.5, BAR: 0.4, MOR: 0.5, TH2: 0.25, IGE: 0.1, ITCH: 0.25 } },
  gly: { label: "당화(AGEs)·내인성 노화", unit: { cell: "AGEs/고당 또는 반복계대", mouse: "D-gal·고당 식이·연령", human: "연령·당화 지수" },
    w: { ECM: 0.4, COL: 1.0, SEN: 0.9, INF: 0.5, OX: 0.6, PIG: 0.3, BAR: 0.5, MOR: 0.8, TH2: 0.05, IGE: 0, ITCH: 0.1 } },
  allergen: { label: "알레르겐·햅텐 감작 (DNCB·HDM·MC903)", unit: { cell: "TNF-α/IFN-γ 또는 IL-4/IL-13 자극 (HaCaT)", mouse: "DNCB 반복 도포 (BALB/c·NC/Nga)", human: "알레르겐 노출·증상 지수" },
    w: { ECM: 0.3, COL: 0.2, SEN: 0.2, INF: 0.8, OX: 0.4, PIG: 0.2, BAR: 0.8, MOR: 0.7, TH2: 1.0, IGE: 1.0, ITCH: 0.9 } },
  barrier: { label: "피부장벽 손상 (계면활성제·테이프 스트리핑)", unit: { cell: "SDS 처리 (HaCaT)", mouse: "SDS·아세톤 반복 처리", human: "세정제·건조 환경 노출 지수" },
    w: { ECM: 0.2, COL: 0.2, SEN: 0.2, INF: 0.6, OX: 0.3, PIG: 0.1, BAR: 1.0, MOR: 0.4, TH2: 0.5, IGE: 0.3, ITCH: 0.6 } }
};
AI.AXES = { ECM: "ECM 분해", COL: "콜라겐 합성", SEN: "세포노화", INF: "염증", OX: "산화스트레스", PIG: "색소침착", BAR: "피부장벽", MOR: "피부 형태", TH2: "Th2 면역반응", IGE: "IgE 감작", ITCH: "가려움" };
AI.BM = [
  ["MMP-1", "ECM", 1, ["cell", "mouse"], 180], ["Type I collagen (COL1A1·PIP)", "COL", -1, ["cell", "mouse", "human"], 55],
  ["p16INK4a", "SEN", 1, ["cell", "mouse"], 150], ["SA-β-gal 양성률", "SEN", 1, ["cell"], 200], ["IL-6", "INF", 1, ["cell", "mouse"], 220],
  ["ROS", "OX", 1, ["cell", "mouse"], 160], ["멜라닌", "PIG", 1, ["cell", "mouse", "human"], 90], ["TEWL", "BAR", 1, ["mouse", "human"], 110],
  ["Filaggrin", "BAR", -1, ["cell", "mouse", "human"], 60], ["주름 깊이", "MOR", 1, ["mouse", "human"], 70], ["피부 탄력 (R2)", "MOR", -1, ["human"], 40],
  ["표피 두께", "MOR", 1, ["mouse"], 120], ["TSLP", "TH2", 1, ["cell", "mouse"], 250], ["IL-4·IL-13", "TH2", 1, ["mouse", "human"], 220],
  ["혈청 IgE", "IGE", 1, ["mouse", "human"], 300], ["긁기 횟수·가려움 점수", "ITCH", 1, ["mouse", "human"], 200]
].map(([name, axis, dir, orgs, max]) => ({ name, axis, dir, orgs, max }));
AI.BMI = Object.fromEntries(AI.BM.map(b => [b.name, b]));
AI.COMP = {
  "후보물질 A (항산화형)": { ec50: 12, hill: 1.2, E: { ECM: .5, COL: .3, SEN: .55, INF: .4, OX: .9, PIG: .45, BAR: .2, MOR: .35, TH2: .2, IGE: .1, ITCH: .2 } },
  "후보물질 B (항염·장벽형)": { ec50: 20, hill: 1.0, E: { ECM: .45, COL: .2, SEN: .25, INF: .9, OX: .3, PIG: .2, BAR: .8, MOR: .4, TH2: .55, IGE: .3, ITCH: .6 } },
  "후보물질 C (콜라겐·항노화형)": { ec50: 8, hill: 1.5, E: { ECM: .7, COL: .85, SEN: .6, INF: .3, OX: .25, PIG: .1, BAR: .3, MOR: .7, TH2: .05, IGE: 0, ITCH: .1 } },
  "후보물질 D (미백형)": { ec50: 15, hill: 1.3, E: { ECM: .15, COL: .1, SEN: .15, INF: .25, OX: .45, PIG: .9, BAR: .1, MOR: .1, TH2: .05, IGE: 0, ITCH: .05 } },
  "후보물질 E (면역조절·항알레르기형)": { ec50: 10, hill: 1.2, E: { ECM: .2, COL: .1, SEN: .1, INF: .7, OX: .2, PIG: .05, BAR: .45, MOR: .45, TH2: .9, IGE: .8, ITCH: .85 } }
};
AI.TIM = { pre: ["전처리 (예방)", 1.0], co: ["동시 처리", 0.8], post: ["후처리 (개선)", 0.5] };

/* ---------- 기전 시뮬레이터 */
AI.truth = (c) => {
  const o = AI.ORG[c.org], s = AI.STR[c.stress], [d0, d1] = o.dur;
  const dn = Math.min(Math.max((c.days - d0) / Math.max(d1 - d0, 1), 0), 1.5), durF = 1 - Math.exp(-2.2 * (0.25 + dn)), levF = c.level / (c.level + 0.35) * 1.35;
  const out = {};
  for (const b of AI.BM) {
    if (!b.orgs.includes(c.org)) continue;
    const dmg = b.max * o.sens * s.w[b.axis] * levF * durF;
    let prot = 0; const cp = AI.COMP[c.comp];
    if (cp && c.conc > 0) { const ec = cp.ec50 * o.ecs, h = c.conc ** cp.hill / (ec ** cp.hill + c.conc ** cp.hill); prot = cp.E[b.axis] * h * AI.TIM[c.timing][1] * o.deliver; }
    out[b.name] = b.dir * dmg * (1 - prot);
  }
  return out;
};
AI.simulate = (n = 1800, seed = 7) => {
  const r = S.rng(seed), rows = [], comps = [...Object.keys(AI.COMP), "없음"], orgs = Object.keys(AI.ORG), strs = Object.keys(AI.STR);
  for (let i = 0; i < n; i++) {
    const org = r.choice(orgs), stress = r.choice(strs, [0.18, 0.2, 0.17, 0.25, 0.2]), level = r.choice([0.25, 0.5, 0.75, 1.0]);
    const [d0, d1] = AI.ORG[org].dur, days = Math.round(r.choice([0, 1, 2, 3].map(k => d0 + (d1 - d0) * k / 3)));
    const comp = r.choice(comps, [0.17, 0.17, 0.17, 0.17, 0.17, 0.15]);
    const conc = comp === "없음" ? 0 : r.choice([0.03, 0.1, 0.3, 1.0]) * AI.ORG[org].cmax;
    const timing = comp === "없음" ? "pre" : r.choice(Object.keys(AI.TIM));
    const c = { org, stress, level, days, comp, conc, timing }, t = AI.truth(c);
    for (const [b, v] of Object.entries(t)) rows.push({ ...c, source: "시뮬레이션", bm: b, value: v + r.norm(0, 0.12 * Math.abs(v) + 2), truth: v });
  }
  return rows;
};

/* ---------- 특징 인코딩 */
AI.encode = (c, comps) => {
  const f = [];
  for (const o of Object.keys(AI.ORG)) f.push(c.org === o ? 1 : 0);
  for (const s of Object.keys(AI.STR)) f.push(c.stress === s ? c.level : 0);
  f.push(c.level);
  const [d0, d1] = AI.ORG[c.org].dur; f.push((c.days - d0) / Math.max(d1 - d0, 1));
  const lc = Math.log10(c.conc / AI.ORG[c.org].cmax * 100 + 1) / 2;
  for (const k of comps) f.push(c.comp === k ? lc : 0);
  for (const t of Object.keys(AI.TIM)) f.push(c.timing === t && c.conc > 0 ? lc : 0);
  return f;
};
AI.groups = comps => { const g = []; const push = (name, n) => { for (let i = 0; i < n; i++) g.push(name); }; push("개체", 3); push("스트레스 종류·강도", 6); push("노출 기간", 1); push("후보물질·농도", comps.length); push("처리 시점", 3); return g; };

/* ---------- 모델 */
AI.Model = class {
  constructor(data, maxN = 350, seed = 0) {
    this.comps = [...new Set(data.map(d => d.comp))].filter(c => c !== "없음").sort();
    this.models = {}; this.train = {};
    const r = S.rng(seed), by = {};
    data.forEach(d => (by[d.bm] = by[d.bm] || []).push(d));
    for (const [b, rows] of Object.entries(by)) {
      if (rows.length < 15) continue;
      let use = rows;
      if (rows.length > maxN) { const real = rows.filter(d => d.source !== "시뮬레이션"), sim = rows.filter(d => d.source === "시뮬레이션"); const p = r.perm(sim.length).slice(0, Math.max(maxN - real.length, 50)); use = real.concat(p.map(i => sim[i])); }
      const g = S.fitGP(use.map(d => AI.encode(d, this.comps)), use.map(d => d.value));
      if (g) { this.models[b] = g; this.train[b] = use; }
    }
  }
  predict(c, b) { return this.models[b].predict(AI.encode(c, this.comps)); }
};

AI.profile = (m, c, disease) => {
  const rows = [];
  for (const b of AI.BM) {
    if (!b.orgs.includes(c.org) || !m.models[b.name]) continue;
    if (disease && !(b.axis in AI.DIS[disease].axes)) continue;
    const [s] = m.predict({ ...c, comp: "없음", conc: 0, timing: "pre" }, b.name), [w, sd] = m.predict(c, b.name);
    rows.push({ bm: b.name, axis: b.axis, s, w, ci: 1.96 * sd, p: Math.abs(s) > 10 ? (s - w) / s * 100 : NaN });
  }
  return rows;
};
AI.score = (prof, disease) => { const ax = AI.DIS[disease].axes; let n = 0, d = 0; prof.forEach(r => { if (Number.isFinite(r.p) && ax[r.axis]) { n += ax[r.axis] * Math.max(-50, Math.min(120, r.p)); d += ax[r.axis]; } }); return d ? Math.max(-50, Math.min(120, n / d)) : NaN; };
AI.dindex = (vals, org, disease) => { const ax = AI.DIS[disease].axes; let n = 0, d = 0; for (const [b, v] of Object.entries(vals)) { const bm = AI.BMI[b]; if (!ax[bm.axis]) continue; n += ax[bm.axis] * Math.min(Math.max(v * bm.dir / (bm.max * AI.ORG[org].sens), 0), 1.2); d += ax[bm.axis]; } return d ? Math.min(100, n / d * 100) : NaN; };
AI.disease = (prof, org, disease) => { const s = AI.dindex(Object.fromEntries(prof.map(r => [r.bm, r.s])), org, disease), w = AI.dindex(Object.fromEntries(prof.map(r => [r.bm, r.w])), org, disease); return { s, w, p: s > 1 ? (s - w) / s * 100 : NaN }; };
AI.dose = (m, c, disease, n = 12) => { const cm = AI.ORG[c.org].cmax, out = [[0, 0]]; for (let i = 0; i < n; i++) { const x = cm * 0.005 * (200 ** (i / (n - 1))); out.push([x, AI.score(AI.profile(m, { ...c, conc: x }, disease), disease)]); } return out; };
AI.fit4pl = pts => { const p = pts.filter(q => q[0] > 0 && Number.isFinite(q[1])), cs = p.map(q => q[0]), rs = p.map(q => q[1]); let best = null; const lo = Math.log10(Math.min(...cs) / 3), hi = Math.log10(Math.max(...cs) * 3);
  for (let i = 0; i <= 100; i++) { const le = lo + (hi - lo) * i / 100; for (let h = 0.5; h <= 3.001; h += 0.1) { const g = cs.map(c => c ** h / ((10 ** le) ** h + c ** h)); const gg = g.reduce((s, v) => s + v * v, 0), gr = g.reduce((s, v, k) => s + v * rs[k], 0), top = gr / gg; const sse = rs.reduce((s, v, k) => s + (v - top * g[k]) ** 2, 0); if (!best || sse < best.sse) best = { sse, ec50: 10 ** le, hill: h, top }; } }
  best.rec80 = best.ec50 * 4 ** (1 / best.hill); return best; };

AI.cv = (data, k = 5, maxN = 250, seed = 1) => {
  const r = S.rng(seed), by = {}, comps = [...new Set(data.map(d => d.comp))].filter(c => c !== "없음").sort(), res = [], pairs = [];
  data.forEach(d => (by[d.bm] = by[d.bm] || []).push(d));
  for (const [b, rows0] of Object.entries(by)) {
    if (rows0.length < 30) continue;
    const rows = rows0.length > maxN ? r.perm(rows0.length).slice(0, maxN).map(i => rows0[i]) : rows0;
    const fold = r.perm(rows.length).map(i => i % k), pred = new Array(rows.length), sd = new Array(rows.length);
    for (let f = 0; f < k; f++) {
      const tr = rows.filter((_, i) => fold[i] !== f), g = S.fitGP(tr.map(d => AI.encode(d, comps)), tr.map(d => d.value));
      rows.forEach((d, i) => { if (fold[i] === f) [pred[i], sd[i]] = g.predict(AI.encode(d, comps)); });
    }
    const y = rows.map(d => d.value), my = S.mean(y), ss = y.reduce((s, v) => s + (v - my) ** 2, 0), se = y.reduce((s, v, i) => s + (v - pred[i]) ** 2, 0);
    let r2t = NaN; if (rows.every(d => Number.isFinite(d.truth))) { const t = rows.map(d => d.truth), mt = S.mean(t); r2t = 1 - t.reduce((s, v, i) => s + (v - pred[i]) ** 2, 0) / t.reduce((s, v) => s + (v - mt) ** 2, 0); }
    res.push({ bm: b, n: rows.length, r2: 1 - se / ss, rmse: Math.sqrt(se / rows.length), cover: y.filter((v, i) => Math.abs(v - pred[i]) <= 1.96 * sd[i]).length / y.length, r2t });
    pairs.push({ bm: b, pts: y.map((v, i) => [v, pred[i]]) });
  }
  return { res, pairs };
};
AI.importance = (m, b, seed = 2) => {
  const r = S.rng(seed), rows = m.train[b], X = rows.map(d => AI.encode(d, m.comps)), y = rows.map(d => d.value), g = m.models[b], groups = AI.groups(m.comps);
  const mse = XX => XX.reduce((s, x, i) => s + (g.predict(x)[0] - y[i]) ** 2, 0) / XX.length, base = mse(X), out = {};
  for (const name of [...new Set(groups)]) { const idx = groups.map((n, i) => n === name ? i : -1).filter(i => i >= 0); let inc = 0; for (let rep = 0; rep < 2; rep++) { const p = r.perm(X.length); inc += mse(X.map((x, i) => x.map((v, j) => idx.includes(j) ? X[p[i]][j] : v))) - base; } out[name] = Math.max(inc / 2, 0); }
  const tot = Object.values(out).reduce((s, v) => s + v, 0) || 1; return Object.entries(out).map(([k, v]) => ({ k, v: v / tot * 100 })).sort((a, b) => b.v - a.v);
};

/* 실데이터·문헌 양식 */
AI.TEMPLATE = [["organism", "model", "stress", "level", "days", "compound", "conc", "timing", "biomarker", "value", "pmid", "note"],
  ["cell", "HDF", "uv", 0.5, 2, "물질X", 10, "pre", "MMP-1", 85, "", "value = 무처리 대조 대비 % 변화 (증가 +, 감소 −)"],
  ["cell", "HDF", "uv", 0.5, 2, "없음", 0, "pre", "MMP-1", 160, "", "스트레스 단독군도 함께 입력"],
  ["mouse", "HR-1", "uv", 0.75, 56, "물질X", 0.5, "co", "주름 깊이", 30, "", "level: 약 0.25 / 중 0.5 / 강 0.75 / 최강 1.0"],
  ["mouse", "BALB/c", "allergen", 0.75, 28, "물질Y", 1.0, "co", "혈청 IgE", 120, "12345678", "문헌값은 PMID 입력 (예시 번호)"]];
AI.parseReal = text => {
  const { header, rows } = S.parseCSV(text), h = header.map(x => x.toLowerCase()), need = ["organism", "stress", "level", "days", "compound", "conc", "timing", "biomarker", "value"];
  const miss = need.filter(x => !h.includes(x)); if (miss.length) throw new Error("필수 컬럼이 없습니다: " + miss.join(", "));
  const ix = Object.fromEntries(h.map((x, i) => [x, i])), out = [], bad = [];
  rows.forEach((r, n) => {
    const g = k => (r[ix[k]] || "").trim(), d = { org: g("organism").toLowerCase(), stress: g("stress").toLowerCase(), level: +g("level"), days: +g("days"), comp: g("compound"), conc: +g("conc"), timing: g("timing").toLowerCase(), bm: g("biomarker"), value: +g("value") };
    const pm = ix.pmid != null ? g("pmid") : "";
    if (!AI.ORG[d.org] || !AI.STR[d.stress] || !AI.TIM[d.timing] || !AI.BMI[d.bm] || ![d.level, d.days, d.conc, d.value].every(Number.isFinite)) { bad.push(n + 2); return; }
    out.push({ ...d, source: /^\d{5,9}$/.test(pm) ? "문헌" : "실데이터", pmid: pm, truth: NaN });
  });
  return { out, bad };
};

/* ---------- 화면 */
AI.state = { disease: "aging", org: "cell", stress: "uv", level: 0.75, days: 4, comp: null, conc: 30, timing: "pre", tab: "sim", extra: [] };
AI.init = root => {
  AI.root = root;
  root.innerHTML = `<div class="loading card">AI 예측 모델을 학습하는 중입니다… (기전 시뮬레이션 데이터 1,800조건 · 지표 16개)</div>`;
  setTimeout(() => { AI.sim = AI.simulate(); AI.retrain(); AI.render(); }, 30);
};
AI.retrain = () => { const t0 = performance.now(); AI.data = AI.sim.concat(AI.state.extra); AI.model = new AI.Model(AI.data); AI.trainMs = performance.now() - t0; AI.cvRes = null; if (!AI.state.comp || !AI.model.comps.includes(AI.state.comp)) AI.state.comp = AI.model.comps[2] || AI.model.comps[0]; };
AI.cond = () => { const s = AI.state; return { org: s.org, stress: s.stress, level: s.level, days: s.days, comp: s.comp, conc: s.conc, timing: s.timing }; };

AI.render = () => {
  const s = AI.state, o = AI.ORG[s.org], m = AI.model;
  if (!AI.DIS[s.disease].stresses.includes(s.stress)) s.stress = AI.DIS[s.disease].stresses[0];
  s.days = Math.min(Math.max(s.days, o.dur[0]), o.dur[1]); s.conc = Math.min(s.conc, o.cmax);
  const opt = (obj, cur, lab = v => v) => Object.keys(obj).map(k => `<option value="${k}" ${k === cur ? "selected" : ""}>${esc(lab(k))}</option>`).join("");
  const nReal = AI.data.filter(d => d.source !== "시뮬레이션").length;
  AI.root.innerHTML = `
  <div class="kpis"><div class="kpi"><span>예측 지표</span><b>${Object.keys(m.models).length}개</b></div><div class="kpi"><span>학습 측정값</span><b>${AI.data.length.toLocaleString()}건</b></div><div class="kpi"><span>실험·문헌 데이터</span><b>${nReal.toLocaleString()}건</b></div><div class="kpi"><span>후보물질</span><b>${m.comps.length}종</b></div></div>
  <div class="split">
   <aside class="card ctrl" aria-label="예측 조건">
    <h3>예측 조건</h3>
    <div class="seg" role="radiogroup" aria-label="질환">${Object.keys(AI.DIS).map(d => `<button data-dis="${d}" aria-pressed="${d === s.disease}">${AI.DIS[d].label}</button>`).join("")}</div>
    <label for="ai-org">개체</label><select id="ai-org">${opt(AI.ORG, s.org, k => AI.ORG[k].label)}</select>
    <label for="ai-str">환경·스트레스 요인</label><select id="ai-str">${AI.DIS[s.disease].stresses.map(k => `<option value="${k}" ${k === s.stress ? "selected" : ""}>${esc(AI.STR[k].label)}</option>`).join("")}</select>
    <p class="hint">강도 기준: ${esc(AI.STR[s.stress].unit[s.org])}</p>
    <label for="ai-lev">스트레스 강도 <b>${{ 0.25: "약", 0.5: "중", 0.75: "강", 1: "최강" }[s.level]}</b></label><input id="ai-lev" type="range" min="1" max="4" step="1" value="${s.level * 4}">
    <label for="ai-days">노출 기간 <b>${s.days}일</b></label><input id="ai-days" type="range" min="${o.dur[0]}" max="${o.dur[1]}" step="1" value="${s.days}">
    <label for="ai-comp">후보물질</label><select id="ai-comp">${m.comps.map(c => `<option ${c === s.comp ? "selected" : ""}>${esc(c)}</option>`).join("")}</select>
    <label for="ai-conc">농도 <b>${C.fmt(s.conc, 3)} ${o.unit}</b></label><input id="ai-conc" type="range" min="0" max="${o.cmax}" step="${o.cmax / 100}" value="${s.conc}">
    <label>처리 시점</label><div class="seg">${Object.keys(AI.TIM).map(t => `<button data-tim="${t}" aria-pressed="${t === s.timing}">${AI.TIM[t][0]}</button>`).join("")}</div>
   </aside>
   <section>
    <nav class="subtabs" role="tablist">${[["sim", "① 예측 시뮬레이터"], ["dose", "② 용량-반응 · 최적 농도"], ["rank", "③ 후보물질 비교"], ["perf", "④ 모델 성능 · 근거"], ["hub", "⑤ 데이터 허브"], ["rev", "⑥ 전사체 기전 역전"]].map(([k, l]) => `<button role="tab" data-tab="${k}" aria-selected="${k === s.tab}">${l}</button>`).join("")}</nav>
    <div id="ai-body"></div>
   </section>
  </div>`;
  const q = sel => AI.root.querySelector(sel);
  AI.root.querySelectorAll("[data-dis]").forEach(b => b.onclick = () => { s.disease = b.dataset.dis; AI.render(); });
  AI.root.querySelectorAll("[data-tim]").forEach(b => b.onclick = () => { s.timing = b.dataset.tim; AI.render(); });
  AI.root.querySelectorAll("[data-tab]").forEach(b => b.onclick = () => { s.tab = b.dataset.tab; AI.render(); });
  q("#ai-org").onchange = e => { const prev = AI.ORG[s.org]; s.org = e.target.value; const n = AI.ORG[s.org]; s.conc = s.conc / prev.cmax * n.cmax; s.days = Math.round(n.dur[0] + (s.days - prev.dur[0]) / (prev.dur[1] - prev.dur[0]) * (n.dur[1] - n.dur[0])); AI.render(); };
  q("#ai-str").onchange = e => { s.stress = e.target.value; AI.render(); };
  q("#ai-comp").onchange = e => { s.comp = e.target.value; AI.render(); };
  const live = (id, f) => { const el = q(id); el.oninput = () => { f(+el.value); el.previousElementSibling.querySelector("b").textContent = id === "#ai-lev" ? { 0.25: "약", 0.5: "중", 0.75: "강", 1: "최강" }[s.level] : id === "#ai-days" ? s.days + "일" : C.fmt(s.conc, 3) + " " + o.unit; }; el.onchange = () => AI.body(); };
  live("#ai-lev", v => s.level = v / 4); live("#ai-days", v => s.days = v); live("#ai-conc", v => s.conc = v);
  AI.body();
};

AI.body = () => {
  const s = AI.state, m = AI.model, o = AI.ORG[s.org], el = AI.root.querySelector("#ai-body"), c = AI.cond(), lab = AI.DIS[s.disease].short;
  const condTxt = `${o.label} · ${AI.STR[s.stress].label} · 강도 ${s.level} · ${s.days}일 · ${s.comp} ${C.fmt(s.conc, 3)} ${o.unit} · ${AI.TIM[s.timing][0]}`;
  if (s.tab === "sim") {
    const prof = AI.profile(m, c, s.disease), dz = AI.disease(prof, s.org, s.disease), sc = AI.score(prof, s.disease);
    const axes = [...new Set(prof.map(r => r.axis))], ra = axes.map(a => S.mean(prof.filter(r => r.axis === a).map(r => Math.abs(r.s)))), rb = axes.map(a => S.mean(prof.filter(r => r.axis === a).map(r => Math.abs(r.w))));
    el.innerHTML = `<div class="kpis"><div class="kpi"><span>${lab} (스트레스 단독)</span><b>${C.fmt(dz.s, 0)} / 100</b></div><div class="kpi"><span>${lab} (+ 후보물질)</span><b class="good">${C.fmt(dz.w, 0)} / 100</b></div><div class="kpi"><span>질환 예방률</span><b class="good">${C.fmt(dz.p, 0)}%</b></div><div class="kpi"><span>기전 축 가중 예방점수</span><b>${C.fmt(sc, 0)}%</b></div></div>
      <p class="hint">${esc(condTxt)}</p>
      <div class="grid2"><div class="card">${C.paired(prof.map(r => ({ label: r.bm, sub: AI.AXES[r.axis], a: r.s, b: r.w, e: r.ci, p: r.p })), { xlab: "무처리 대조 대비 변화 (%) · 오차막대 = 95% 예측구간", title: "기전지표 예측" })}</div>
      <div class="card">${axes.length >= 3 ? C.radar(axes.map(a => AI.AXES[a]), ra, rb) : "<p>기전 축이 3개 미만이라 레이더를 생략합니다.</p>"}<p class="hint">기전 축별 평균 손상 크기(상대값). 초록 영역이 작을수록 보호효과가 큽니다.</p></div></div>
      <p class="hint">예방률 = (스트레스 단독 − 병용) ÷ 스트레스 단독 × 100. 질환 지수는 질환 관련 기전 축 손상을 가중 평균해 0~100으로 환산. 스트레스 영향이 10% 미만인 지표는 예방률을 계산하지 않습니다.</p>`;
  } else if (s.tab === "dose") {
    const pts = AI.dose(m, c, s.disease), f = AI.fit4pl(pts), curve = []; for (let i = 0; i <= 60; i++) { const x = o.cmax * 0.005 * 200 ** (i / 60); curve.push([x, f.top * x ** f.hill / (f.ec50 ** f.hill + x ** f.hill)]); }
    el.innerHTML = `<div class="kpis"><div class="kpi"><span>EC₅₀</span><b>${C.fmt(f.ec50, 3)} ${o.unit}</b></div><div class="kpi"><span>최대 예방점수 (외삽)</span><b>${C.fmt(f.top, 0)}%</b></div><div class="kpi"><span>권장 농도 (최대효과 80%)</span><b class="good">${C.fmt(f.rec80, 3)} ${o.unit}</b></div></div>
      <div class="card">${C.xy({ logx: true, series: [{ pts: curve, color: "var(--accent)", line: true, noDots: true, name: "4PL 적합" }, { pts: pts.filter(p => p[0] > 0), color: "var(--accent-2)", name: "AI 예측 예방점수" }], vlines: [{ x: f.ec50, label: "EC₅₀" }, { x: f.rec80, label: "권장", c: "var(--warn)" }], xlab: `농도 (${o.unit}, log)`, ylab: `${AI.DIS[s.disease].label} 예방점수 (%)` })}</div>
      <p class="hint">${esc(s.comp)} · ${esc(condTxt)}. 권장 농도는 최대 예방효과의 80%에 도달하는 최소 농도입니다.${f.rec80 > o.cmax ? " <b>시험 범위를 넘으므로 참고값입니다.</b>" : ""}</p>`;
  } else if (s.tab === "rank") {
    const rows = []; m.comps.forEach(cp => Object.keys(AI.TIM).forEach(t => { const p = AI.profile(m, { ...c, comp: cp, timing: t }, s.disease); rows.push({ cp, t, sc: AI.score(p, s.disease), dz: AI.disease(p, s.org, s.disease).p }); }));
    const best = m.comps.map(cp => rows.filter(r => r.cp === cp).sort((a, b) => b.sc - a.sc)[0]).sort((a, b) => b.sc - a.sc);
    el.innerHTML = `<div class="callout">이 조건의 1순위: <b>${esc(best[0].cp)}</b> (${AI.TIM[best[0].t][0]}) · 예방점수 ${C.fmt(best[0].sc, 0)}%</div>
      <div class="grid2"><div class="card">${C.hbars(best.map(b => ({ k: b.cp, v: b.sc })), { unit: "%", xlab: `${AI.DIS[s.disease].label} 예방점수 (물질별 최적 시점)`, left: 230 })}</div>
      <div class="card tblwrap"><table><thead><tr><th>후보물질</th>${Object.keys(AI.TIM).map(t => `<th class="n">${AI.TIM[t][0]}</th>`).join("")}<th class="n">질환 예방률(최적)</th></tr></thead><tbody>${m.comps.map(cp => `<tr><td>${esc(cp)}</td>${Object.keys(AI.TIM).map(t => `<td class="n">${C.fmt(rows.find(r => r.cp === cp && r.t === t).sc, 1)}</td>`).join("")}<td class="n">${C.fmt(best.find(b => b.cp === cp).dz, 0)}%</td></tr>`).join("")}</tbody></table></div></div>
      <p class="hint">현재 개체·스트레스·강도·기간·농도(${C.fmt(s.conc, 3)} ${o.unit}) 조건. 후보물질 A~E는 시연용 가상 물질이며, 데이터 허브에 실제 물질 데이터를 올리면 함께 비교됩니다.</p>`;
  } else if (s.tab === "perf") {
    const bms = Object.keys(m.models).sort(), sel = AI.impSel && m.models[AI.impSel] ? AI.impSel : bms[0];
    el.innerHTML = `<div class="card"><p><b>모델</b>: 지표별 가우시안 프로세스 회귀 (RBF 커널, 주변우도 최대화로 하이퍼파라미터 선택). 입력 = 개체 · 스트레스 종류×강도 · 노출 기간 · 후보물질×log 농도 · 처리 시점, 출력 = 지표 변화(%)와 95% 예측구간. 학습 ${C.fmt(AI.trainMs / 1000, 1)}초.</p>
      <button class="btn primary" id="ai-cv">5-fold 교차검증 실행</button><div id="ai-cvout">${AI.cvRes ? "" : '<p class="hint">버튼을 누르면 지표별 예측 정확도를 계산합니다 (10~30초).</p>'}</div></div>
      <div class="card"><h3>예측에 영향을 준 요인 (순열 중요도)</h3><label for="ai-imp">지표</label> <select id="ai-imp">${bms.map(b => `<option ${b === sel ? "selected" : ""}>${esc(b)}</option>`).join("")}</select><div id="ai-impout"></div></div>`;
    const showCV = () => { const { res, pairs } = AI.cvRes; const mr = S.mean(res.map(r => r.r2)), mc = S.mean(res.map(r => r.cover));
      el.querySelector("#ai-cvout").innerHTML = `<div class="kpis"><div class="kpi"><span>평균 R²</span><b>${C.fmt(mr, 3)}</b></div><div class="kpi"><span>95% 예측구간 포함률</span><b>${C.fmt(mc * 100, 0)}%</b></div><div class="kpi"><span>평균 R² (정답 대비)</span><b>${C.fmt(S.mean(res.filter(r => Number.isFinite(r.r2t)).map(r => r.r2t)), 3)}</b></div></div>
      <div class="tblwrap"><table><thead><tr><th>지표</th><th class="n">n</th><th class="n">R²</th><th class="n">RMSE</th><th class="n">95% 구간 포함률</th><th class="n">R² (정답 대비)</th></tr></thead><tbody>${res.map(r => `<tr><td>${esc(r.bm)}</td><td class="n">${r.n}</td><td class="n">${C.fmt(r.r2, 3)}</td><td class="n">${C.fmt(r.rmse, 1)}</td><td class="n">${C.fmt(r.cover * 100, 0)}%</td><td class="n">${C.fmt(r.r2t, 3)}</td></tr>`).join("")}</tbody></table></div>
      <div class="minis">${pairs.map(p => `<div><h4>${esc(p.bm)}</h4>${C.xy({ W: 240, H: 200, series: [{ pts: p.pts, color: "var(--accent)", r: 1.8, op: .5 }], diag: true, xlab: "관측", ylab: "예측" })}</div>`).join("")}</div>`; };
    if (AI.cvRes) showCV();
    el.querySelector("#ai-cv").onclick = () => { el.querySelector("#ai-cvout").innerHTML = '<p class="loading">교차검증 중…</p>'; setTimeout(() => { AI.cvRes = AI.cv(AI.data); showCV(); }, 30); };
    const imp = () => { AI.impSel = el.querySelector("#ai-imp").value; el.querySelector("#ai-impout").innerHTML = C.hbars(AI.importance(m, AI.impSel).map(d => ({ k: d.k, v: d.v, c: "var(--accent-2)" })), { unit: "%", xlab: "예측 기여도 (%)", left: 160 }); };
    el.querySelector("#ai-imp").onchange = imp; imp();
  } else if (s.tab === "hub") {
    const src = {}; AI.data.forEach(d => { const k = d.source + "|" + d.org; src[k] = (src[k] || 0) + 1; });
    el.innerHTML = `<div class="card"><p>실험값이나 <b>논문(PubMed/PMC)에서 추출한 값</b>을 표준 양식으로 올리면 시뮬레이션 데이터와 합쳐 모델을 다시 학습합니다. <code>pmid</code> 칸에 번호를 넣으면 '문헌', 비워 두면 '실데이터'로 구분됩니다.</p>
      <div class="row"><button class="btn" id="ai-tpl">입력 양식 내려받기 (CSV)</button><label class="btn primary" for="ai-up">실험·문헌 데이터 올리기</label><input id="ai-up" type="file" accept=".csv,.txt,.tsv" hidden>${AI.state.extra.length ? '<button class="btn" id="ai-clr">올린 데이터 제거</button>' : ""}</div><div id="ai-upmsg"></div>
      <details><summary>양식 설명</summary><ul><li><b>organism</b>: cell / mouse / human · <b>stress</b>: uv / ox / gly / allergen / barrier</li><li><b>level</b>: 0.25(약)·0.5(중)·0.75(강)·1.0(최강) · <b>days</b>: 노출 기간(일)</li><li><b>compound</b>: 물질명(스트레스 단독군은 <code>없음</code>) · <b>conc</b>: 세포 µM, 마우스·인체 도포 % · <b>timing</b>: pre / co / post</li><li><b>biomarker</b>: ${AI.BM.map(b => esc(b.name)).join(", ")}</li><li><b>value</b>: 무처리 대조 대비 % 변화 · <b>pmid</b>: 문헌 번호(선택)</li></ul></details></div>
      <div class="card"><h3>현재 학습 데이터 구성 (측정값 수)</h3><table><thead><tr><th>출처</th>${Object.keys(AI.ORG).map(k => `<th class="n">${AI.ORG[k].label}</th>`).join("")}</tr></thead><tbody>${["시뮬레이션", "실데이터", "문헌"].map(sr => `<tr><td>${sr}</td>${Object.keys(AI.ORG).map(k => `<td class="n">${(src[sr + "|" + k] || 0).toLocaleString()}</td>`).join("")}</tr>`).join("")}</tbody></table>
      <p class="hint">시뮬레이션 데이터는 스트레스 강도·기간의 포화형 손상, 개체 민감도, 후보물질 Hill 용량-반응과 처리 시점 효과를 반영한 기전 모델에서 생성한 검증용 데이터입니다. 후보물질 A~E는 시연용 가상 물질입니다.</p></div>
      ${AI.state.extra.length ? `<div class="card tblwrap"><h3>올린 데이터 (${AI.state.extra.length}건)</h3><table><thead><tr><th>출처</th><th>개체</th><th>스트레스</th><th>물질</th><th class="n">농도</th><th>지표</th><th class="n">값(%)</th><th>PMID</th></tr></thead><tbody>${AI.state.extra.slice(0, 200).map(d => `<tr><td>${d.source}</td><td>${d.org}</td><td>${d.stress}</td><td>${esc(d.comp)}</td><td class="n">${d.conc}</td><td>${esc(d.bm)}</td><td class="n">${C.fmt(d.value, 1)}</td><td>${esc(d.pmid || "")}</td></tr>`).join("")}</tbody></table></div>` : ""}`;
    el.querySelector("#ai-tpl").onclick = () => App.download("SkinAge-AI_입력양식.csv", S.toCSV(AI.TEMPLATE[0], AI.TEMPLATE.slice(1)), "text/csv");
    el.querySelector("#ai-up").onchange = ev => { const f = ev.target.files[0]; if (!f) return; f.text().then(t => { try { const { out, bad } = AI.parseReal(t); if (!out.length) throw new Error("읽을 수 있는 행이 없습니다."); AI.state.extra = out; AI.root.querySelector("#ai-body").innerHTML = '<p class="loading card">올린 데이터로 다시 학습하는 중…</p>'; setTimeout(() => { AI.retrain(); AI.render(); AI.root.querySelector("#ai-upmsg") && (AI.root.querySelector("#ai-upmsg").innerHTML = `<p class="ok">${out.length}건 반영 (실데이터 ${out.filter(d => d.source === "실데이터").length}, 문헌 ${out.filter(d => d.source === "문헌").length})${bad.length ? ` · 형식이 맞지 않아 제외한 행: ${bad.slice(0, 10).join(", ")}${bad.length > 10 ? " 외" : ""}` : ""}</p>`); }, 30); } catch (e) { el.querySelector("#ai-upmsg").innerHTML = `<p class="err">${esc(e.message)}</p>`; } }); };
    const clr = el.querySelector("#ai-clr"); if (clr) clr.onclick = () => { AI.state.extra = []; AI.retrain(); AI.render(); };
  } else if (s.tab === "rev") {
    el.innerHTML = `<div class="card"><p>후보물질 처리 전사체(DEG)가 <b>피부노화 유전자 패널의 노화 방향을 얼마나 되돌리는지</b> 계산합니다 (Connectivity Map 방식 역전 점수, −1~+1, 양수 = 노화 역전).</p>
      <div class="row"><button class="btn primary" id="rv-exo">실데이터 예시: 엑소좀 EV17 (마우스 색소세포)</button><label class="btn" for="rv-up">DEG 결과표 올리기</label><input id="rv-up" type="file" accept=".csv,.txt,.tsv" hidden></div><div id="rv-out"></div></div>`;
    const run = deg => { const r = Omics.reversal(deg), cats = {}; r.used.forEach(u => (cats[u.cat] = cats[u.cat] || []).push(u.rev));
      el.querySelector("#rv-out").innerHTML = `<div class="kpis"><div class="kpi"><span>피부노화 역전 점수</span><b class="${r.score > 0 ? "good" : ""}">${(r.score >= 0 ? "+" : "") + C.fmt(r.score, 3)}</b></div><div class="kpi"><span>사용된 패널 유전자</span><b>${r.used.length}개</b></div></div>
      ${C.hbars(Object.entries(cats).map(([k, v]) => ({ k, v: S.mean(v), c: S.mean(v) > 0 ? "var(--accent)" : "var(--down)" })).sort((a, b) => b.v - a.v), { dec: 2, xlab: "범주별 역전 점수 (+ 노화 역전 / − 노화 방향)", left: 170 })}
      <div class="tblwrap"><table><thead><tr><th>유전자</th><th>범주</th><th>노화 방향</th><th class="n">log2FC</th><th class="n">역전</th></tr></thead><tbody>${r.used.map(u => `<tr><td>${u.gene}</td><td>${u.cat}</td><td>${u.dir === "up" ? "▲" : "▼"}</td><td class="n">${C.fmt(u.lfc, 2)}</td><td class="n">${C.fmt(u.rev, 2)}</td></tr>`).join("")}</tbody></table></div>
      <p class="hint">마우스 유전자 기호는 대문자로 사람 기호에 대응했습니다. 역전 점수는 기전 가설을 세우기 위한 탐색 지표입니다.</p>`; };
    el.querySelector("#rv-exo").onclick = () => Omics.loadExoDEG().then(run);
    el.querySelector("#rv-up").onchange = ev => { const f = ev.target.files[0]; if (f) f.text().then(t => { try { run(Omics.standardizeDEG(t)); } catch (e) { el.querySelector("#rv-out").innerHTML = `<p class="err">${esc(e.message)}</p>`; } }); };
  }
};
