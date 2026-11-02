/* SVG 차트 (테마 토큰 사용) */
"use strict";
const C = {};
const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
C.fmt = (x, d = 1) => Number.isFinite(x) ? (Math.abs(x) >= 1000 ? x.toLocaleString("ko-KR", { maximumFractionDigits: 0 }) : (+x.toFixed(d)).toLocaleString("ko-KR", { maximumFractionDigits: d })) : "–";
C.nice = v => { if (!(v > 0)) return 1; const p = 10 ** Math.floor(Math.log10(v)), m = v / p; return (m <= 1 ? 1 : m <= 2 ? 2 : m <= 2.5 ? 2.5 : m <= 5 ? 5 : 10) * p; };

/* 가로 막대 2계열 (스트레스 단독 vs 병용) + 오차막대 + 우측 라벨 */
C.paired = (rows, o = {}) => {
  const W = o.W || 720, rowH = 34, P = { l: 230, r: 86, t: 30, b: 40 }, H = P.t + P.b + rows.length * rowH;
  const vals = rows.flatMap(r => [r.a, r.b + (r.e || 0), r.b - (r.e || 0), 0]);
  let lo = Math.min(...vals), hi = Math.max(...vals); const span = C.nice((hi - lo) / 4);
  lo = Math.floor(lo / span) * span; hi = Math.ceil(hi / span) * span;
  const iw = W - P.l - P.r, X = v => P.l + (v - lo) / (hi - lo) * iw;
  let s = `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.title || "")}">`;
  for (let v = lo; v <= hi + 1e-9; v += span) s += `<line x1="${X(v)}" x2="${X(v)}" y1="${P.t}" y2="${H - P.b}" stroke="var(--grid)"/><text x="${X(v)}" y="${H - P.b + 15}" text-anchor="middle">${C.fmt(v, 0)}</text>`;
  s += `<line x1="${X(0)}" x2="${X(0)}" y1="${P.t}" y2="${H - P.b}" stroke="var(--sub)"/>`;
  rows.forEach((r, i) => {
    const y = P.t + i * rowH;
    s += `<text x="${P.l - 8}" y="${y + 14}" text-anchor="end" class="lab">${esc(r.label)}</text><text x="${P.l - 8}" y="${y + 27}" text-anchor="end" class="sm">${esc(r.sub || "")}</text>`;
    const bar = (v, yy, cl) => `<rect x="${Math.min(X(0), X(v))}" y="${yy}" width="${Math.abs(X(v) - X(0))}" height="11" rx="2" class="${cl}"/>`;
    s += bar(r.a, y + 4, "b-base") + bar(r.b, y + 17, "b-main");
    if (r.e) s += `<line x1="${X(r.b - r.e)}" x2="${X(r.b + r.e)}" y1="${y + 22.5}" y2="${y + 22.5}" stroke="var(--ink)" stroke-width="1"/>`;
    if (Number.isFinite(r.p)) s += `<text x="${W - P.r + 8}" y="${y + 20}" class="pv">예방 ${C.fmt(r.p, 0)}%</text>`;
  });
  s += `<text x="${P.l + iw / 2}" y="${H - 6}" text-anchor="middle">${esc(o.xlab || "")}</text>`;
  s += `<rect x="${P.l}" y="8" width="11" height="11" class="b-base"/><text x="${P.l + 16}" y="18">${esc(o.la || "스트레스 단독")}</text><rect x="${P.l + 120}" y="8" width="11" height="11" class="b-main"/><text x="${P.l + 136}" y="18">${esc(o.lb || "스트레스 + 후보물질")}</text>`;
  return s + "</svg>";
};

/* 레이더 */
C.radar = (labels, a, b, o = {}) => {
  const W = 420, H = 380, cx = 210, cy = 195, R = 140, n = labels.length, mx = Math.max(...a, ...b, 1);
  const pt = (i, v) => { const t = -Math.PI / 2 + i * 2 * Math.PI / n; return [cx + Math.cos(t) * R * v / mx, cy + Math.sin(t) * R * v / mx]; };
  let s = `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img" aria-label="기전 축 레이더">`;
  for (const f of [0.25, 0.5, 0.75, 1]) s += `<polygon points="${labels.map((_, i) => pt(i, mx * f).join(",")).join(" ")}" fill="none" stroke="var(--grid)"/>`;
  labels.forEach((l, i) => { const [x, y] = pt(i, mx * 1.13); s += `<line x1="${cx}" y1="${cy}" x2="${pt(i, mx)[0]}" y2="${pt(i, mx)[1]}" stroke="var(--grid)"/><text x="${x}" y="${y + 4}" text-anchor="middle" class="sm">${esc(l)}</text>`; });
  s += `<polygon points="${a.map((v, i) => pt(i, v).join(",")).join(" ")}" class="r-base"/><polygon points="${b.map((v, i) => pt(i, v).join(",")).join(" ")}" class="r-main"/>`;
  s += `<rect x="12" y="10" width="11" height="11" class="b-base"/><text x="28" y="20">스트레스 단독</text><rect x="12" y="28" width="11" height="11" class="b-main"/><text x="28" y="38">+ 후보물질</text>`;
  return s + "</svg>";
};

/* 산점/선 (log x 옵션, 곡선, 세로선) */
C.xy = (o) => {
  const W = o.W || 640, H = o.H || 300, P = { l: 56, r: 16, t: 16, b: 44 }, iw = W - P.l - P.r, ih = H - P.t - P.b;
  const tx = v => o.logx ? Math.log10(v) : v;
  const allx = [].concat(...o.series.map(s => s.pts.map(p => tx(p[0])))), ally = [].concat(...o.series.map(s => s.pts.map(p => p[1])));
  let x0 = o.x0 ?? Math.min(...allx), x1 = o.x1 ?? Math.max(...allx), y0 = o.y0 ?? Math.min(0, ...ally), y1 = o.y1 ?? C.nice(Math.max(...ally) * 1.08);
  if (o.logx) { x0 = Math.floor(x0); x1 = Math.ceil(x1); }
  const X = v => P.l + (tx(v) - x0) / (x1 - x0 || 1) * iw, Y = v => P.t + ih - (v - y0) / (y1 - y0 || 1) * ih;
  let s = `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.ylab || "")}">`;
  for (let i = 0; i <= 4; i++) { const v = y0 + (y1 - y0) * i / 4; s += `<line x1="${P.l}" x2="${W - P.r}" y1="${Y(v)}" y2="${Y(v)}" stroke="var(--grid)"/><text x="${P.l - 6}" y="${Y(v) + 4}" text-anchor="end">${C.fmt(v, Math.abs(y1 - y0) < 5 ? 2 : 0)}</text>`; }
  if (o.logx) for (let e = x0; e <= x1; e++) s += `<text x="${X(10 ** e)}" y="${H - P.b + 16}" text-anchor="middle">${C.fmt(10 ** e, 3)}</text>`;
  else for (let i = 0; i <= 5; i++) { const v = x0 + (x1 - x0) * i / 5; s += `<text x="${X(v)}" y="${H - P.b + 16}" text-anchor="middle">${C.fmt(v, 1)}</text>`; }
  s += `<line x1="${P.l}" x2="${W - P.r}" y1="${P.t + ih}" y2="${P.t + ih}" stroke="var(--sub)"/>`;
  (o.vlines || []).forEach(v => { if (tx(v.x) >= x0 && tx(v.x) <= x1) s += `<line x1="${X(v.x)}" x2="${X(v.x)}" y1="${P.t}" y2="${P.t + ih}" stroke="${v.c || "var(--sub)"}" stroke-dasharray="4 4"/><text x="${X(v.x) + 4}" y="${P.t + 12}" class="sm" style="fill:${v.c || "var(--sub)"}">${esc(v.label || "")}</text>`; });
  (o.hlines || []).forEach(v => s += `<line x1="${P.l}" x2="${W - P.r}" y1="${Y(v)}" y2="${Y(v)}" stroke="var(--sub)" stroke-dasharray="4 4"/>`);
  o.series.forEach(se => {
    if (se.line) s += `<path d="${se.pts.map((p, i) => (i ? "L" : "M") + X(p[0]).toFixed(1) + " " + Y(p[1]).toFixed(1)).join("")}" fill="none" stroke="${se.color}" stroke-width="2"/>`;
    if (!se.noDots) se.pts.forEach(p => s += `<circle cx="${X(p[0]).toFixed(1)}" cy="${Y(p[1]).toFixed(1)}" r="${se.r || 3.5}" fill="${se.color}" fill-opacity="${se.op || 0.9}"><title>${esc(p[2] || "")}</title></circle>`);
  });
  if (o.diag) s += `<line x1="${X(Math.max(x0, y0))}" y1="${Y(Math.max(x0, y0))}" x2="${X(Math.min(x1, y1))}" y2="${Y(Math.min(x1, y1))}" stroke="var(--sub)" stroke-dasharray="4 4"/>`;
  s += `<text x="${P.l + iw / 2}" y="${H - 6}" text-anchor="middle">${esc(o.xlab || "")}</text><text transform="translate(13 ${P.t + ih / 2}) rotate(-90)" text-anchor="middle">${esc(o.ylab || "")}</text>`;
  const leg = o.series.filter(s => s.name);
  leg.forEach((se, i) => s += `<circle cx="${P.l + 12 + i * 170}" cy="${P.t + 6}" r="4.5" fill="${se.color}"/><text x="${P.l + 22 + i * 170}" y="${P.t + 10}" style="fill:var(--ink)">${esc(se.name)}</text>`);
  return s + "</svg>";
};

/* 가로 막대 (값 하나) */
C.hbars = (items, o = {}) => {
  const W = o.W || 600, rowH = 26, P = { l: o.left || 220, r: 60, t: 10, b: 32 }, H = P.t + P.b + items.length * rowH;
  const vals = items.map(d => d.v), lo = Math.min(0, ...vals), hi = Math.max(0, ...vals) || 1, iw = W - P.l - P.r, X = v => P.l + (v - lo) / (hi - lo) * iw;
  let s = `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.title || "막대")}">`;
  s += `<line x1="${X(0)}" x2="${X(0)}" y1="${P.t}" y2="${H - P.b}" stroke="var(--sub)"/>`;
  items.forEach((d, i) => { const y = P.t + i * rowH; s += `<text x="${P.l - 8}" y="${y + 16}" text-anchor="end" class="lab">${esc(d.k)}</text><rect x="${Math.min(X(0), X(d.v))}" y="${y + 5}" width="${Math.abs(X(d.v) - X(0))}" height="15" rx="3" fill="${d.c || "var(--accent)"}"/><text x="${d.v >= 0 ? X(d.v) + 5 : X(d.v) - 5}" y="${y + 17}" text-anchor="${d.v >= 0 ? "start" : "end"}" style="fill:var(--ink)">${C.fmt(d.v, o.dec ?? 0)}${o.unit || ""}</text>`; });
  s += `<text x="${P.l + iw / 2}" y="${H - 8}" text-anchor="middle">${esc(o.xlab || "")}</text>`;
  return s + "</svg>";
};

/* 화산 그림 (NS 점은 표본추출) */
C.volcano = (deg, o) => {
  const W = 640, H = 420, P = { l: 52, r: 14, t: 14, b: 42 }, iw = W - P.l - P.r, ih = H - P.t - P.b;
  const pts = deg.filter(d => Number.isFinite(d.lfc) && Number.isFinite(d.padj));
  const ys = pts.map(d => -Math.log10(Math.max(d.padj, 1e-300)));
  const absl = pts.map(d => Math.abs(d.lfc)).sort((a, b) => a - b); const xl = Math.max(o.lfc * 2, absl[Math.floor(absl.length * 0.999)] * 1.1 || 4);
  const y1 = C.nice(Math.min(Math.max(...ys), 60) * 1.05);
  const X = v => P.l + (Math.max(-xl, Math.min(xl, v)) + xl) / (2 * xl) * iw, Y = v => P.t + ih - Math.min(v, y1) / y1 * ih;
  let s = `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img" aria-label="Volcano plot">`;
  for (let i = 0; i <= 4; i++) { const v = y1 * i / 4; s += `<line x1="${P.l}" x2="${W - P.r}" y1="${Y(v)}" y2="${Y(v)}" stroke="var(--grid)"/><text x="${P.l - 6}" y="${Y(v) + 4}" text-anchor="end">${C.fmt(v, 0)}</text>`; }
  for (let i = -2; i <= 2; i++) { const v = xl * i / 2; s += `<text x="${X(v)}" y="${H - P.b + 16}" text-anchor="middle">${C.fmt(v, 1)}</text>`; }
  let k = 0;
  pts.forEach((d, i) => { if (d.status === "NS" && (k++ % 3)) return; const cl = d.status === "Up" ? "var(--up)" : d.status === "Down" ? "var(--down)" : "var(--ns)"; s += `<circle cx="${X(d.lfc).toFixed(1)}" cy="${Y(ys[i]).toFixed(1)}" r="${d.status === "NS" ? 1.6 : 2.4}" fill="${cl}" fill-opacity="${d.status === "NS" ? .5 : .85}"><title>${esc(d.gene)} log2FC ${C.fmt(d.lfc, 2)} FDR ${d.padj.toExponential(1)}</title></circle>`; });
  s += `<line x1="${P.l}" x2="${W - P.r}" y1="${Y(-Math.log10(o.padj))}" y2="${Y(-Math.log10(o.padj))}" stroke="var(--sub)" stroke-dasharray="4 4"/>`;
  [-o.lfc, o.lfc].forEach(v => s += `<line x1="${X(v)}" x2="${X(v)}" y1="${P.t}" y2="${P.t + ih}" stroke="var(--sub)" stroke-dasharray="4 4"/>`);
  const lab = pts.filter(d => d.status !== "NS").sort((a, b) => a.padj - b.padj).slice(0, 12).concat(pts.filter(d => (o.highlight || []).includes(d.gene)));
  const seen = new Set(); lab.forEach(d => { if (seen.has(d.gene)) return; seen.add(d.gene); s += `<text x="${X(d.lfc) + 4}" y="${Y(-Math.log10(Math.max(d.padj, 1e-300))) - 3}" class="sm" style="fill:var(--ink)">${esc(d.gene)}</text>`; });
  s += `<text x="${P.l + iw / 2}" y="${H - 6}" text-anchor="middle">log2 Fold Change</text><text transform="translate(13 ${P.t + ih / 2}) rotate(-90)" text-anchor="middle">−log10(FDR)</text>`;
  return s + "</svg>";
};
