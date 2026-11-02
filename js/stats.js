/* SkinAge-AI 통계·수치 계산 (브라우저용, 의존성 없음) */
"use strict";
const S = {};

/* 재현 가능한 난수 (mulberry32 + Box–Muller) */
S.rng = function (seed) {
  let a = seed >>> 0;
  const u = () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
  return {
    u, norm: (m = 0, s = 1) => { let x = 0, y = 0; while (x === 0) x = u(); y = u(); return m + s * Math.sqrt(-2 * Math.log(x)) * Math.cos(2 * Math.PI * y); },
    choice: (arr, p) => { if (!p) return arr[Math.floor(u() * arr.length)]; let r = u(), c = 0; for (let i = 0; i < arr.length; i++) { c += p[i]; if (r <= c) return arr[i]; } return arr[arr.length - 1]; },
    perm: n => { const a = [...Array(n).keys()]; for (let i = n - 1; i > 0; i--) { const j = Math.floor(u() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; } return a; }
  };
};

S.mean = a => a.reduce((s, x) => s + x, 0) / a.length;
S.sd = a => { const m = S.mean(a); return Math.sqrt(a.reduce((s, x) => s + (x - m) ** 2, 0) / Math.max(a.length - 1, 1)); };
S.median = a => { const b = [...a].sort((x, y) => x - y), n = b.length; return n % 2 ? b[(n - 1) / 2] : (b[n / 2 - 1] + b[n / 2]) / 2; };
S.pearson = (x, y) => { const mx = S.mean(x), my = S.mean(y); let a = 0, b = 0, c = 0; for (let i = 0; i < x.length; i++) { a += (x[i] - mx) * (y[i] - my); b += (x[i] - mx) ** 2; c += (y[i] - my) ** 2; } return a / Math.sqrt(b * c); };
S.rank = a => { const idx = a.map((v, i) => [v, i]).sort((p, q) => p[0] - q[0]); const r = new Array(a.length); let i = 0; while (i < idx.length) { let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++; for (let k = i; k <= j; k++) r[idx[k][1]] = (i + j) / 2 + 1; i = j + 1; } return r; };
S.spearman = (x, y) => S.pearson(S.rank(x), S.rank(y));

/* Benjamini–Hochberg */
S.bh = p => {
  const n = p.length, idx = p.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0]), out = new Array(n);
  let prev = 1;
  for (let k = n - 1; k >= 0; k--) { const [v, i] = idx[k]; prev = Math.min(prev, v * n / (k + 1)); out[i] = Math.min(prev, 1); }
  return out;
};

/* 정규분포 양측 p */
S.normSf = z => { const t = 1 / (1 + 0.2316419 * Math.abs(z)); const d = 0.3989422804014327 * Math.exp(-z * z / 2); const p = d * t * (0.31938153 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))); return z > 0 ? p : 1 - p; };
S.p2 = z => Math.min(1, 2 * S.normSf(Math.abs(z)));

/* log-gamma (Lanczos) · 초기하분포 상측확률 */
S.lgamma = x => { const g = 7, c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61503916999185, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]; if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - S.lgamma(1 - x); x -= 1; let a = c[0]; const t = x + g + 0.5; for (let i = 1; i < 9; i++) a += c[i] / (x + i); return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a); };
S.lchoose = (n, k) => S.lgamma(n + 1) - S.lgamma(k + 1) - S.lgamma(n - k + 1);
S.hyperSf = (k, N, K, n) => { /* P(X >= k) */
  let s = 0; const hi = Math.min(K, n), denom = S.lchoose(N, n);
  for (let i = k; i <= hi; i++) s += Math.exp(S.lchoose(K, i) + S.lchoose(N - K, n - i) - denom);
  return Math.min(1, s);
};

/* Welch t-test p (t 분포 근사: 정규 + 소표본 보정) */
S.tSf2 = (t, df) => { /* 양측: 정규 근사에 Cornish-Fisher 보정 */
  const x = df / (df + t * t); return S.ibeta(df / 2, 0.5, x);
};
S.ibeta = (a, b, x) => { /* 정규화 불완전 베타 (연분수) */
  if (x <= 0) return 0; if (x >= 1) return 1;
  const bt = Math.exp(S.lgamma(a + b) - S.lgamma(a) - S.lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  const cf = (a, b, x) => { let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap; if (Math.abs(d) < 1e-30) d = 1e-30; d = 1 / d; let h = d; for (let m = 1; m <= 200; m++) { const m2 = 2 * m; let aa = m * (b - m) * x / ((qam + m2) * (a + m2)); d = 1 + aa * d; if (Math.abs(d) < 1e-30) d = 1e-30; c = 1 + aa / c; if (Math.abs(c) < 1e-30) c = 1e-30; d = 1 / d; h *= d * c; aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2)); d = 1 + aa * d; if (Math.abs(d) < 1e-30) d = 1e-30; c = 1 + aa / c; if (Math.abs(c) < 1e-30) c = 1e-30; d = 1 / d; const del = d * c; h *= del; if (Math.abs(del - 1) < 3e-12) break; } return h; };
  return x < (a + 1) / (a + b + 2) ? bt * cf(a, b, x) / a : 1 - bt * cf(b, a, 1 - x) / b;
};

/* 선형대수: 촐레스키 */
S.cholesky = (K, n) => { /* K: Float64Array n*n (row-major) → L (in place lower) */
  const L = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let s = K[i * n + j];
      for (let k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
      if (i === j) { if (s <= 0) return null; L[i * n + i] = Math.sqrt(s); }
      else L[i * n + j] = s / L[j * n + j];
    }
  }
  return L;
};
S.solveL = (L, n, b) => { const x = new Float64Array(n); for (let i = 0; i < n; i++) { let s = b[i]; for (let k = 0; k < i; k++) s -= L[i * n + k] * x[k]; x[i] = s / L[i * n + i]; } return x; };
S.solveLT = (L, n, b) => { const x = new Float64Array(n); for (let i = n - 1; i >= 0; i--) { let s = b[i]; for (let k = i + 1; k < n; k++) s -= L[k * n + i] * x[k]; x[i] = s / L[i * n + i]; } return x; };

/* 가우시안 프로세스 회귀 (RBF, 표준화 타깃) */
S.GP = class {
  constructor(ls, noise) { this.ls = ls; this.noise = noise; }
  k(a, b) { let d = 0; for (let i = 0; i < a.length; i++) { const t = (a[i] - b[i]) / this.ls; d += t * t; } return Math.exp(-0.5 * d); }
  fit(X, y) {
    const n = X.length; this.X = X;
    this.ym = S.mean(y); this.ys = S.sd(y) || 1;
    const yt = y.map(v => (v - this.ym) / this.ys);
    const K = new Float64Array(n * n);
    for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) { const v = this.k(X[i], X[j]); K[i * n + j] = v; K[j * n + i] = v; }
    for (let i = 0; i < n; i++) K[i * n + i] += this.noise ** 2 + 1e-8;
    const L = S.cholesky(K, n); if (!L) return null;
    this.L = L; this.n = n;
    this.alpha = S.solveLT(L, n, S.solveL(L, n, yt));
    let ld = 0; for (let i = 0; i < n; i++) ld += Math.log(L[i * n + i]);
    let q = 0; for (let i = 0; i < n; i++) q += yt[i] * this.alpha[i];
    this.lml = -0.5 * q - ld - 0.5 * n * Math.log(2 * Math.PI);
    return this;
  }
  predict(x) {
    const n = this.n, ks = new Float64Array(n);
    let mu = 0; for (let i = 0; i < n; i++) { ks[i] = this.k(x, this.X[i]); mu += ks[i] * this.alpha[i]; }
    const v = S.solveL(this.L, n, ks); let vv = 0; for (let i = 0; i < n; i++) vv += v[i] * v[i];
    const varr = Math.max(1 - vv, 1e-12) + this.noise ** 2;
    return [mu * this.ys + this.ym, Math.sqrt(varr) * this.ys];
  }
};
S.fitGP = (X, y, lsGrid = [0.8, 1.4, 2.2], nzGrid = [0.1, 0.25]) => {
  let best = null;
  for (const ls of lsGrid) for (const nz of nzGrid) { const g = new S.GP(ls, nz).fit(X, y); if (g && (!best || g.lml > best.lml)) best = g; }
  return best;
};

/* CSV 파서 (따옴표 지원, 콤마/탭 자동) */
S.parseCSV = text => {
  text = text.replace(/^﻿/, "");
  const first = text.split(/\r?\n/, 1)[0];
  const sep = (first.match(/\t/g) || []).length > (first.match(/,/g) || []).length ? "\t" : ",";
  const rows = []; let row = [], cur = "", q = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (q) { if (ch === '"') { if (text[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += ch; }
    else if (ch === '"') q = true;
    else if (ch === sep) { row.push(cur); cur = ""; }
    else if (ch === "\n" || ch === "\r") { if (ch === "\r" && text[i + 1] === "\n") i++; row.push(cur); rows.push(row); row = []; cur = ""; }
    else cur += ch;
  }
  if (cur.length || row.length) { row.push(cur); rows.push(row); }
  const header = rows.shift().map(h => h.trim());
  return { header, rows: rows.filter(r => r.length > 1 || (r[0] || "").trim() !== "") };
};
S.toCSV = (header, rows) => "﻿" + [header, ...rows].map(r => r.map(v => { const s = v == null ? "" : String(v); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; }).join(",")).join("\n");
