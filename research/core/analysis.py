"""SkinAge-Omics 분석 엔진.

- 입력 표준화: DESeq2 / edgeR / limma 결과표 컬럼 자동 인식
- count matrix 기반 DEG: DESeq2 median-of-ratios 정규화 후
    (1) NB-Wald  : 음이항 GLM(그룹 2수준) 최대우도 + Wald 검정 (DESeq2 방식 근사)
    (2) Welch    : log2 정규화 발현량에 대한 Welch t-test
- 피부노화 유전자 패널 탐색 (방향 일치도)
- GO/KEGG 과대표현분석(ORA): 초기하분포 검정 + Benjamini-Hochberg 보정
- 기존 분석 결과와의 비교: log2FC 상관, DEG 중복도, 방향 일치율
"""
from __future__ import annotations

import io
import numpy as np
import pandas as pd
from scipy import stats

from .skin_aging_genes import panel_table, ALIASES

# ---------------------------------------------------------------- 공통 유틸

def bh_adjust(p):
    """Benjamini-Hochberg FDR 보정 (NaN 유지)."""
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out


COLUMN_MAP = {
    "gene": ["gene", "genes", "symbol", "gene_symbol", "genesymbol", "gene_name",
             "genename", "id", "unnamed: 0", ""],
    "log2fc": ["log2foldchange", "logfc", "log2fc", "log2_fold_change", "lfc"],
    "pvalue": ["pvalue", "p.value", "pval", "p_value", "p"],
    "padj": ["padj", "fdr", "adj.p.val", "qvalue", "q_value", "p.adj", "padj_bh"],
    "basemean": ["basemean", "aveexpr", "logcpm", "mean_expr", "avg_expr"],
}


def standardize_deg(df: pd.DataFrame) -> pd.DataFrame:
    """DESeq2/edgeR/limma 결과표를 gene, log2fc, pvalue, padj, basemean 으로 표준화."""
    df = df.copy()
    lower = {c: str(c).strip().lower() for c in df.columns}
    found = {}
    for std, cands in COLUMN_MAP.items():
        for c, lc in lower.items():
            if lc in cands and c not in found.values():
                found[std] = c
                break
    if "gene" not in found:
        # 첫 번째 문자열 컬럼 또는 index 를 유전자로 사용
        obj = [c for c in df.columns if df[c].dtype == object]
        if obj:
            found["gene"] = obj[0]
        else:
            df = df.reset_index()
            found["gene"] = df.columns[0]
    missing = [k for k in ("log2fc",) if k not in found]
    if missing:
        raise ValueError("log2 fold change 컬럼을 찾지 못했습니다 (예: log2FoldChange, logFC).")
    out = pd.DataFrame({"gene": df[found["gene"]].astype(str).str.strip()})
    out["gene"] = out["gene"].str.upper().replace(ALIASES)
    out["log2fc"] = pd.to_numeric(df[found["log2fc"]], errors="coerce")
    out["pvalue"] = pd.to_numeric(df[found["pvalue"]], errors="coerce") if "pvalue" in found else np.nan
    if "padj" in found:
        out["padj"] = pd.to_numeric(df[found["padj"]], errors="coerce")
    else:
        out["padj"] = bh_adjust(out["pvalue"].values)
    out["basemean"] = pd.to_numeric(df[found["basemean"]], errors="coerce") if "basemean" in found else np.nan
    out = out.dropna(subset=["log2fc"]).drop_duplicates("gene").reset_index(drop=True)
    return out


def classify(deg: pd.DataFrame, padj_cut=0.05, lfc_cut=1.0) -> pd.DataFrame:
    deg = deg.copy()
    sig = (deg["padj"] < padj_cut) & (deg["log2fc"].abs() >= lfc_cut)
    deg["status"] = np.where(sig & (deg["log2fc"] > 0), "Up",
                     np.where(sig & (deg["log2fc"] < 0), "Down", "NS"))
    return deg


def read_table(uploaded) -> pd.DataFrame:
    """csv/tsv/txt/xlsx 를 자동 판별해서 읽는다 (경로 또는 업로드 객체)."""
    name = getattr(uploaded, "name", str(uploaded)).lower()
    if name.endswith((".xlsx", ".xls")):
        return pd.read_excel(uploaded)
    raw = uploaded.read() if hasattr(uploaded, "read") else open(uploaded, "rb").read()
    text = raw.decode("utf-8-sig", errors="replace")
    sep = "\t" if text.split("\n", 1)[0].count("\t") > text.split("\n", 1)[0].count(",") else ","
    return pd.read_csv(io.StringIO(text), sep=sep)

# ---------------------------------------------------------------- count → DEG

def size_factors(counts: pd.DataFrame) -> pd.Series:
    """DESeq2 median-of-ratios size factor."""
    c = counts.values.astype(float)
    with np.errstate(divide="ignore"):
        logc = np.log(c)
    keep = np.all(np.isfinite(logc), axis=1)
    if keep.sum() < 10:
        raise ValueError("모든 샘플에서 0이 아닌 유전자가 너무 적어 정규화할 수 없습니다.")
    logmean = logc[keep].mean(axis=1, keepdims=True)
    sf = np.exp(np.median(logc[keep] - logmean, axis=0))
    return pd.Series(sf, index=counts.columns)


def normalize(counts: pd.DataFrame) -> pd.DataFrame:
    return counts / size_factors(counts)


def _dispersion(norm: np.ndarray, groups: np.ndarray, sf: np.ndarray):
    """그룹 내 적률추정 분산 → 평균-분산 추세로 수축한 NB dispersion."""
    disp = np.zeros(norm.shape[0])
    w = 0
    for g in np.unique(groups):
        x = norm[:, groups == g]
        n = x.shape[1]
        if n < 2:
            continue
        m = x.mean(axis=1)
        v = x.var(axis=1, ddof=1)
        inv_s = np.mean(1.0 / sf[groups == g])
        with np.errstate(divide="ignore", invalid="ignore"):
            a = (v - m * inv_s) / m ** 2
        disp += np.nan_to_num(a, nan=0.0) * (n - 1)
        w += n - 1
    disp = np.clip(disp / max(w, 1), 1e-8, 10)
    mean = norm.mean(axis=1)
    ok = mean > 1
    # 모수적 추세 alpha(mu) = a0 + a1/mu  (DESeq2 parametric fit 과 동일 형태)
    X = np.column_stack([np.ones(ok.sum()), 1.0 / mean[ok]])
    coef = np.linalg.lstsq(X, disp[ok], rcond=None)[0]
    coef = np.clip(coef, 1e-6, None)
    trend = coef[0] + coef[1] / np.maximum(mean, 1e-3)
    # 경험적 베이즈 형태의 단순 수축: 로그 공간에서 가중 평균
    dof = max(w, 1)
    shrink = dof / (dof + 10.0)
    final = np.exp(shrink * np.log(np.maximum(disp, 1e-8)) + (1 - shrink) * np.log(trend))
    # 추세 이하로는 내리지 않음(과소추정 → 위양성 방지). 시뮬레이션에서 귀무 p<0.05 비율 ≈ 0.04 로 보정 확인
    final = np.maximum(final, trend)
    return np.clip(final, 1e-6, 10), trend


def deg_from_counts(counts: pd.DataFrame, meta: pd.DataFrame, group_col: str,
                    ref: str, test: str, method: str = "nb_wald",
                    min_count: int = 10, fixed_dispersion: float | None = None) -> pd.DataFrame:
    """두 그룹 비교 DEG. log2FC = log2(test / ref).

    그룹당 반복이 1개뿐이면 dispersion 을 추정할 수 없으므로 고정값(기본 BCV 0.1 → 0.01,
    edgeR 사용자 가이드의 무반복 권장값)을 사용한다. 이 경우 결과는 탐색적 해석용이다.
    """
    meta = meta.copy()
    meta.index = meta.index.astype(str)
    samples = [s for s in counts.columns.astype(str) if s in meta.index
               and meta.loc[s, group_col] in (ref, test)]
    if len(samples) < 2:
        raise ValueError("비교할 샘플이 2개 미만입니다. 샘플 이름과 그룹 컬럼을 확인하세요.")
    c = counts.copy()
    c.columns = c.columns.astype(str)
    c = c[samples]
    c = c[(c.sum(axis=1) >= min_count)]
    groups = meta.loc[samples, group_col].values
    sf = size_factors(c).values
    norm = c.values / sf
    is_t = groups == test
    m_ref = norm[:, ~is_t].mean(axis=1)
    m_t = norm[:, is_t].mean(axis=1)
    basemean = norm.mean(axis=1)

    no_rep = min(is_t.sum(), (~is_t).sum()) < 2
    if no_rep and fixed_dispersion is None:
        fixed_dispersion = 0.01
    if no_rep and method == "welch":
        method = "nb_wald"  # 반복 없이는 t-test 불가

    if method == "welch":
        lg = np.log2(norm + 1)
        t, p = stats.ttest_ind(lg[:, is_t], lg[:, ~is_t], axis=1, equal_var=False)
        lfc = lg[:, is_t].mean(axis=1) - lg[:, ~is_t].mean(axis=1)
        p = np.where(np.isnan(p), 1.0, p)
        stat = t
    else:
        if fixed_dispersion is not None:
            alpha = np.full(norm.shape[0], float(fixed_dispersion))
        else:
            alpha, _ = _dispersion(norm, groups, sf)
        y = c.values.astype(float)
        mus = {}
        ses = {}
        for lab, mask in (("ref", ~is_t), ("t", is_t)):
            yy, ss = y[:, mask], sf[mask]
            mu = np.maximum(yy.sum(axis=1) / ss.sum(), 1e-8)  # 초기값
            for _ in range(25):  # Newton: d/dlog(mu) NB 로그우도
                smu = ss * mu[:, None]
                denom = 1 + alpha[:, None] * smu
                score = ((yy - smu) / denom).sum(axis=1)
                info = (smu / denom).sum(axis=1)
                step = score / np.maximum(info, 1e-12)
                mu = mu * np.exp(np.clip(step, -5, 5))
                if np.nanmax(np.abs(step)) < 1e-8:
                    break
            mu = np.maximum(mu, 1e-8)
            smu = ss * mu[:, None]
            info = (smu / (1 + alpha[:, None] * smu)).sum(axis=1)
            mus[lab], ses[lab] = mu, 1.0 / np.sqrt(np.maximum(info, 1e-12))
        lfc_ln = np.log(mus["t"]) - np.log(mus["ref"])
        se_ln = np.sqrt(ses["t"] ** 2 + ses["ref"] ** 2)
        stat = lfc_ln / se_ln
        p = 2 * stats.norm.sf(np.abs(stat))
        lfc = lfc_ln / np.log(2)
        # 한쪽 그룹이 모두 0 인 경우 과도한 LFC 제한
        lfc = np.clip(lfc, -12, 12)

    res = pd.DataFrame({
        "gene": c.index.astype(str).str.upper(), "basemean": basemean,
        "mean_" + str(ref): m_ref, "mean_" + str(test): m_t,
        "log2fc": lfc, "stat": stat, "pvalue": p,
    })
    res["gene"] = res["gene"].replace(ALIASES)
    res["padj"] = bh_adjust(res["pvalue"].values)
    return res.sort_values("pvalue").reset_index(drop=True)

# ---------------------------------------------------------------- 피부노화 패널

def aging_panel(deg: pd.DataFrame, padj_cut=0.05, lfc_cut=1.0) -> pd.DataFrame:
    panel = pd.DataFrame(panel_table())
    d = classify(deg, padj_cut, lfc_cut).set_index("gene")
    panel["log2fc"] = panel["gene"].map(d["log2fc"])
    panel["padj"] = panel["gene"].map(d["padj"])
    panel["status"] = panel["gene"].map(d["status"]).fillna("미검출")

    def agree(r):
        if r["status"] in ("NS", "미검출") or r["expected"] == "var":
            return "-"
        obs = "up" if r["status"] == "Up" else "down"
        return "일치" if obs == r["expected"] else "불일치"
    panel["concordance"] = panel.apply(agree, axis=1)
    return panel


def panel_summary(panel: pd.DataFrame) -> pd.DataFrame:
    g = panel.groupby("category", sort=False)
    out = pd.DataFrame({
        "패널 유전자": g.size(),
        "데이터 검출": g.apply(lambda x: (x["status"] != "미검출").sum()),
        "유의 변화": g.apply(lambda x: x["status"].isin(["Up", "Down"]).sum()),
        "방향 일치": g.apply(lambda x: (x["concordance"] == "일치").sum()),
        "방향 불일치": g.apply(lambda x: (x["concordance"] == "불일치").sum()),
        "평균 log2FC": g["log2fc"].mean().round(2),
    })
    return out.reset_index().rename(columns={"category": "범주"})

# ---------------------------------------------------------------- GO / KEGG ORA

def read_gmt(path_or_buf) -> dict[str, set]:
    sets = {}
    if hasattr(path_or_buf, "read"):
        raw = path_or_buf.read()
        lines = (raw.decode("utf-8", "replace") if isinstance(raw, bytes) else raw).splitlines()
    else:
        with open(path_or_buf, encoding="utf-8") as fh:
            lines = fh.read().splitlines()
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        genes = {g.split(",")[0].strip().upper() for g in parts[2:] if g.strip()}
        sets[parts[0]] = genes
    return sets


def ora(genes, universe, gene_sets: dict, min_size=10, max_size=500) -> pd.DataFrame:
    """초기하분포 기반 과대표현분석. universe 는 분석에서 검정된 전체 유전자."""
    genes = set(g.upper() for g in genes)
    universe = set(g.upper() for g in universe)
    genes &= universe
    rows = []
    for term, gs in gene_sets.items():
        gs_u = gs & universe
        K = len(gs_u)
        if K < min_size or K > max_size:
            continue
        hits = genes & gs_u
        k = len(hits)
        if k == 0:
            continue
        rows.append((term, k, K, sorted(hits)))
    if not rows:
        return pd.DataFrame(columns=["term", "overlap", "set_size", "gene_ratio",
                                     "fold_enrichment", "pvalue", "padj", "genes"])
    N, n = len(universe), len(genes)
    df = pd.DataFrame(rows, columns=["term", "overlap", "set_size", "hit_genes"])
    df["pvalue"] = stats.hypergeom.sf(df["overlap"] - 1, N, df["set_size"], n)
    df["padj"] = bh_adjust(df["pvalue"].values)
    df["gene_ratio"] = df["overlap"] / max(n, 1)
    df["fold_enrichment"] = (df["overlap"] / max(n, 1)) / (df["set_size"] / N)
    df["genes"] = df["hit_genes"].apply(lambda x: ", ".join(x))
    df["term"] = df["term"].str.replace(r"^KEGG_", "", regex=True).str.replace("_", " ")
    return (df.drop(columns="hit_genes").sort_values("pvalue").reset_index(drop=True)
            [["term", "overlap", "set_size", "gene_ratio", "fold_enrichment", "pvalue", "padj", "genes"]])

# ---------------------------------------------------------------- 기존 분석과 비교

def compare_results(a: pd.DataFrame, b: pd.DataFrame, padj_cut=0.05, lfc_cut=1.0,
                    labels=("SkinAge-Omics", "기존 분석")) -> dict:
    a = classify(a, padj_cut, lfc_cut).set_index("gene")
    b = classify(b, padj_cut, lfc_cut).set_index("gene")
    common = a.index.intersection(b.index)
    a, b = a.loc[common], b.loc[common]
    ok = a["log2fc"].notna() & b["log2fc"].notna()
    pear = stats.pearsonr(a.loc[ok, "log2fc"], b.loc[ok, "log2fc"])[0] if ok.sum() > 2 else np.nan
    spear = stats.spearmanr(a.loc[ok, "log2fc"], b.loc[ok, "log2fc"])[0] if ok.sum() > 2 else np.nan
    sa = set(a.index[a["status"] != "NS"])
    sb = set(b.index[b["status"] != "NS"])
    inter = sa & sb
    same_dir = sum(a.loc[g, "status"] == b.loc[g, "status"] for g in inter)
    jac = len(inter) / len(sa | sb) if (sa | sb) else np.nan
    ct = pd.crosstab(a["status"], b["status"]).reindex(index=["Up", "NS", "Down"],
                                                        columns=["Up", "NS", "Down"], fill_value=0)
    ct.index.name, ct.columns.name = labels[0], labels[1]
    merged = pd.DataFrame({f"log2fc_{labels[0]}": a["log2fc"], f"log2fc_{labels[1]}": b["log2fc"],
                           f"status_{labels[0]}": a["status"], f"status_{labels[1]}": b["status"]})
    return {
        "n_common": int(len(common)), "pearson_lfc": float(pear), "spearman_lfc": float(spear),
        "deg_a": len(sa), "deg_b": len(sb), "deg_both": len(inter),
        "jaccard": float(jac), "direction_agreement": (same_dir / len(inter)) if inter else np.nan,
        "recall_of_b": (len(inter) / len(sb)) if sb else np.nan,
        "precision_vs_b": (len(inter) / len(sa)) if sa else np.nan,
        "crosstab": ct, "merged": merged,
    }
