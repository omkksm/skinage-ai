"""시각화 함수 (앱과 CLI 공용, matplotlib 기반)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.cluster.hierarchy import linkage, leaves_list

_KO_FONTS = ["Malgun Gothic", "AppleGothic", "NanumGothic", "Noto Sans CJK KR", "Noto Sans KR"]
# 동봉한 한글 폰트(브라우저 실행 환경처럼 시스템 폰트가 없을 때 사용)
from pathlib import Path as _P
_bundled = _P(__file__).resolve().parents[1] / "fonts" / "NotoSansKR-subset.otf"
if _bundled.exists():
    try:
        font_manager.fontManager.addfont(str(_bundled))
    except Exception:
        pass
_avail = {f.name for f in font_manager.fontManager.ttflist}
plt.rcParams["font.family"] = [f for f in _KO_FONTS if f in _avail] + ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

UP, DOWN, NS = "#C8463D", "#2F6DB5", "#BDBDBD"
COLORS = {"Up": UP, "Down": DOWN, "NS": NS}


def volcano(deg, padj_cut=0.05, lfc_cut=1.0, label_genes=None, top_n=12, title=None):
    d = deg.copy()
    d["nlp"] = -np.log10(d["padj"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(7, 5.6))
    for s in ["NS", "Down", "Up"]:
        x = d[d["status"] == s]
        ax.scatter(x["log2fc"], x["nlp"], s=7 if s == "NS" else 11, c=COLORS[s],
                   alpha=0.5 if s == "NS" else 0.8, lw=0,
                   label=f"{s} ({len(x):,})")
    ax.axhline(-np.log10(padj_cut), ls="--", c="#777", lw=0.8)
    for v in (-lfc_cut, lfc_cut):
        ax.axvline(v, ls="--", c="#777", lw=0.8)
    lab = d[d["status"] != "NS"].nsmallest(top_n, "padj")
    if label_genes:
        lab = pd.concat([lab, d[d["gene"].isin(label_genes) & (d["status"] != "NS")]]).drop_duplicates("gene")
    for _, r in lab.iterrows():
        ax.annotate(r["gene"], (r["log2fc"], r["nlp"]), fontsize=7.5, xytext=(3, 2),
                    textcoords="offset points")
    lim = np.nanpercentile(np.abs(d["log2fc"]), 99.95) * 1.15
    ax.set_xlim(-max(lim, lfc_cut * 2), max(lim, lfc_cut * 2))
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(adjusted p)")
    ax.set_title(title or "Volcano plot")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    return fig


def ma_plot(deg, title=None):
    d = deg.dropna(subset=["basemean"])
    fig, ax = plt.subplots(figsize=(7, 4.6))
    for s in ["NS", "Down", "Up"]:
        x = d[d["status"] == s]
        ax.scatter(np.log10(x["basemean"].clip(lower=0.1)), x["log2fc"], s=6, c=COLORS[s],
                   alpha=0.5 if s == "NS" else 0.8, lw=0, label=s)
    ax.axhline(0, c="#555", lw=0.8)
    ax.set_xlabel("log10 평균 발현량 (baseMean)")
    ax.set_ylabel("log2 Fold Change")
    ax.set_title(title or "MA plot")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig


def heatmap(norm_counts: pd.DataFrame, genes, groups: pd.Series, title=None):
    genes = [g for g in genes if g in norm_counts.index]
    x = np.log2(norm_counts.loc[genes] + 1)
    z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1).replace(0, 1), axis=0)
    if len(z) > 2:
        z = z.iloc[leaves_list(linkage(z.values, "average", metric="correlation"))]
    order = groups.sort_values(kind="stable").index
    z = z[order]
    h = max(3.5, 0.17 * len(z) + 1.4)
    fig, ax = plt.subplots(figsize=(7, h))
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_yticks(range(len(z)))
    ax.set_yticklabels(z.index, fontsize=6.5 if len(z) > 40 else 8)
    ax.set_xticks(range(z.shape[1]))
    ax.set_xticklabels(order, rotation=90, fontsize=7.5)
    gvals = groups[order].values
    pal = dict(zip(sorted(set(gvals)), ["#7BA05B", "#B5835A", "#8E7CC3", "#E0A030"]))
    for i, g in enumerate(gvals):
        ax.add_patch(plt.Rectangle((i - 0.5, -1.6), 1, 0.9, color=pal[g], clip_on=False))
    for g, c in pal.items():
        ax.bar(0, 0, color=c, label=g)
    ax.legend(frameon=False, fontsize=8, loc="lower left", bbox_to_anchor=(1.12, 0))
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="z-score")
    ax.set_title(title or "Heatmap", pad=18)
    fig.tight_layout()
    return fig


def panel_plot(panel: pd.DataFrame, title=None):
    p = panel[panel["status"] != "미검출"].copy()
    p = p.iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.2, max(4, 0.19 * len(p) + 1.2)))
    cols = [COLORS.get(s, NS) for s in p["status"]]
    y = np.arange(len(p))
    ax.barh(y, p["log2fc"], color=cols, height=0.7)
    marks = {"up": "▲", "down": "▼", "var": "◆"}
    labels = [f"{g}  {marks[e]}" for g, e in zip(p["gene"], p["expected"])]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(0, c="#555", lw=0.8)
    # 범주 경계
    cats = p["category"].values
    for i in range(1, len(cats)):
        if cats[i] != cats[i - 1]:
            ax.axhline(i - 0.5, c="#ddd", lw=0.8)
    starts = {}
    for i, c in enumerate(cats):
        starts.setdefault(c, []).append(i)
    xmax = np.nanmax(np.abs(p["log2fc"])) if len(p) else 1
    for c, idx in starts.items():
        ax.text(xmax * 1.05, np.mean(idx), c, fontsize=7.5, va="center", color="#444")
    ax.set_xlim(-xmax * 1.1, xmax * 1.9)
    ax.set_xlabel("log2 Fold Change (노화 / 젊은 피부)")
    ax.set_title(title or "피부노화 유전자 패널 (▲증가·▼감소 문헌 경향, ◆상황의존)")
    fig.tight_layout()
    return fig


def ora_dotplot(df: pd.DataFrame, top_n=15, title=None, color="#C8463D"):
    d = df.nsmallest(top_n, "padj").iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.6, max(3, 0.32 * len(d) + 1.2)))
    if d.empty:
        ax.text(0.5, 0.5, "유의한 경로 없음", ha="center", va="center")
        ax.axis("off")
        return fig
    nl = -np.log10(d["padj"].clip(lower=1e-300))
    sc = ax.scatter(d["fold_enrichment"], range(len(d)), s=d["overlap"] * 12 + 15, c=nl,
                    cmap="Reds" if color == "#C8463D" else "Blues", edgecolor="#333", lw=0.4)
    terms = [t if len(t) < 55 else t[:52] + "…" for t in d["term"]]
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(terms, fontsize=7.5)
    ax.set_xlabel("Fold enrichment")
    fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02, label="-log10(FDR)")
    ax.set_title(title or "Enrichment")
    fig.tight_layout()
    return fig


def pca_plot(pcs, var, groups, title=None):
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    pal = ["#7BA05B", "#B5835A", "#8E7CC3", "#E0A030"]
    for i, g in enumerate(sorted(groups.unique())):
        s = groups[groups == g].index
        ax.scatter(pcs.loc[s, "PC1"], pcs.loc[s, "PC2"], s=60, c=pal[i % 4], label=g, edgecolor="#333", lw=0.5)
        for k in s:
            ax.annotate(k, (pcs.loc[k, "PC1"], pcs.loc[k, "PC2"]), fontsize=7, xytext=(4, 3),
                        textcoords="offset points")
    ax.set_xlabel(f"PC1 ({var[0]:.1%})")
    ax.set_ylabel(f"PC2 ({var[1]:.1%})")
    ax.set_title(title or "PCA (분산 상위 1,000 유전자)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig


def score_plot(scores, title=None):
    d = scores.sort_values("노화점수_LOOCV")
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    pal = ["#7BA05B", "#B5835A", "#8E7CC3", "#E0A030"]
    cmap = {g: pal[i % 4] for i, g in enumerate(sorted(d["group"].unique()))}
    ax.bar(d["sample"], d["노화점수_LOOCV"], color=[cmap[g] for g in d["group"]])
    ax.axhline(0.5, ls="--", c="#777", lw=0.8)
    ax.set_ylim(0, 1)
    ax.set_ylabel("전사체 피부노화 점수 (LOOCV)")
    for g, c in cmap.items():
        ax.bar(0, 0, color=c, label=g)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(title or "AI 피부노화 점수 (샘플별, 학습에 쓰지 않은 상태에서 예측)")
    ax.tick_params(axis="x", rotation=90, labelsize=7.5)
    fig.tight_layout()
    return fig


def importance_plot(imp, top_n=20, title=None):
    d = imp.head(top_n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(6, max(3.5, 0.25 * len(d) + 1)))
    ax.barh(d["gene"], d["weight"], color=[UP if w > 0 else DOWN for w in d["weight"]])
    ax.axvline(0, c="#555", lw=0.8)
    for lab in ax.get_yticklabels():
        if lab.get_text() in set(imp.loc[imp["피부노화 패널"], "gene"]):
            lab.set_fontweight("bold")
    ax.set_xlabel("표준화 계수 (+ 노화 쪽 / − 젊음 쪽)")
    ax.set_title(title or f"모델 중요 유전자 상위 {top_n} (굵게 = 피부노화 패널)")
    fig.tight_layout()
    return fig


def compare_scatter(merged: pd.DataFrame, labels, r=None, title=None):
    xa, xb = merged.columns[0], merged.columns[1]
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    both = (merged.iloc[:, 2] != "NS") & (merged.iloc[:, 3] != "NS")
    ax.scatter(merged[xb], merged[xa], s=5, c="#BDBDBD", lw=0, alpha=0.5, label="기타")
    ax.scatter(merged.loc[both, xb], merged.loc[both, xa], s=8, c="#5B3E96", lw=0,
               alpha=0.8, label="양쪽 모두 DEG")
    lim = np.nanpercentile(np.abs(merged[[xa, xb]].values), 99.8) * 1.1
    ax.plot([-lim, lim], [-lim, lim], ls="--", c="#777", lw=0.8)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel(f"log2FC — {labels[1]}")
    ax.set_ylabel(f"log2FC — {labels[0]}")
    ax.set_title(title or (f"log2FC 비교 (Pearson r = {r:.3f})" if r is not None else "log2FC 비교"))
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------- SkinAge-AI 예측 시각화
def prevention_bars(profile, title=None):
    d = profile.iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.2, max(3.2, 0.42 * len(d) + 1.2)))
    y = np.arange(len(d))
    ax.barh(y + 0.18, d["stress_only"], height=0.36, color="#B0B0B0", label="스트레스 단독")
    ax.barh(y - 0.18, d["with_compound"], height=0.36, color="#2F8F6A", label="스트레스 + 후보물질",
            xerr=d["ci95"], error_kw={"elinewidth": 0.8, "ecolor": "#555"})
    ax.set_yticks(y)
    ax.set_yticklabels([f"{b}  ({a})" for b, a in zip(d.biomarker, d.axis)], fontsize=8)
    ax.axvline(0, c="#555", lw=0.8)
    for yi, p in zip(y, d["prevention"]):
        if np.isfinite(p):
            ax.text(ax.get_xlim()[1], yi, f" 예방 {p:.0f}%", va="center", fontsize=7.5, color="#2F8F6A")
    ax.set_xlabel("무처리 대조 대비 변화 (%)  ·  오차막대 = 95% 예측구간")
    ax.set_title(title or "기전지표 예측")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    return fig


def radar(profile, title=None):
    axes_ = profile.groupby("axis", sort=False)
    labels, s_vals, c_vals = [], [], []
    for a, g in axes_:
        labels.append(a)
        s_vals.append(np.mean(np.abs(g.stress_only)))
        c_vals.append(np.mean(np.abs(g.with_compound)))
    mx = max(max(s_vals), 1)
    s_vals = [v / mx * 100 for v in s_vals]
    c_vals = [v / mx * 100 for v in c_vals]
    n = len(labels)
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    ang += ang[:1]
    s_vals += s_vals[:1]
    c_vals += c_vals[:1]
    fig = plt.figure(figsize=(5.2, 5.0))
    ax = fig.add_subplot(111, polar=True)
    ax.plot(ang, s_vals, color="#9A9A9A", lw=1.5, label="스트레스 단독")
    ax.fill(ang, s_vals, color="#9A9A9A", alpha=0.15)
    ax.plot(ang, c_vals, color="#2F8F6A", lw=2, label="+ 후보물질")
    ax.fill(ang, c_vals, color="#2F8F6A", alpha=0.25)
    ax.set_xticks(ang[:-1])
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yticklabels([])
    ax.set_title(title or "기전 축별 손상 프로파일 (상대값)", pad=18)
    ax.legend(frameon=False, fontsize=8, loc="upper right", bbox_to_anchor=(1.25, 1.12))
    fig.tight_layout()
    return fig


def dose_curve(dr, fit, unit, title=None):
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    c = dr.conc.values
    m = c > 0
    ax.scatter(c[m], dr.score.values[m], s=28, c="#2F8F6A", zorder=3, label="AI 예측 예방점수")
    xs = np.geomspace(c[m].min(), c[m].max(), 200)
    ys = fit["top"] * xs ** fit["hill"] / (fit["ec50"] ** fit["hill"] + xs ** fit["hill"])
    ax.plot(xs, ys, c="#1F5F47", lw=1.6, label=f"4PL 적합 (EC₅₀ {fit['ec50']:.3g} {unit})")
    ax.axvline(fit["ec50"], ls="--", c="#888", lw=0.8)
    if fit["rec80"] <= c.max() * 1.5:
        ax.axvline(fit["rec80"], ls=":", c="#C2780C", lw=1.2, label=f"권장 농도(최대효과 80%) {fit['rec80']:.3g} {unit}")
    ax.set_xscale("log")
    ax.set_xlabel(f"농도 ({unit}, log)")
    ax.set_ylabel("예방점수 (%)")
    ax.set_title(title or "용량-반응")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig


def ranking_bars(df, title=None):
    d = df.sort_values("예방점수").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(6.6, max(2.6, 0.5 * len(d) + 1)))
    ax.barh(d["후보물질"], d["예방점수"], color="#2F8F6A")
    for i, v in enumerate(d["예방점수"]):
        ax.text(v, i, f" {v:.0f}%", va="center", fontsize=8)
    ax.set_xlabel("질환 관련 기전 축 가중 예방점수 (%)")
    ax.set_title(title or "후보물질 비교")
    ax.tick_params(axis="y", labelsize=8)
    fig.tight_layout()
    return fig


def parity_plot(pairs, title=None):
    bms = list(pairs.biomarker.unique())
    cols = 4
    rows = int(np.ceil(len(bms) / cols))
    fig, axs = plt.subplots(rows, cols, figsize=(11, 2.6 * rows))
    for ax, b in zip(np.ravel(axs), bms):
        d = pairs[pairs.biomarker == b]
        ax.scatter(d.observed, d.predicted, s=5, alpha=0.5, c="#2F8F6A", lw=0)
        lo, hi = min(d.observed.min(), d.predicted.min()), max(d.observed.max(), d.predicted.max())
        ax.plot([lo, hi], [lo, hi], ls="--", c="#888", lw=0.8)
        ax.set_title(b, fontsize=8)
        ax.tick_params(labelsize=6)
    for ax in np.ravel(axs)[len(bms):]:
        ax.axis("off")
    fig.suptitle(title or "교차검증: 관측값(x) vs AI 예측값(y)", fontsize=10)
    fig.tight_layout()
    return fig


def importance_group_bars(s, title=None):
    fig, ax = plt.subplots(figsize=(5.6, 2.8))
    d = s.iloc[::-1]
    ax.barh(d.index, d.values, color="#4B5BAB")
    ax.set_xlabel("예측 기여도 (%) · 순열 중요도")
    ax.set_title(title or "예측에 영향을 준 요인")
    fig.tight_layout()
    return fig
