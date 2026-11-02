"""앱과 동일한 분석을 명령줄에서 한 번에 실행해 결과 그림/표를 저장한다.

사용 예)
  python scripts/run_pipeline.py                                   # 데모 데이터
  python scripts/run_pipeline.py --counts counts.csv --meta meta.csv --group group --ref young --test aged
  python scripts/run_pipeline.py ... --compare 기존_DESeq2_결과.csv  # 기존 분석과 비교
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from core import analysis as A, plots as P, ml as ML  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--counts", default=str(ROOT / "data/demo/demo_counts.csv"))
ap.add_argument("--meta", default=str(ROOT / "data/demo/demo_metadata.csv"))
ap.add_argument("--group", default="group")
ap.add_argument("--ref", default="young")
ap.add_argument("--test", default="aged")
ap.add_argument("--method", default="nb_wald", choices=["nb_wald", "welch"])
ap.add_argument("--padj", type=float, default=0.05)
ap.add_argument("--lfc", type=float, default=1.0)
ap.add_argument("--compare", default=None, help="기존 DEG 결과표 (없으면 Welch 결과와 비교)")
ap.add_argument("--out", default=str(ROOT / "results"))
a = ap.parse_args()

out = Path(a.out); (out / "figures").mkdir(parents=True, exist_ok=True); (out / "tables").mkdir(exist_ok=True)
counts = A.read_table(a.counts); counts = counts.set_index(counts.columns[0])
counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0).round().astype(int)
meta = A.read_table(a.meta); meta = meta.set_index(meta.columns[0]); meta[a.group] = meta[a.group].astype(str)
cols = [c for c in counts.columns if str(c) in meta.index and meta.loc[str(c), a.group] in (a.ref, a.test)]
norm = A.normalize(counts[cols]); norm.index = norm.index.astype(str).str.upper()
groups = meta.loc[cols, a.group]


def save(fig, name):
    fig.savefig(out / "figures" / name, dpi=200, bbox_inches="tight")


deg = A.classify(A.deg_from_counts(counts, meta, a.group, a.ref, a.test, a.method), a.padj, a.lfc)
deg.to_csv(out / "tables/deg_results.csv", index=False)
save(P.volcano(deg, a.padj, a.lfc, label_genes=["MMP1", "COL1A1", "CDKN2A", "IL6"],
               title=f"Volcano plot ({a.test} vs {a.ref})"), "01_volcano.png")
save(P.ma_plot(deg), "02_ma_plot.png")
top = deg[deg.status != "NS"].nsmallest(40, "padj").gene
save(P.heatmap(norm, list(top), groups, title="상위 40개 DEG (z-score)"), "03_heatmap_top40.png")

panel = A.aging_panel(deg, a.padj, a.lfc)
panel.to_csv(out / "tables/skin_aging_panel.csv", index=False)
A.panel_summary(panel).to_csv(out / "tables/skin_aging_panel_summary.csv", index=False)
save(P.panel_plot(panel), "04_skin_aging_panel.png")

summary = {"n_genes_tested": len(deg), "up": int((deg.status == "Up").sum()), "down": int((deg.status == "Down").sum())}
for key, fn in [("KEGG", "KEGG_MSigDB_v7.5.1.gmt"), ("GO_BP", "GO_Biological_Process_2021.gmt"),
                ("Hallmark", "Hallmark_v7.0.gmt")]:
    gs = A.read_gmt(ROOT / "data/genesets" / fn)
    for direc, col in (("Up", "#C8463D"), ("Down", "#2F6DB5")):
        r = A.ora(deg.gene[deg.status == direc], deg.gene, gs)
        r.to_csv(out / f"tables/ora_{key}_{direc}.csv", index=False)
        summary[f"{key}_{direc}_sig_terms"] = int((r.padj < 0.05).sum())
        summary[f"{key}_{direc}_top"] = r.head(5).term.tolist()
        if key in ("KEGG", "GO_BP"):
            save(P.ora_dotplot(r, 15, f"{key} · {direc} DEG", col), f"05_ora_{key}_{direc}.png")

if a.compare:
    other = A.standardize_deg(A.read_table(a.compare)); lab = "기존 분석"
else:
    other = A.deg_from_counts(counts, meta, a.group, a.ref, a.test,
                              "welch" if a.method == "nb_wald" else "nb_wald")
    lab = "Welch t-test" if a.method == "nb_wald" else "NB-Wald"
cmp = A.compare_results(deg, other, a.padj, a.lfc, ("SkinAge-Omics", lab))
save(P.compare_scatter(cmp["merged"], ("SkinAge-Omics", lab), cmp["pearson_lfc"]), "06_compare_scatter.png")
cmp["crosstab"].to_csv(out / "tables/compare_crosstab.csv")
summary["compare"] = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in cmp.items()
                      if k not in ("crosstab", "merged")}
summary["compare"]["against"] = lab

scores, imp, met, _ = ML.aging_classifier(norm, groups, a.test)
pcs, var = ML.pca(norm)
save(P.pca_plot(pcs, var, groups), "07_pca.png")
save(P.score_plot(scores), "08_ai_aging_score.png")
save(P.importance_plot(imp), "09_ai_feature_importance.png")
scores.to_csv(out / "tables/ai_aging_scores.csv", index=False)
imp.to_csv(out / "tables/ai_feature_importance.csv", index=False)
summary["ai"] = met

truth_f = ROOT / "data/demo/demo_truth.csv"
if a.counts.endswith("demo_counts.csv") and truth_f.exists():
    t = pd.read_csv(truth_f).set_index("gene").loc[deg.gene]
    sig = (deg.status != "NS").values
    big = (t.true_log2fc.abs() >= 1).values
    null = (t.true_log2fc == 0).values
    summary["simulation_truth"] = {
        "precision": round(float((sig & t.is_true_deg.values).sum() / sig.sum()), 4),
        "recall_true_lfc_ge1": round(float((sig & big).sum() / big.sum()), 4),
        "null_p_lt_0.05_rate": round(float((deg.pvalue.values[null] < 0.05).mean()), 4),
        "lfc_pearson_vs_truth_trueDEG": round(float(pd.Series(deg.log2fc.values[t.is_true_deg.values]).corr(
            pd.Series(t.true_log2fc.values[t.is_true_deg.values]))), 4),
    }
(out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2, default=str))
