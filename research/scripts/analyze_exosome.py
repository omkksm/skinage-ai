"""엑소좀(EV17) 처리 B16 세포 RNA-seq 실데이터 분석 + 기존 edgeR/clusterProfiler 결과와 비교.

입력(압축 해제 폴더):
  24.11.21_RAW_data_*.csv            edgeR 전체 결과 + raw count (α-MSH vs EV17 5ppm, 각 n=1)
  24.11_21_1.5_fold_change_FDR_0.05.csv   기존 DEG (|FC|≥1.5, FDR<0.05)
  GO_and_KEGG/...                    기존 clusterProfiler GO·KEGG 결과
사용:  python scripts/analyze_exosome.py <압축해제폴더> <출력폴더>
"""
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from core import analysis as A, plots as P  # noqa: E402

src, out = Path(sys.argv[1]), Path(sys.argv[2])
(out / "그림").mkdir(parents=True, exist_ok=True)
(out / "표").mkdir(exist_ok=True)
LFC = np.log2(1.5)
PADJ = 0.05


def save(fig, name):
    fig.savefig(out / "그림" / name, dpi=200, bbox_inches="tight")


raw = pd.read_csv(glob.glob(str(src / "24.11.21_RAW_data_*.csv"))[0])
old_deg = pd.read_csv(glob.glob(str(src / "24.11_21_1.5_fold_change*.csv"))[0])
raw = raw.dropna(subset=["external_gene_name"])
raw["gene"] = raw["external_gene_name"].astype(str)
counts = raw.groupby("gene")[["RawCount_a.MSH", "RawCount_EV17_5ppm"]].sum()
counts.columns = ["aMSH", "EV17"]
meta = pd.DataFrame({"group": ["aMSH", "EV17"]}, index=["aMSH", "EV17"])

# 1) SkinAge-Omics 로 재분석 (무반복 → 고정 dispersion 0.01, BCV 0.1)
ours = A.deg_from_counts(counts, meta, "group", "aMSH", "EV17", "nb_wald", min_count=10)
ours = A.classify(ours, PADJ, LFC)
ours.to_csv(out / "표" / "SkinAge-Omics_DEG_전체.csv", index=False, encoding="utf-8-sig")

# 2) 기존 edgeR 결과 표준화
edger = raw.rename(columns={"logFC": "log2fc", "PValue": "pvalue", "FDR": "padj"})
edger = edger[["gene", "log2fc", "pvalue", "padj", "logCPM"]].drop_duplicates("gene")
edger["gene"] = edger["gene"].str.upper()
edger["basemean"] = 2 ** edger["logCPM"]
edger = A.classify(edger.drop(columns="logCPM"), PADJ, LFC)
cmp = A.compare_results(ours, edger, PADJ, LFC, ("SkinAge-Omics", "기존 edgeR"))
save(P.compare_scatter(cmp["merged"], ("SkinAge-Omics", "기존 edgeR"), cmp["pearson_lfc"]), "01_기존edgeR_비교_log2FC.png")
cmp["crosstab"].to_csv(out / "표" / "기존edgeR_비교_교차표.csv", encoding="utf-8-sig")

# 기존 DEG 목록(파일)과 edgeR 전체결과 기준 재분류가 같은지 확인
old_set = set(old_deg["ID"].astype(str).str.upper())
edger_set = set(edger.gene[edger.status != "NS"])

# 3) 시각화
save(P.volcano(ours, PADJ, LFC, title="EV17 5ppm vs α-MSH (SkinAge-Omics)"), "02_volcano.png")
save(P.ma_plot(ours, title="MA plot (EV17 vs α-MSH)"), "03_ma_plot.png")
norm = A.normalize(counts)
norm.index = norm.index.str.upper()
top = ours[ours.status != "NS"].nsmallest(40, "padj").gene
fig = P.heatmap(norm, list(top), meta["group"], title="상위 40개 DEG (2 샘플, 행 z-score)")
save(fig, "04_heatmap_top40.png")

# 4) 피부노화·색소 패널 (마우스 유전자 기호 → 대문자로 사람 기호와 대응)
panel = A.aging_panel(ours, PADJ, LFC)
panel.to_csv(out / "표" / "피부노화_색소_패널.csv", index=False, encoding="utf-8-sig")
save(P.panel_plot(panel, "피부노화·색소 유전자 패널 (EV17 / α-MSH)"), "05_panel.png")
melano = ["TYR", "TYRP1", "DCT", "MITF", "PMEL", "MLANA", "MC1R", "SLC45A2", "OCA2", "RAB27A", "GPNMB", "SOX10"]
mel = ours.set_index("gene").reindex(melano)[["log2fc", "padj", "status", "basemean"]].reset_index()
mel.to_csv(out / "표" / "멜라닌생성_핵심유전자.csv", index=False, encoding="utf-8-sig")

# 5) GO·KEGG ORA (앱 내장 유전자세트) 를 기존 DEG 목록에 적용 → 기존 clusterProfiler 결과와 비교
gs_go = A.read_gmt(ROOT / "data/genesets/GO_Biological_Process_2021.gmt")
gs_kegg = A.read_gmt(ROOT / "data/genesets/KEGG_MSigDB_v7.5.1.gmt")
universe = edger.gene
res = {}
for direc in ("Up", "Down"):
    genes = edger.gene[edger.status == direc]
    for key, gs in (("GO_BP", gs_go), ("KEGG", gs_kegg)):
        r = A.ora(genes, universe, gs)
        r.to_csv(out / "표" / f"ORA_{key}_{direc}.csv", index=False, encoding="utf-8-sig")
        res[(key, direc)] = r
        save(P.ora_dotplot(r, 15, f"{key} · {direc} DEG (EV17 vs α-MSH)",
                           "#C8463D" if direc == "Up" else "#2F6DB5"), f"06_ORA_{key}_{direc}.png")


def read_cp(path):
    return pd.read_csv(path)


go_cmp = {}
for direc, pat in (("Up", "GO_and_KEGG/GO/up/*bp*"), ("Down", "GO_and_KEGG/GO/down/*bp*")):
    cp = read_cp(glob.glob(str(src / pat))[0])
    cp_sig = cp[cp["p.adjust"] < 0.05]
    cp_ids = list(cp_sig["ID"])
    ours_go = res[("GO_BP", direc)].copy()
    ours_go["GOID"] = ours_go["term"].str.extract(r"\((GO:\d+)\)")
    ours_sig = ours_go[ours_go.padj < 0.05]
    ours_ids = set(ours_sig.GOID.dropna())
    tested = set(ours_go.GOID.dropna())
    top20 = cp_ids[:20]
    hit_top20 = [g for g in top20 if g in set(ours_go.GOID[ours_go.pvalue < 0.05])]
    go_cmp[direc] = {
        "기존_유의GO수": len(cp_ids), "앱_유의GO수": len(ours_ids),
        "공통_유의GO수": len(set(cp_ids) & ours_ids),
        "기존_상위20중_앱에서_p<0.05": f"{len(hit_top20)}/{len(top20)}",
        "기존_상위5": list(cp_sig["Description"].head(5)),
        "앱_상위5": [re.sub(r" \(GO:\d+\)", "", t) for t in ours_sig.term.head(5)],
    }

kegg_cp = read_cp(glob.glob(str(src / "GO_and_KEGG/KEGG/*kegg_down*"))[0])
kegg_ours = res[("KEGG", "Down")]
kegg_check = []
for d in kegg_cp["Description"]:
    name = d.split(" - ")[0].upper().replace("-", " ").replace(",", "")
    hit = kegg_ours[kegg_ours.term.str.upper().str.replace("-", " ") == name]
    kegg_check.append({"기존_KEGG(down)": d.split(" - ")[0],
                       "앱_FDR": float(hit.padj.iloc[0]) if len(hit) else None,
                       "앱_p": float(hit.pvalue.iloc[0]) if len(hit) else None,
                       "비고": "" if len(hit) else "앱 KEGG(MSigDB 2011판)에 없는 경로"})
pd.DataFrame(kegg_check).to_csv(out / "표" / "KEGG_down_기존대비.csv", index=False, encoding="utf-8-sig")

summary = {
    "샘플": "α-MSH 대조 1 vs EV17 5ppm 1 (B16 계열 마우스 세포로 추정, 무반복)",
    "검정유전자": int(len(ours)),
    "앱_DEG": {"Up": int((ours.status == "Up").sum()), "Down": int((ours.status == "Down").sum())},
    "기존_DEG_파일": {"총": len(old_deg), "Up": int((old_deg.logFC > 0).sum()), "Down": int((old_deg.logFC < 0).sum())},
    "기존DEG파일_vs_edgeR전체결과_재분류_일치": len(old_set & edger_set) / max(len(old_set | edger_set), 1),
    "기존edgeR_비교": {k: (round(v, 4) if isinstance(v, float) else v) for k, v in cmp.items() if k not in ("crosstab", "merged")},
    "GO_BP_비교": go_cmp,
    "KEGG_down_비교": kegg_check,
    "멜라닌생성유전자": mel.round(4).to_dict(orient="records"),
}
(out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2, default=str))
