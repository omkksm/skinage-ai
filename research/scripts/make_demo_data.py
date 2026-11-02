"""데모용 합성 RNA-seq count 데이터 생성 (젊은 피부 6 vs 노화 피부 6).

※ 실제 환자/공개 데이터가 아닌 '정답을 알고 있는' 시뮬레이션 데이터입니다.
   파이프라인 동작 확인과 정답 복원율(민감도/정밀도) 검증용으로만 사용합니다.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from core.analysis import read_gmt  # noqa: E402
from core.skin_aging_genes import SKIN_AGING_PANEL  # noqa: E402

rng = np.random.default_rng(20261102)
GS = ROOT / "data" / "genesets"
go = read_gmt(GS / "GO_Biological_Process_2021.gmt")
hall = read_gmt(GS / "Hallmark_v7.0.gmt")
kegg = read_gmt(GS / "KEGG_MSigDB_v7.5.1.gmt")

panel = {g: d for cat in SKIN_AGING_PANEL.values() for g, d in cat.items()}
pool = sorted(set().union(*go.values()) | set().union(*kegg.values()))
n_genes = 12000
others = [g for g in pool if g not in panel]
genes = sorted(set(panel) | set(rng.choice(others, n_genes - len(panel), replace=False)))

base = pd.Series(np.exp(rng.normal(4.3, 1.9, len(genes))), index=genes).clip(0.5, 2e5)
lfc = pd.Series(0.0, index=genes)

# 1) 피부노화 패널 유전자: 문헌 경향 방향으로 효과 부여 (일부는 효과 없음으로 남겨 현실성 확보)
for g, d in panel.items():
    if d == "var" or rng.random() < 0.15:
        continue
    lfc[g] = (1 if d == "up" else -1) * rng.uniform(1.0, 2.4)
    if g in ("COL1A1", "COL1A2", "COL3A1", "ELN"):
        base[g] *= 30  # 진피에서 고발현
# 2) 경로 수준 신호
def plant(geneset, n, sign, lo, hi):
    cand = [g for g in geneset if g in lfc.index and lfc[g] == 0]
    for g in rng.choice(cand, min(n, len(cand)), replace=False):
        lfc[g] = sign * rng.uniform(lo, hi)
plant(hall["HALLMARK_INFLAMMATORY_RESPONSE"], 45, +1, 0.9, 2.0)
plant(hall["HALLMARK_INTERFERON_GAMMA_RESPONSE"], 25, +1, 0.8, 1.6)
plant(kegg["KEGG_ECM_RECEPTOR_INTERACTION"], 30, -1, 0.9, 2.0)
plant(hall["HALLMARK_OXIDATIVE_PHOSPHORYLATION"], 40, -1, 0.8, 1.5)
plant(kegg["KEGG_CELL_CYCLE"], 20, +1, 0.8, 1.6)
# 3) 무작위 배경 DEG
rest = lfc.index[lfc == 0]
for g in rng.choice(rest, 220, replace=False):
    lfc[g] = rng.choice([-1, 1]) * rng.uniform(1.0, 2.5)

samples = [f"Y{i}" for i in range(1, 7)] + [f"A{i}" for i in range(1, 7)]
group = ["young"] * 6 + ["aged"] * 6
age = list(rng.integers(21, 35, 6)) + list(rng.integers(62, 78, 6))
sf = rng.uniform(0.7, 1.4, 12)
alpha = 0.06 + 1.5 / base.values  # 평균-분산 관계
counts = np.zeros((len(genes), 12), dtype=int)
for j in range(12):
    mu = base.values * sf[j] * (2 ** lfc.values if group[j] == "aged" else 1)
    mu = mu * np.exp(rng.normal(0, 0.12, len(genes)))  # 샘플 간 생물학적 잡음
    r = 1 / alpha
    counts[:, j] = rng.negative_binomial(r, r / (r + mu))

out = ROOT / "data" / "demo"
out.mkdir(parents=True, exist_ok=True)
pd.DataFrame(counts, index=pd.Index(genes, name="gene"), columns=samples).to_csv(out / "demo_counts.csv")
pd.DataFrame({"sample": samples, "group": group, "age": age, "site": "forearm (sim)"}).to_csv(
    out / "demo_metadata.csv", index=False)
truth = pd.DataFrame({"gene": genes, "true_log2fc": lfc.values,
                      "is_true_deg": (lfc.abs() >= 0.8).values})
truth.to_csv(out / "demo_truth.csv", index=False)
print(f"genes={len(genes)} true DEG={int(truth.is_true_deg.sum())} -> {out}")
