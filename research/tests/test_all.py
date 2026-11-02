"""SkinAge-Omics 자동 테스트.  실행:  python tests/test_all.py"""
import io
import runpy
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

from core import analysis as A, ml  # noqa: E402

DEMO = ROOT / "data" / "demo"
counts = pd.read_csv(DEMO / "demo_counts.csv", index_col=0)
meta = pd.read_csv(DEMO / "demo_metadata.csv", index_col=0)
truth = pd.read_csv(DEMO / "demo_truth.csv").set_index("gene")
results = []


def check(name, cond, detail=""):
    results.append((name, bool(cond), detail))
    print(("PASS" if cond else "FAIL"), name, detail)


# 1. BH 보정이 R p.adjust(method="BH") 과 같은 값을 주는지 (손계산 예시)
p = np.array([0.01, 0.04, 0.03, 0.005])
check("BH 보정 정확성", np.allclose(A.bh_adjust(p), [0.02, 0.04, 0.04, 0.02]))

# 2. size factor: 샘플을 2배로 만들면 size factor 도 2배
c2 = counts.copy(); c2["Y1"] = c2["Y1"] * 2
sf1, sf2 = A.size_factors(counts), A.size_factors(c2)
check("median-of-ratios 정규화 비례성", abs(sf2["Y1"] / sf1["Y1"] / (sf2["Y2"] / sf1["Y2"]) - 2) < 1e-6)

# 3. DEG: 정답 복원 + 귀무 유전자 p-value 보정
for meth in ("nb_wald", "welch"):
    d = A.classify(A.deg_from_counts(counts, meta, "group", "young", "aged", meth)).set_index("gene")
    t = truth.loc[d.index]
    sig = d.status != "NS"
    null = t.true_log2fc == 0
    prec = (sig & t.is_true_deg).sum() / max(sig.sum(), 1)
    big = t.true_log2fc.abs() >= 1
    rec = (sig & big).sum() / big.sum()
    fpr = (d.pvalue[null] < 0.05).mean()
    check(f"DEG[{meth}] 정밀도 ≥ 0.9", prec >= 0.9, f"precision={prec:.3f}")
    check(f"DEG[{meth}] 귀무 p<0.05 비율 ≤ 0.06", fpr <= 0.06, f"rate={fpr:.4f}")
    check(f"DEG[{meth}] 방향 정확(MMP1↑, COL1A1↓)", d.loc["MMP1", "log2fc"] > 0 and d.loc["COL1A1", "log2fc"] < 0)
    print(f"      recall(|true lfc|≥1)={rec:.3f}, DEG={sig.sum()}")

# 4. 결과표 표준화: DESeq2/edgeR/limma 형식
deg = A.deg_from_counts(counts, meta, "group", "young", "aged")
fmt = {
    "DESeq2": deg.rename(columns={"gene": "Unnamed: 0", "log2fc": "log2FoldChange", "basemean": "baseMean"}),
    "edgeR": deg.rename(columns={"gene": "genes", "log2fc": "logFC", "pvalue": "PValue", "padj": "FDR"}),
    "limma": deg.rename(columns={"gene": "ID", "log2fc": "logFC", "pvalue": "P.Value", "padj": "adj.P.Val"}),
}
for k, v in fmt.items():
    s = A.standardize_deg(v)
    check(f"입력 형식 인식: {k}", len(s) == len(deg) and np.allclose(s.log2fc.values, deg.log2fc.values))

# 5. ORA: 심어 둔 경로가 검출되는지
gs = A.read_gmt(ROOT / "data/genesets/KEGG_MSigDB_v7.5.1.gmt")
dd = A.classify(deg)
dn = A.ora(dd.gene[dd.status == "Down"], dd.gene, gs)
check("KEGG ORA: ECM-receptor interaction (Down) FDR<0.05",
      (dn[dn.term == "ECM RECEPTOR INTERACTION"].padj < 0.05).any())
hs = A.read_gmt(ROOT / "data/genesets/Hallmark_v7.0.gmt")
up = A.ora(dd.gene[dd.status == "Up"], dd.gene, hs)
check("Hallmark ORA: INFLAMMATORY RESPONSE (Up) FDR<0.05",
      (up[up.term == "HALLMARK INFLAMMATORY RESPONSE"].padj < 0.05).any())

# 6. 비교 모듈: 자기 자신과 비교하면 완전 일치
r = A.compare_results(deg, deg)
check("비교 모듈 자기일치 (r=1, Jaccard=1)", abs(r["pearson_lfc"] - 1) < 1e-9 and r["jaccard"] == 1)

# 7. AI 모듈
s, imp, met, model = ml.aging_classifier(A.normalize(counts), meta["group"], "aged")
check("AI LOOCV 정확도 ≥ 0.9", met["loocv_accuracy"] >= 0.9, f"acc={met['loocv_accuracy']:.3f}, auc={met['loocv_auc']:.3f}")
# 라벨을 섞으면 성능이 우연 수준으로 떨어져야 함(누수 없음 확인)
rng = np.random.default_rng(1)
perm = pd.Series(rng.permutation(meta["group"].values), index=meta.index)
_, _, metp, _ = ml.aging_classifier(A.normalize(counts), perm, "aged")
check("AI 라벨 셔플 대조군 (정보누수 없음: 정확도 ≤ 0.75)", metp["loocv_accuracy"] <= 0.75,
      f"acc={metp['loocv_accuracy']:.3f}, auc={metp['loocv_auc']:.3f}")

# 8. 앱 스모크 테스트 (가짜 streamlit 으로 모든 탭 실행)
import st_stub  # noqa: E402
for mode_name, extra in [("데모 데이터", {}),
                         ("DEG 결과표 업로드", {"DEG 결과표 (DESeq2/edgeR/limma)": "deseq"})]:
    st_stub.OVERRIDES.clear(); st_stub.FIGS.clear()
    st_stub.OVERRIDES["입력 데이터"] = mode_name
    if extra:
        b = io.BytesIO(fmt["DESeq2"].to_csv(index=False).encode()); b.name = "deseq2.csv"
        st_stub.OVERRIDES["DEG 결과표 (DESeq2/edgeR/limma)"] = b
    else:
        b = io.BytesIO(deg.to_csv(index=False).encode()); b.name = "old.csv"
        st_stub.OVERRIDES["기존 DEG 결과표"] = b
    st_stub.install()
    try:
        runpy.run_path(str(ROOT / "app.py"), run_name="__main__")
        ok = True
    except st_stub.mod.Stop:
        ok = True
    except Exception as e:  # noqa: BLE001
        ok = False
        import traceback; traceback.print_exc()
    check(f"앱 실행: {mode_name}", ok, f"그림 {len(st_stub.FIGS)}개 생성")

# 8b. AI 예측 엔진
from core import predict as PR  # noqa: E402
_d = PR.simulate(900, seed=3)
_cv, _ = PR.cross_validate(_d, k=3, max_n=200)
check("AI 예측 엔진 교차검증 평균 R² ≥ 0.6", _cv.R2.mean() >= 0.6, f"mean R²={_cv.R2.mean():.3f}")
_ai = PR.SkinAgeAI(_d, max_n=250)
_p = PR.prevention_profile(_ai, "cell", "uv", 0.75, 4, "후보물질 C (콜라겐·항노화형)", 50, "pre")
_t0 = PR.true_effect("cell", "uv", 0.75, 4, "없음", 0, "pre")
_t1 = PR.true_effect("cell", "uv", 0.75, 4, "후보물질 C (콜라겐·항노화형)", 50, "pre")
_tp = np.array([(_t0[b] - _t1[b]) / _t0[b] * 100 for b in _p.biomarker])
check("예방률 예측이 정답과 같은 경향 (상관 ≥ 0.6)", np.corrcoef(_tp, _p.prevention)[0, 1] >= 0.6,
      f"r={np.corrcoef(_tp, _p.prevention)[0, 1]:.2f}")
_f = PR.fit_4pl([0, 1, 3, 10, 30, 100], [0, 9, 23, 50, 75, 90])
check("4PL 적합 EC50 합리성 (5~20)", 5 < _f["ec50"] < 20, f"EC50={_f['ec50']:.2f}")
check("실데이터·문헌 양식 인식", (PR.standardize_real(PR.template_df()).source == "문헌").sum() == 1)

# 9. 통합 앱(streamlit_app.py): 모든 모듈 + 실데이터 예시 모드
for module, extra in [("🤖 SkinAge-AI · 예방효과 예측", {}),
                      ("🤖 SkinAge-AI · 예방효과 예측", {"예측할 질환": "atopic", "개체": "mouse"}),
                      ("🧬 SkinAge-Omics · 전사체 분석", {"입력 데이터": "데모 데이터"}),
                      ("🧫 SkinCell ProtoLab · Wet-lab", {}), ("🎨 PigmentLab · 색소 분석", {}),
                      ("ⓘ 플랫폼 소개", {})]:
    st_stub.OVERRIDES.clear(); st_stub.FIGS.clear(); st_stub.CALLS.clear()
    st_stub.OVERRIDES["모듈 선택"] = module
    st_stub.OVERRIDES.update(extra)
    st_stub.install()
    try:
        runpy.run_path(str(ROOT / "streamlit_app.py"), run_name="__main__")
        ok = True
    except st_stub.mod.Stop:
        ok = True
    except Exception:  # noqa: BLE001
        ok = False
        import traceback; traceback.print_exc()
    label = module + (" / " + extra.get("입력 데이터", "") if extra else "")
    check(f"통합 앱 실행: {label}", ok, f"그림 {len(st_stub.FIGS)}개, html {st_stub.CALLS.count('components.html')}개")

n_fail = sum(not r[1] for r in results)
print(f"\n총 {len(results)}개 중 통과 {len(results) - n_fail}개")
pd.DataFrame(results, columns=["test", "pass", "detail"]).to_csv(ROOT / "results" / "test_report.csv", index=False)
sys.exit(1 if n_fail else 0)
