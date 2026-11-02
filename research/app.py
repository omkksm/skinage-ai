"""SkinAge-Omics: RNA-seq 기반 피부노화 전사체 분석 앱 (Streamlit)

실행:  streamlit run app.py
개발:  문경미
"""
from pathlib import Path
import io

import numpy as np
import pandas as pd
import streamlit as st

from core import analysis as A
from core import plots as P
from core import ml as ML
from core.skin_aging_genes import SKIN_AGING_PANEL

ROOT = Path(__file__).resolve().parent
DEMO = ROOT / "data" / "demo"
GS_DIR = ROOT / "data" / "genesets"
GENESETS = {
    "KEGG (MSigDB c2.cp.kegg v7.5.1)": GS_DIR / "KEGG_MSigDB_v7.5.1.gmt",
    "GO Biological Process 2021": GS_DIR / "GO_Biological_Process_2021.gmt",
    "MSigDB Hallmark v7.0": GS_DIR / "Hallmark_v7.0.gmt",
}
VERSION = "1.0.0"

EXO = ROOT / "data" / "exosome"
if not globals().get("_SUITE"):
    st.set_page_config(page_title="SkinAge-Omics", page_icon="🧬", layout="wide")


@st.cache_data(show_spinner=False)
def load_gmt(path: str):
    return A.read_gmt(path)


@st.cache_data(show_spinner=False)
def run_deg(counts, meta, group_col, ref, test, method, min_count):  # noqa: D401
    return A.deg_from_counts(counts, meta, group_col, ref, test, method, min_count)


def fig_download(fig, name, key):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    st.download_button("그림 PNG 저장", buf.getvalue(), file_name=name, mime="image/png", key=key)
    P.plt.close(fig)


def csv_download(df, name, key, label="표 CSV 저장"):
    st.download_button(label, df.to_csv(index=False).encode("utf-8-sig"), file_name=name,
                       mime="text/csv", key=key)


# ------------------------------------------------------------------ 사이드바: 데이터 입력
if not globals().get("_SUITE"):
    st.sidebar.title("🧬 SkinAge-Omics")
    st.sidebar.caption(f"피부노화 전사체 분석 v{VERSION}")

EXO_MODE = "실데이터 예시: 엑소좀 RNA-seq"
_modes = ["데모 데이터", "Count matrix 업로드", "DEG 결과표 업로드"]
if (EXO / "exosome_counts.csv").exists():  # 연구실 원자료는 배포본에 포함하지 않음(원자료 보호)
    _modes.insert(1, EXO_MODE)
mode = st.sidebar.radio("입력 데이터", _modes,
                        help="Count matrix를 올리면 앱 안에서 DEG 분석을 수행합니다. "
                             "DESeq2/edgeR/limma 결과표를 바로 올려 시각화만 할 수도 있습니다.")

st.sidebar.subheader("DEG 기준")
padj_cut = st.sidebar.number_input("FDR (padj) 기준", 0.0001, 0.5, 0.05, 0.01, format="%.4f")
lfc_cut = st.sidebar.number_input("|log2FC| 기준", 0.0, 5.0, 0.585 if mode == EXO_MODE else 1.0, 0.25,
                                  help="1.0 = 2배, 0.585 = 1.5배")

deg = None
norm = None
groups = None
counts = meta = None
no_rep = False

if mode in ("데모 데이터", EXO_MODE, "Count matrix 업로드"):
    if mode == EXO_MODE:
        counts = pd.read_csv(EXO / "exosome_counts.csv", index_col=0)
        meta = pd.read_csv(EXO / "exosome_metadata.csv", index_col=0)
        st.sidebar.info("본 연구실 실데이터: 마우스 색소세포, α-MSH 대조 vs EV17 엑소좀 5 ppm (각 1샘플). "
                        "기존 edgeR 결과(2024.11)와 ④ 탭에서 비교할 수 있습니다.")
    elif mode == "데모 데이터":
        counts = pd.read_csv(DEMO / "demo_counts.csv", index_col=0)
        meta = pd.read_csv(DEMO / "demo_metadata.csv", index_col=0)
        st.sidebar.info("합성 데모 데이터: 젊은 피부 6 vs 노화 피부 6 (정답이 알려진 시뮬레이션)")
    else:
        fc = st.sidebar.file_uploader("Count matrix (행=유전자, 열=샘플)", type=["csv", "tsv", "txt", "xlsx"])
        fm = st.sidebar.file_uploader("샘플 정보 (첫 열=샘플명)", type=["csv", "tsv", "txt", "xlsx"])
        if fc is not None and fm is not None:
            counts = A.read_table(fc)
            counts = counts.set_index(counts.columns[0])
            counts.index = counts.index.astype(str)
            counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0).round().astype(int)
            counts = counts.groupby(level=0).sum()
            meta = A.read_table(fm)
            meta = meta.set_index(meta.columns[0])
    if counts is not None:
        cand = [c for c in meta.columns if meta[c].nunique() >= 2 and meta[c].nunique() <= 10]
        if not cand:
            st.error("샘플 정보에 그룹 컬럼(2~10개 수준)이 없습니다.")
            st.stop()
        group_col = st.sidebar.selectbox("그룹 컬럼", cand, index=cand.index("group") if "group" in cand else 0)
        levels = sorted(meta[group_col].astype(str).unique())
        meta[group_col] = meta[group_col].astype(str)
        dref = levels.index("young") if "young" in levels else (levels.index("aMSH") if "aMSH" in levels else 0)
        ref = st.sidebar.selectbox("기준 그룹 (예: 젊은 피부)", levels, index=dref)
        others = [l for l in levels if l != ref]
        test = st.sidebar.selectbox("비교 그룹 (예: 노화 피부)", others,
                                    index=others.index("aged") if "aged" in others else 0)
        _n = meta[group_col].value_counts()
        no_rep = min(_n.get(ref, 0), _n.get(test, 0)) < 2
        method = st.sidebar.selectbox("DEG 방법", ["nb_wald", "welch"],
                                      format_func=lambda m: {"nb_wald": "NB-Wald (DESeq2 방식)",
                                                             "welch": "Welch t-test (log 정규화)"}[m])
        min_count = st.sidebar.number_input("최소 총 read 수 (필터)", 0, 1000, 10, 5)
        if no_rep:
            st.sidebar.warning("그룹당 반복이 1개뿐이라 고정 dispersion(BCV 0.1, edgeR 무반복 권장값)으로 "
                               "NB 검정을 수행합니다. 탐색적 결과로 해석하세요.")
        try:
            with st.spinner("DEG 분석 중..."):
                deg = run_deg(counts, meta, group_col, ref, test, method, min_count)
        except ValueError as e:
            st.error(str(e))
            st.stop()
        norm = A.normalize(counts[[c for c in counts.columns if str(c) in meta.index
                                   and meta.loc[str(c), group_col] in (ref, test)]])
        groups = meta.loc[norm.columns.astype(str), group_col]
else:
    fd = st.sidebar.file_uploader("DEG 결과표 (DESeq2/edgeR/limma)", type=["csv", "tsv", "txt", "xlsx"])
    if fd is not None:
        try:
            deg = A.standardize_deg(A.read_table(fd))
        except ValueError as e:
            st.error(str(e))
            st.stop()

st.title("SkinAge-Omics · 피부노화 전사체 분석")
st.caption("RNA-seq 기반 피부노화 관련 유전자 및 신호전달경로 분석 | DEG 시각화 · 피부노화 유전자 탐색 · GO/KEGG 경로 분석 · 기존 분석 비교")

if deg is None:
    st.info("왼쪽에서 데이터를 선택하거나 업로드하세요. 처음이라면 '데모 데이터'로 기능을 확인할 수 있습니다.")
    st.stop()

deg = A.classify(deg, padj_cut, lfc_cut)
n_up, n_down = (deg.status == "Up").sum(), (deg.status == "Down").sum()

c1, c2, c3, c4 = st.columns(4)
c1.metric("검정 유전자", f"{len(deg):,}")
c2.metric("증가 DEG (Up)", f"{n_up:,}")
c3.metric("감소 DEG (Down)", f"{n_down:,}")
panel = A.aging_panel(deg, padj_cut, lfc_cut)
c4.metric("피부노화 패널 유의 변화", f"{panel.status.isin(['Up', 'Down']).sum()} / {len(panel)}")

tab1, tab2, tab3, tab4, tab6, tab5 = st.tabs(["① DEG 시각화", "② 피부노화 유전자 탐색",
                                              "③ GO·KEGG 경로 분석", "④ 기존 분석과 비교",
                                              "⑤ AI 피부노화 점수", "ⓘ 정보"])

# ------------------------------------------------------------------ ① DEG 시각화
with tab1:
    a, b = st.columns(2)
    with a:
        hl = st.multiselect("Volcano에 표시할 유전자", sorted(deg.gene),
                            default=[g for g in ["MMP1", "COL1A1", "CDKN2A", "IL6"] if g in set(deg.gene)])
        f = P.volcano(deg, padj_cut, lfc_cut, label_genes=hl, title="Volcano plot")
        st.pyplot(f); fig_download(f, "volcano.png", "dl_vol")
    with b:
        if deg["basemean"].notna().any():
            f = P.ma_plot(deg, title="MA plot")
            st.pyplot(f); fig_download(f, "ma_plot.png", "dl_ma")
        else:
            st.info("baseMean 컬럼이 없어 MA plot을 생략합니다.")
    if norm is not None:
        st.subheader("상위 DEG 히트맵")
        topn = st.slider("표시할 상위 DEG 수", 10, 100, 40, 5)
        top = deg[deg.status != "NS"].nsmallest(topn, "padj")["gene"]
        idx_upper = {g.upper(): g for g in norm.index}
        nm = norm.rename(index=lambda g: str(g).upper())
        if len(top) >= 2:
            f = P.heatmap(nm, list(top), groups, title=f"상위 {len(top)}개 DEG (z-score, log2 정규화)")
            st.pyplot(f); fig_download(f, "heatmap_top_deg.png", "dl_hm")
    st.subheader("DEG 결과표")
    show = st.radio("표시", ["유의 DEG만", "전체"], horizontal=True)
    tbl = deg[deg.status != "NS"] if show == "유의 DEG만" else deg
    st.dataframe(tbl.sort_values("padj"), use_container_width=True, height=320)
    csv_download(deg, "deg_results.csv", "dl_deg", "전체 DEG 결과 CSV 저장")

# ------------------------------------------------------------------ ② 피부노화 유전자 탐색
with tab2:
    st.markdown("피부노화 관련 **10개 기능 범주, %d개 유전자** 큐레이션 패널과 비교합니다. "
                "'예상 방향'은 노화/광노화 피부에서 문헌상 보고되는 일반적 경향이며 참고용입니다."
                % len(panel))
    cats = st.multiselect("범주 선택", list(SKIN_AGING_PANEL), default=list(SKIN_AGING_PANEL))
    sub = panel[panel.category.isin(cats)]
    st.dataframe(A.panel_summary(sub), use_container_width=True, hide_index=True)
    a, b = st.columns([1.15, 1])
    with a:
        f = P.panel_plot(sub)
        st.pyplot(f); fig_download(f, "skin_aging_panel.png", "dl_panel")
    with b:
        st.dataframe(sub.rename(columns={"gene": "유전자", "category": "범주", "expected": "예상 방향",
                                         "status": "관측", "concordance": "방향 일치"}),
                     use_container_width=True, height=520, hide_index=True)
        csv_download(sub, "skin_aging_panel.csv", "dl_panel_csv")
    st.subheader("개별 유전자 검색")
    q = st.text_input("유전자 기호 (쉼표 구분, 예: MMP1, COL1A1, SIRT1)", "MMP1, COL1A1, SIRT1")
    qs = [s.strip().upper() for s in q.split(",") if s.strip()]
    hit = deg[deg.gene.isin(qs)]
    st.dataframe(hit, use_container_width=True, hide_index=True)
    miss = sorted(set(qs) - set(hit.gene))
    if miss:
        st.caption("데이터에 없는 유전자: " + ", ".join(miss))
    if norm is not None and len(hit):
        g = st.selectbox("그룹별 발현 보기", list(hit.gene))
        nm = norm.rename(index=lambda x: str(x).upper())
        if g in nm.index:
            vals = pd.DataFrame({"sample": nm.columns, "log2(정규화 count+1)": np.log2(nm.loc[g].values + 1),
                                 "group": groups.values})
            st.bar_chart(vals.set_index("sample")["log2(정규화 count+1)"])
            st.dataframe(vals.groupby("group")["log2(정규화 count+1)"].agg(["mean", "std"]).round(3))

# ------------------------------------------------------------------ ③ GO·KEGG
with tab3:
    a, b, c = st.columns(3)
    db = a.selectbox("유전자세트 DB", list(GENESETS) + ["사용자 GMT 업로드"])
    min_s = b.number_input("최소 세트 크기", 5, 100, 10)
    max_s = c.number_input("최대 세트 크기", 50, 2000, 500)
    if db == "사용자 GMT 업로드":
        up = st.file_uploader("GMT 파일 (예: MSigDB에서 내려받은 .gmt)", type=["gmt", "txt"])
        gsets = A.read_gmt(up) if up is not None else None
    else:
        gsets = load_gmt(str(GENESETS[db]))
    if gsets:
        universe = deg.gene
        res_up = A.ora(deg.gene[deg.status == "Up"], universe, gsets, min_s, max_s)
        res_dn = A.ora(deg.gene[deg.status == "Down"], universe, gsets, min_s, max_s)
        st.caption(f"배경 유전자 = 검정된 전체 {len(universe):,}개 | 초기하분포 검정 + BH-FDR")
        a, b = st.columns(2)
        with a:
            st.markdown(f"**증가 DEG ({n_up}) 과대표현 경로**")
            f = P.ora_dotplot(res_up, 15, f"Up · {db.split(' (')[0]}", "#C8463D")
            st.pyplot(f); fig_download(f, "ora_up.png", "dl_oraup")
            st.dataframe(res_up.head(50), use_container_width=True, height=260)
            csv_download(res_up, "ora_up.csv", "dl_oraup_csv")
        with b:
            st.markdown(f"**감소 DEG ({n_down}) 과대표현 경로**")
            f = P.ora_dotplot(res_dn, 15, f"Down · {db.split(' (')[0]}", "#2F6DB5")
            st.pyplot(f); fig_download(f, "ora_down.png", "dl_oradn")
            st.dataframe(res_dn.head(50), use_container_width=True, height=260)
            csv_download(res_dn, "ora_down.csv", "dl_oradn_csv")

# ------------------------------------------------------------------ ④ 기존 분석과 비교
with tab4:
    st.markdown("기존에 수행한 DEG 결과(DESeq2/edgeR/limma 등)를 올리면, 현재 결과와 "
                "**log2FC 상관, DEG 중복도, 방향 일치율**을 계산합니다.")
    src = ["파일 업로드"]
    if mode == EXO_MODE:
        src.insert(0, "기존 edgeR 결과 (2024.11, 동봉)")
    if counts is not None and not no_rep:
        src.append("앱 내 다른 방법(Welch ↔ NB-Wald)과 비교")
    how = st.radio("비교 대상", src, horizontal=True)
    other = None
    lab_b = "기존 분석"
    if how.startswith("기존 edgeR"):
        other = A.standardize_deg(pd.read_csv(EXO / "exosome_edgeR_2024-11.csv"))
        lab_b = "기존 edgeR"
    elif how == "파일 업로드":
        fo = st.file_uploader("기존 DEG 결과표", type=["csv", "tsv", "txt", "xlsx"], key="cmp")
        if fo is not None:
            other = A.standardize_deg(A.read_table(fo))
    else:
        alt = "welch" if method == "nb_wald" else "nb_wald"
        other = run_deg(counts, meta, group_col, ref, test, alt, min_count)
        lab_b = {"welch": "Welch t-test", "nb_wald": "NB-Wald"}[alt]
    if other is not None:
        r = A.compare_results(deg, other, padj_cut, lfc_cut, ("SkinAge-Omics", lab_b))
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("공통 유전자", f"{r['n_common']:,}")
        m2.metric("log2FC Pearson r", f"{r['pearson_lfc']:.3f}")
        m3.metric("DEG 중복 (Jaccard)", f"{r['jaccard']:.3f}")
        m4.metric("공통 DEG 방향 일치율", f"{r['direction_agreement']:.1%}")
        a, b = st.columns([1, 1])
        with a:
            f = P.compare_scatter(r["merged"], ("SkinAge-Omics", lab_b), r["pearson_lfc"])
            st.pyplot(f); fig_download(f, "compare_scatter.png", "dl_cmp")
        with b:
            st.markdown("**DEG 판정 교차표**")
            st.dataframe(r["crosstab"], use_container_width=True)
            st.markdown(f"- SkinAge-Omics DEG: **{r['deg_a']}** / {lab_b} DEG: **{r['deg_b']}** / 공통: **{r['deg_both']}**\n"
                        f"- {lab_b} DEG 재현율: **{r['recall_of_b']:.1%}**\n"
                        f"- Spearman ρ (log2FC): **{r['spearman_lfc']:.3f}**")
            csv_download(r["merged"].reset_index(), "comparison_merged.csv", "dl_cmp_csv")

# ------------------------------------------------------------------ ⑤ AI 피부노화 점수
with tab6:
    st.markdown("**머신러닝(L2 로지스틱 회귀)**으로 젊은/노화 피부 전사체를 구분하는 모델을 학습하고, "
                "샘플별 **전사체 피부노화 점수(0~1)**를 계산합니다. 성능은 Leave-one-out 교차검증으로 "
                "평가하며, 특징 선택(피부노화 패널 + 분산 상위 유전자)은 각 fold 안에서만 수행해 정보 누수를 막습니다.")
    if norm is None:
        st.info("AI 모듈은 count matrix 입력(데모 또는 업로드)에서 사용할 수 있습니다.")
    elif groups.value_counts().min() < 3:
        st.info("AI 분류 모델은 그룹마다 3개 이상의 샘플이 필요합니다. 현재 데이터는 샘플 수가 부족해 "
                "모델을 학습하지 않습니다. 데모 데이터에서 기능을 확인할 수 있습니다.")
    else:
        a, b = st.columns(2)
        n_var = a.slider("분산 상위 유전자 수 (+ 피부노화 패널)", 50, 1000, 200, 50)
        lam = b.select_slider("L2 정규화 강도 λ", [0.1, 1.0, 10.0, 100.0], value=10.0)
        with st.spinner("모델 학습·교차검증 중..."):
            scores, imp, met, model = ML.aging_classifier(norm, groups, test, n_var, lam)
            pcs, var = ML.pca(norm)
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("샘플 수", met["n_samples"])
        m2.metric("특징 유전자 수", met["n_features"])
        m3.metric("LOOCV 정확도", f"{met['loocv_accuracy']:.1%}")
        m4.metric("LOOCV AUC", f"{met['loocv_auc']:.3f}")
        if met["n_samples"] < 10:
            st.warning("샘플 수가 적어 성능 추정의 불확실성이 큽니다. 독립 코호트 검증을 권장합니다.")
        a, b = st.columns(2)
        with a:
            f = P.pca_plot(pcs, var, groups)
            st.pyplot(f); fig_download(f, "pca.png", "dl_pca")
        with b:
            f = P.score_plot(scores)
            st.pyplot(f); fig_download(f, "ai_aging_score.png", "dl_score")
        a, b = st.columns([1, 1])
        with a:
            f = P.importance_plot(imp, 20)
            st.pyplot(f); fig_download(f, "ai_feature_importance.png", "dl_imp")
        with b:
            st.dataframe(scores, use_container_width=True, hide_index=True)
            csv_download(scores, "ai_aging_scores.csv", "dl_score_csv")
            csv_download(imp, "ai_feature_importance.csv", "dl_imp_csv", "특징 중요도 CSV 저장")
        st.subheader("새 샘플 점수 예측")
        fn = st.file_uploader("새 샘플 count matrix (같은 형식)", type=["csv", "tsv", "txt", "xlsx"], key="newscore")
        if fn is not None:
            new = A.read_table(fn)
            new = new.set_index(new.columns[0]).apply(pd.to_numeric, errors="coerce").fillna(0)
            sc = ML.score_new(model, A.normalize(new)).round(3).rename("전사체 피부노화 점수")
            st.dataframe(sc.reset_index().rename(columns={"index": "sample"}), hide_index=True)

# ------------------------------------------------------------------ 정보
with tab5:
    st.markdown(f"""
**SkinAge-Omics v{VERSION}** · 개발자: 문경미

| 단계 | 방법 |
|---|---|
| 정규화 | DESeq2 median-of-ratios size factor |
| DEG (NB-Wald) | 음이항 GLM 최대우도(그룹 2수준) + Wald 검정, 평균-분산 추세로 dispersion 수축 |
| DEG (Welch) | log2(정규화 count+1)에 Welch t-test |
| 다중검정 | Benjamini-Hochberg FDR |
| AI 모듈 | PCA(비지도), L2 로지스틱 회귀(지도) + LOOCV 성능 평가, 표준화 계수 기반 중요 유전자 |
| 경로 분석 | 초기하분포 ORA (배경 = 검정 유전자 전체), Up/Down 분리 |
| 유전자세트 | KEGG (MSigDB c2.cp.kegg v7.5.1), GO BP 2021 (Enrichr), MSigDB Hallmark v7.0 — GSEApy 저장소(BSD-3) 동봉본 |

입력 형식: count matrix(csv/tsv/xlsx, 첫 열 유전자 기호), 샘플 정보(첫 열 샘플명 + 그룹 컬럼),
또는 DESeq2(`log2FoldChange, pvalue, padj, baseMean`)·edgeR(`logFC, PValue, FDR, logCPM`)·limma(`logFC, P.Value, adj.P.Val, AveExpr`) 결과표.
""")
