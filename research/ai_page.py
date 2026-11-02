"""SkinAge-AI 예측 모듈 (통합 앱의 첫 화면).

피부노화·아토피 피부염 개체(세포·마우스·인체)에서 환경·스트레스와 후보물질 처리에 따른
기전지표, 질환 지수, 예방효과를 정량 예측한다.
"""
from pathlib import Path
import hashlib
import io

import numpy as np
import pandas as pd
import streamlit as st

from core import predict as PR
from core import plots as P
from core import analysis as A
from core.skin_aging_genes import panel_table

ROOT = Path(__file__).resolve().parent


@st.cache_data(show_spinner=False)
def sim_data(n=1800, seed=7):
    return PR.simulate(n, seed)


@st.cache_resource(show_spinner=False)
def train_model(key, data):
    return PR.SkinAgeAI(data, max_n=400)


@st.cache_data(show_spinner=False)
def run_cv(key, data):
    return PR.cross_validate(data, k=5, max_n=300)


def data_key(df):
    return hashlib.md5(pd.util.hash_pandas_object(df[["organism", "stress", "level", "days", "compound", "conc",
                                                      "timing", "biomarker", "value"]], index=False).values).hexdigest()


def show(fig, name, key):
    st.pyplot(fig)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    st.download_button("그림 PNG 저장", buf.getvalue(), file_name=name, mime="image/png", key=key)
    P.plt.close(fig)


# ------------------------------------------------------------------ 학습 데이터 구성
sim = sim_data()
extra = st.session_state.get("ai_extra_data")
data = pd.concat([sim, extra], ignore_index=True) if extra is not None and len(extra) else sim
key = data_key(data)
with st.spinner("AI 예측 모델 학습 중… (처음 한 번, 브라우저 실행 시 1분 정도)"):
    ai = train_model(key, data)

st.title("SkinAge-AI · 피부 질환 예방효과 정량예측 플랫폼")
st.caption("피부노화·아토피 피부염 개체(세포·마우스·인체)에서 환경·스트레스와 후보물질 처리에 따른 "
           "기전지표 · 질환 지수 · 예방효과를 AI로 정량 예측합니다.")

c1, c2, c3, c4 = st.columns(4)
c1.metric("예측 지표", f"{len(ai.models)}개")
c2.metric("학습 측정값", f"{len(data):,}건")
c3.metric("실험·문헌 데이터", f"{int((data.source != '시뮬레이션').sum()):,}건")
c4.metric("후보물질", f"{len(ai.compounds)}종")

tabs = st.tabs(["① 예측 시뮬레이터", "② 용량-반응 · 최적 농도", "③ 후보물질 비교", "④ 모델 성능 · 근거",
                "⑤ 데이터 허브 (실험·문헌)", "⑥ 전사체 기전 역전"])

# ------------------------------------------------------------------ 공통 조건 선택
with st.sidebar:
    st.subheader("예측 조건")
    disease = st.radio("예측할 질환", list(PR.DISEASES), format_func=lambda d: PR.DISEASES[d]["label"],
                       horizontal=True, key="ai_disease")
    org = st.selectbox("개체", list(PR.ORGANISMS), format_func=lambda o: PR.ORGANISMS[o]["label"], key="ai_org")
    stress = st.selectbox("환경·스트레스 요인", PR.DISEASES[disease]["stresses"],
                          format_func=lambda s: PR.STRESSES[s]["label"], key="ai_stress")
    st.caption("강도 기준: " + PR.STRESSES[stress]["unit"][org])
    level = st.select_slider("스트레스 강도", [0.25, 0.5, 0.75, 1.0], value=0.75,
                             format_func=lambda v: {0.25: "약", 0.5: "중", 0.75: "강", 1.0: "최강"}[v], key="ai_level")
    d0, d1 = PR.ORGANISMS[org]["dur"]
    days = st.slider("노출 기간 (일)", d0, d1, int((d0 + d1) / 2), key=f"ai_days_{org}")
    compound = st.selectbox("후보물질", ai.compounds, key="ai_comp")
    unit = PR.ORGANISMS[org]["conc_unit"]
    cmax = PR.ORGANISMS[org]["conc_max"]
    conc = st.slider(f"농도 ({unit})", 0.0, float(cmax), float(cmax * 0.3), float(cmax / 100), key=f"ai_conc_{org}")
    timing = st.radio("처리 시점", list(PR.TIMING), format_func=lambda t: PR.TIMING[t][0], key="ai_timing")

# ------------------------------------------------------------------ ① 시뮬레이터
with tabs[0]:
    prof = PR.prevention_profile(ai, org, stress, level, days, compound, conc, timing, disease)
    dz = PR.disease_prediction(prof, org, disease)
    score = PR.composite_score(prof, disease)
    lab = PR.DISEASES[disease]["short"]
    m1, m2, m3, m4 = st.columns(4)
    m1.metric(f"{lab} (스트레스 단독)", f"{dz['stress_only']:.0f} / 100")
    m2.metric(f"{lab} (+ 후보물질)", f"{dz['with_compound']:.0f} / 100",
              delta=f"{dz['with_compound'] - dz['stress_only']:.0f}", delta_color="inverse")
    m3.metric("질환 예방률", f"{dz['prevention']:.0f}%")
    m4.metric("기전 축 가중 예방점수", f"{score:.0f}%")
    st.caption(f"{PR.ORGANISMS[org]['label']} · {PR.STRESSES[stress]['label']} {level} · {days}일 · "
               f"{compound} {conc:.3g} {unit} · {PR.TIMING[timing][0]}")
    a, b = st.columns([1.35, 1])
    with a:
        show(P.prevention_bars(prof, f"기전지표 예측 ({PR.ORGANISMS[org]['label']})"), "ai_prevention.png", "ai_dl1")
    with b:
        if prof.axis.nunique() >= 3:
            show(P.radar(prof), "ai_radar.png", "ai_dl2")
        else:
            st.info("기전 축이 3개 미만이라 레이더 차트를 생략합니다.")
    tbl = prof.rename(columns={"biomarker": "지표", "axis": "기전 축", "stress_only": "스트레스 단독(%)",
                               "with_compound": "후보물질 병용(%)", "ci95": "±95% 구간", "prevention": "예방률(%)"})
    st.dataframe(tbl.drop(columns="axis_key").round(1), use_container_width=True, hide_index=True)
    st.caption("값은 무처리 대조 대비 % 변화. 예방률 = (스트레스 단독 − 병용) ÷ 스트레스 단독 × 100. "
               "질환 지수는 질환 관련 기전 축의 손상 정도를 가중 평균해 0~100으로 환산한 값입니다. "
               "스트레스 영향이 10% 미만인 지표는 예방률을 계산하지 않습니다.")

# ------------------------------------------------------------------ ② 용량-반응
with tabs[1]:
    st.markdown(f"**{compound}** · {PR.ORGANISMS[org]['label']} · {PR.STRESSES[stress]['label']} · "
                f"{PR.TIMING[timing][0]} 조건에서 농도별 **{PR.DISEASES[disease]['label']} 예방점수**를 예측합니다.")
    with st.spinner("농도별 예측 중…"):
        dr = PR.dose_response(ai, org, stress, level, days, compound, timing, n=12, disease=disease)
        fit = PR.fit_4pl(dr.conc, dr.score)
    k1, k2, k3 = st.columns(3)
    k1.metric("EC₅₀", f"{fit['ec50']:.3g} {unit}")
    k2.metric("최대 예방점수(외삽)", f"{fit['top']:.0f}%")
    k3.metric("권장 농도 (최대효과 80%)", f"{fit['rec80']:.3g} {unit}",
              help="최대 예방효과의 80%에 도달하는 최소 농도. 시험 범위를 넘으면 참고값입니다.")
    a, b = st.columns([1.2, 1])
    with a:
        show(P.dose_curve(dr, fit, unit, f"용량-반응 · {compound}"), "ai_dose.png", "ai_dl3")
    with b:
        st.dataframe(dr.rename(columns={"conc": f"농도({unit})", "score": "예방점수(%)"}).round(3),
                     use_container_width=True, hide_index=True, height=420)
    if fit["rec80"] > cmax:
        st.warning("권장 농도가 시험 범위를 넘습니다. 이 조건에서는 농도를 높여도 효과가 제한적일 수 있습니다.")

# ------------------------------------------------------------------ ③ 후보물질 비교
with tabs[2]:
    st.markdown("현재 조건(개체·스트레스·강도·기간·농도)에서 모든 후보물질과 처리 시점을 비교합니다.")
    rows = []
    for c in ai.compounds:
        for t in PR.TIMING:
            p = PR.prevention_profile(ai, org, stress, level, days, c, conc, t, disease)
            dzc = PR.disease_prediction(p, org, disease)
            rows.append({"후보물질": c, "처리 시점": PR.TIMING[t][0], "예방점수": PR.composite_score(p, disease),
                         "질환 예방률": dzc["prevention"], f"{lab}(병용)": dzc["with_compound"]})
    rk = pd.DataFrame(rows)
    best = rk.loc[rk.groupby("후보물질")["예방점수"].idxmax()]
    a, b = st.columns([1, 1.1])
    with a:
        show(P.ranking_bars(best[["후보물질", "예방점수"]], f"{PR.DISEASES[disease]['label']} 예방 · 물질별 최적 시점"),
             "ai_ranking.png", "ai_dl4")
    with b:
        piv = rk.pivot(index="후보물질", columns="처리 시점", values="예방점수").round(1)
        st.dataframe(piv, use_container_width=True)
        top = best.sort_values("예방점수", ascending=False).iloc[0]
        st.success(f"이 조건의 1순위: **{top['후보물질']}** ({top['처리 시점']}), 예방점수 {top['예방점수']:.0f}%")

# ------------------------------------------------------------------ ④ 모델 성능
with tabs[3]:
    st.markdown("**모델**: 지표별 가우시안 프로세스 회귀(RBF 커널, 주변우도 최대화로 하이퍼파라미터 선택). "
                "입력 = 개체 · 스트레스 종류×강도 · 노출 기간 · 후보물질×log 농도 · 처리 시점. "
                "출력 = 지표 변화(%)와 95% 예측구간.")
    if st.button("5-fold 교차검증 실행", key="ai_cv_btn") or st.session_state.get("ai_cv_done") == key:
        with st.spinner("교차검증 중… (브라우저 실행 시 1~2분)"):
            cv, pairs = run_cv(key, data)
        st.session_state["ai_cv_done"] = key
        show_cv = cv.rename(columns={"biomarker": "지표", "R2": "R²", "PI95_coverage": "95% 구간 포함률",
                                     "R2_vs_truth": "R² (정답 대비)"})
        st.dataframe(show_cv.round(3), use_container_width=True, hide_index=True)
        st.caption(f"지표 {len(cv)}개 평균 R² {cv.R2.mean():.3f} · 95% 예측구간 포함률 평균 {cv.PI95_coverage.mean():.1%}. "
                   "'정답 대비 R²'는 시뮬레이션 데이터의 노이즈 없는 정답값과 비교한 값입니다.")
        show(P.parity_plot(pairs), "ai_cv_parity.png", "ai_dl5")
    else:
        st.info("버튼을 누르면 교차검증으로 지표별 예측 정확도를 계산합니다.")
    st.subheader("예측에 영향을 준 요인")
    bsel = st.selectbox("지표", sorted(ai.models), key="ai_imp_bm")
    show(P.importance_group_bars(PR.permutation_importance(ai, bsel, n_rep=2), f"{bsel} · 순열 중요도"),
         "ai_importance.png", "ai_dl6")

# ------------------------------------------------------------------ ⑤ 데이터 허브
with tabs[4]:
    st.markdown("실험값이나 **논문(PubMed/PMC)에서 추출한 값**을 표준 양식으로 올리면 시뮬레이션 데이터와 합쳐 "
                "모델을 다시 학습합니다. `pmid` 칸에 번호를 넣으면 '문헌', 비워 두면 '실데이터'로 구분됩니다.")
    st.download_button("입력 양식 CSV 내려받기", PR.template_df().to_csv(index=False).encode("utf-8-sig"),
                       "SkinAge-AI_입력양식.csv", "text/csv", key="ai_tpl")
    with st.expander("양식 설명"):
        st.markdown(f"""
- **organism**: cell / mouse / human · **stress**: uv / ox / gly / allergen / barrier
- **level**: 스트레스 강도 0.25(약)·0.5(중)·0.75(강)·1.0(최강) · **days**: 노출 기간(일)
- **compound**: 물질명 (스트레스 단독군은 `없음`) · **conc**: 세포 µM, 마우스·인체 도포 % · **timing**: pre / co / post
- **biomarker**: {", ".join(PR.BM.biomarker)}
- **value**: 무처리 대조 대비 % 변화 (증가 +, 감소 −) · **pmid**: 문헌 출처 번호(선택)
""")
    up = st.file_uploader("실험·문헌 데이터 업로드 (csv/xlsx)", type=["csv", "xlsx", "txt", "tsv"], key="ai_up")
    if up is not None:
        try:
            real = PR.standardize_real(A.read_table(up))
            st.session_state["ai_extra_data"] = real
            st.success(f"{len(real)}건을 반영했습니다 (실데이터 {int((real.source == '실데이터').sum())}, "
                       f"문헌 {int((real.source == '문헌').sum())}). 모델이 다시 학습됩니다.")
            st.rerun()
        except ValueError as e:
            st.error(str(e))
    if extra is not None and len(extra):
        st.dataframe(extra, use_container_width=True, height=240)
        if st.button("업로드 데이터 제거", key="ai_clear"):
            st.session_state.pop("ai_extra_data", None)
            st.rerun()
    summ = data.groupby(["source", "organism"]).size().unstack(fill_value=0)
    summ.columns = [PR.ORGANISMS[c]["label"] for c in summ.columns]
    st.markdown("**현재 학습 데이터 구성 (측정값 수)**")
    st.dataframe(summ, use_container_width=True)
    st.caption("시뮬레이션 데이터: 스트레스 강도·기간에 대한 포화형 손상, 개체 민감도, 후보물질 Hill 용량-반응과 "
               "처리 시점 효과를 반영한 기전 모델에서 생성했습니다(정답이 알려진 검증용). 후보물질 A~E는 시연용 가상 물질입니다.")

# ------------------------------------------------------------------ ⑥ 전사체 기전 역전
with tabs[5]:
    st.markdown("후보물질 처리 전사체(DEG)가 **피부노화 유전자 패널의 노화 방향을 얼마나 되돌리는지** 계산합니다 "
                "(Connectivity Map 방식의 역전 점수, −1~+1, 양수 = 노화 역전).")
    srcsel = st.radio("DEG 데이터", ["실데이터 예시: 엑소좀 EV17 (마우스 색소세포)", "DEG 결과표 업로드"], key="ai_rev_src")
    deg = None
    if srcsel.startswith("실데이터"):
        raw = ROOT / "data" / "exosome" / "exosome_edgeR_2024-11.csv"
        if raw.exists():
            deg = A.standardize_deg(pd.read_csv(raw))
        else:  # 원자료 보호 배포본: 미리 계산한 요약(유의 DEG·패널 유전자)만 사용
            import json
            sm = json.loads((ROOT / "data" / "exosome_summary.json").read_text(encoding="utf-8"))
            deg = pd.DataFrame(sm["genes"], columns=["gene", "log2fc", "padj", "status"])
    else:
        f = st.file_uploader("DEG 결과표 (DESeq2/edgeR/limma)", type=["csv", "tsv", "txt", "xlsx"], key="ai_rev_up")
        if f is not None:
            deg = A.standardize_deg(A.read_table(f))
    if deg is not None:
        sc, used = PR.reversal_score(deg, panel_table())
        r1, r2 = st.columns(2)
        r1.metric("피부노화 역전 점수", f"{sc:+.3f}")
        r2.metric("사용된 패널 유전자", f"{len(used)}개")
        if len(used):
            cat = used.groupby("category")["reversal"].mean().sort_values()
            fig, ax = P.plt.subplots(figsize=(6.4, 3.4))
            ax.barh(cat.index, cat.values, color=["#2F8F6A" if v > 0 else "#C8463D" for v in cat.values])
            ax.axvline(0, c="#555", lw=0.8)
            ax.set_xlabel("범주별 역전 점수 (+ 노화 역전 / − 노화 방향)")
            ax.set_title("기전 범주별 역전 점수")
            fig.tight_layout()
            show(fig, "ai_reversal.png", "ai_dl7")
            st.dataframe(used.round(3), use_container_width=True, hide_index=True, height=300)
        st.caption("마우스 유전자 기호는 대문자로 사람 기호에 대응했습니다. 역전 점수는 기전 가설 생성을 위한 탐색 지표입니다.")
