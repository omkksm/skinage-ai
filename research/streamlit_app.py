"""SkinAge-AI — 피부노화·아토피 피부염 예방효과 정량예측 AI 플랫폼

하나의 앱에서 네 모듈을 사용합니다.
  0. SkinAge-AI 예측  환경·스트레스 × 후보물질 → 기전지표·질환 지수·예방효과 정량예측 (핵심)
  1. SkinAge-Omics    RNA-seq 피부노화 전사체 분석 (DEG · 피부노화 유전자 · GO/KEGG · 기존 분석 비교 · AI 점수)
  2. SkinCell ProtoLab 피부세포 노화연구 Wet-lab 프로토콜 · 계산기 · 실험실 구축/운영
  3. PigmentLab        멜라닌 · 티로시나제 · AI 이미지 색소 분석 (기타)

실행:  streamlit run streamlit_app.py
개발:  문경미
"""
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

ROOT = Path(__file__).resolve().parent
VERSION = "3.0.0"

st.set_page_config(page_title="SkinAge-AI", page_icon="🧬", layout="wide")

MODULES = {
    "🤖 SkinAge-AI · 예방효과 예측": "ai",
    "🧬 SkinAge-Omics · 전사체 분석": "omics",
    "🧫 SkinCell ProtoLab · Wet-lab": "protolab",
    "🎨 PigmentLab · 색소 분석": "pigment",
    "ⓘ 플랫폼 소개": "about",
}

st.sidebar.title("SkinAge-AI")
st.sidebar.caption(f"피부노화·아토피 예방효과 정량예측 AI 플랫폼 v{VERSION} · 개발 문경미")
choice = st.sidebar.radio("모듈 선택", list(MODULES), key="module")
st.sidebar.divider()
page = MODULES[choice]


def html_module(fname, height):
    html = (ROOT / "modules" / fname).read_text(encoding="utf-8")
    components.html(html, height=height, scrolling=True)


def run_script(name):
    code = (ROOT / name).read_text(encoding="utf-8")
    g = {"__name__": "__suite__", "__file__": str(ROOT / name), "_SUITE": True}
    exec(compile(code, str(ROOT / name), "exec"), g)


if page == "ai":
    run_script("ai_page.py")

elif page == "omics":
    run_script("app.py")

elif page == "protolab":
    st.title("SkinCell ProtoLab · 피부세포 노화연구 Wet-lab 가이드")
    st.caption("배양 → 노화모델 유도 → 평가 15단계 프로토콜 · 실험 계산기 10종 · 실험실 구축·운영 플래너")
    html_module("skincell_protolab.html", 1500)

elif page == "pigment":
    st.title("PigmentLab · 멜라닌 생성 · 미백 효능 분석")
    st.caption("멜라닌 표준곡선 정량 · 티로시나제 동역학과 IC₅₀ · k-means 기반 AI 이미지 색소 분석")
    html_module("pigmentlab.html", 1500)

else:
    st.title("SkinAge-AI · 피부 질환 예방효과 정량예측 AI 플랫폼")
    st.markdown(f"""
**목표**: 피부노화·아토피 피부염 개체(세포·마우스·인체)에서 **환경·스트레스와 후보물질 처리**에 따른
**기전지표 · 질환 지수 · 예방효과**를 AI로 정량 예측합니다.

| 모듈 | 역할 | 주요 기능 |
|---|---|---|
| 🤖 SkinAge-AI 예측 | **예측 엔진 (핵심)** | 조건 입력 → 16개 기전지표 변화와 95% 예측구간, 질환 지수, 예방률, 용량-반응·권장 농도, 후보물질 순위, 실험·문헌(PMID) 데이터 반영 재학습, 전사체 기전 역전 점수 |
| 🧬 SkinAge-Omics | 기전 발굴 (Dry-lab) | RNA-seq DEG, 피부노화 유전자 패널, GO·KEGG, 기존 분석 비교, AI 노화 점수 |
| 🧫 SkinCell ProtoLab | 데이터 생산 (Wet-lab) | 세포 노화모델 15단계 프로토콜, 계산기 10종, 실험실 구축·운영 |
| 🎨 PigmentLab | 효능 정량 | 멜라닌·티로시나제·IC₅₀, AI 이미지 색소 분석 |

**순환 구조**: ProtoLab·PigmentLab에서 표준화된 실험 데이터 생산 → 데이터 허브로 학습 → SkinAge-AI가 다음 실험 조건(물질·농도·시점) 추천 → SkinAge-Omics로 기전 확인

**질환 · 요인 · 지표**
- 피부노화: 광노화(UV), 산화스트레스, 당화·내인성 노화 → MMP-1, 콜라겐, p16, SA-β-gal, IL-6, ROS, 멜라닌, 주름, 탄력, 표피 두께
- 아토피 피부염: 알레르겐·햅텐(DNCB·HDM·MC903), 장벽 손상, 산화스트레스 → TSLP, IL-4·IL-13, IgE, Filaggrin, TEWL, 가려움

버전 {VERSION} · 개발 문경미 · 2026.11
""")
