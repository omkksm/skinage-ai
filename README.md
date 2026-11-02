# SkinAge-AI · 피부질환 예방효과 정량예측 AI 플랫폼

피부노화·아토피 피부염 개체(세포·마우스·인체)에서 **환경·스트레스와 후보물질 처리**에 따른
**기전지표 · 질환 지수 · 예방효과**를 AI로 정량 예측하는 플랫폼입니다.

**실행 주소 (설치 불필요)**: https://omkksm.github.io/skinage-ai/

개발: 문경미 · 2026.11 · v3.1

## 모듈
| 모듈 | 역할 | 주요 기능 |
|---|---|---|
| AI 예측 (핵심) | 예측 엔진 | 16개 기전지표 예측(95% 예측구간), 질환 지수·예방률, 용량-반응·권장 농도, 후보물질 순위, 교차검증·중요도, 실험·문헌(PMID) 데이터 재학습, 전사체 기전 역전 점수 |
| SkinAge-Omics | 기전 발굴 | RNA-seq DEG(NB-Wald), 피부노화 유전자 패널 70개, GO·KEGG·Hallmark ORA, 기존 분석 비교 |
| SkinCell ProtoLab | Wet-lab | 피부세포 노화모델 15단계 프로토콜, 계산기 10종, 실험실 구축·운영 |
| PigmentLab | 효능 정량 | 멜라닌·티로시나제·IC50(4PL), k-means 이미지 색소 분석 |

## 구성
- `index.html`, `js/`, `data/`, `modules/` : 웹판 (HTML·JavaScript, 브라우저에서 바로 실행)
- `research/` : 연구판 (Python·Streamlit, 동일 알고리즘, 자동 테스트 포함) — `streamlit run streamlit_app.py`

## 검증
- AI 예측: 기전 시뮬레이션 데이터 5-fold 교차검증 평균 R² 0.79, 95% 예측구간 포함률 90%
- 전사체: 합성 데이터 DEG 정밀도 98.5%·위양성률 4.2%, 연구실 RNA-seq 실데이터에서 기존 edgeR 결과와 순위상관 1.000·기존 DEG 89.7% 재현·방향 100% 일치
- 한계: 예측 모델은 기전 가정 기반 시뮬레이션 데이터로 학습된 방법론 검증 단계이며, 후보물질 A~E는 시연용 가상 물질입니다.

## 데이터
- 유전자세트: MSigDB KEGG v7.5.1·Hallmark v7.0, Enrichr GO BP 2021 (GSEApy 저장소 동봉본)
- 연구실 실데이터 원자료는 포함하지 않으며, 요약 결과(`data/exosome_summary.json`)만 포함합니다.
