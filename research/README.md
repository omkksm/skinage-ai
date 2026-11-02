# SkinAge-AI — 피부노화·아토피 피부염 예방효과 정량예측 AI 플랫폼

개발자: 문경미 · 버전 3.1.0 · Python / Streamlit

하나의 앱에서 네 모듈을 사용합니다: ⓪ SkinAge-AI 예측(핵심) ① SkinAge-Omics(전사체 분석) ② SkinCell ProtoLab(Wet-lab 가이드) ③ PigmentLab(색소 분석).
웹 실행 주소(설치 불필요): https://claude.ai/artifact/T9vMeKuPErX6ABjuTni2Vg

## 주요 기능
1. **DEG 분석·시각화**: count matrix → DESeq2 방식 정규화 → NB-Wald 또는 Welch 검정 → Volcano, MA, 히트맵, 결과표
   (DESeq2/edgeR/limma 결과표를 바로 올려 시각화도 가능)
2. **피부노화 유전자 탐색**: 10개 범주, 70개 유전자 큐레이션 패널(MMP, 콜라겐, 세포노화, SASP, 항산화, 장벽, 장수신호 등)과 방향 일치도
3. **GO·KEGG 경로 분석**: 초기하분포 ORA + BH-FDR (KEGG, GO BP 2021, Hallmark 내장 / 사용자 GMT 업로드)
4. **기존 분석과 비교**: 기존 DEG 결과표와 log2FC 상관, DEG 중복도(Jaccard), 방향 일치율, 교차표
5. **AI 피부노화 점수**: PCA + L2 로지스틱 회귀 분류기, Leave-one-out 교차검증, 중요 유전자, 새 샘플 점수 예측

## 실행
```bash
pip install -r requirements.txt
streamlit run streamlit_app.py
```
(SkinAge-Omics 단독 실행: streamlit run app.py)
Windows: `run_app_windows.bat` 더블클릭.

## 명령줄 일괄 분석 (앱과 동일한 엔진)
```bash
python scripts/run_pipeline.py --counts counts.csv --meta meta.csv --group group --ref young --test aged --compare 기존_DESeq2.csv
```
→ `results/figures`, `results/tables`, `results/summary.json`

## 검증
`python tests/test_all.py` — 18개 자동 테스트(정규화, BH 보정, 정답을 아는 시뮬레이션 데이터에서 DEG 정밀도와 위양성률, 입력 형식 인식, 경로 검출, 비교 모듈, AI 교차검증과 라벨 셔플 대조, 앱 전체 실행).

## 폴더
```
app.py                  Streamlit 앱
core/analysis.py        정규화·DEG·패널·ORA·비교
core/ml.py              AI 모듈 (PCA, 로지스틱 회귀, LOOCV)
core/plots.py           시각화
core/skin_aging_genes.py 피부노화 유전자 패널
data/genesets/          KEGG / GO BP / Hallmark GMT
data/demo/              합성 데모 데이터 (젊은 6 vs 노화 6, 정답 포함)
scripts/                데모 생성, 일괄 분석, 화면 캡처
tests/                  자동 테스트
```

## 데이터 출처·라이선스
- 유전자세트: GSEApy 저장소(BSD-3-Clause) 동봉본 — MSigDB c2.cp.kegg v7.5.1, MSigDB Hallmark v7.0, Enrichr GO_Biological_Process_2021.
  KEGG/MSigDB 데이터는 각 제공처의 이용 조건을 따릅니다(학술·비영리 연구용).
- 데모 데이터는 합성 데이터이며 실제 사람 데이터가 아닙니다.

## 원자료 보호
연구실 엑소좀 RNA-seq 원자료(count·edgeR 전체 결과)는 배포본에 포함하지 않았습니다. 기전 역전 예시는 미리 계산한 요약(data/exosome_summary.json: 유의 DEG·패널 유전자·요약 통계)만 사용합니다.
