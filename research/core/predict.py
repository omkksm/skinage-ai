"""SkinAge-AI: 환경·스트레스 × 후보물질 처리에 따른 피부노화·아토피 피부염 기전지표·질환 지수·예방효과 정량예측 엔진.

구성
- 기전 축(ECM 분해, 콜라겐, 세포노화, 염증, 산화, 색소, 장벽, 형태)과 개체별(세포·마우스·인체) 측정 지표 정의
- 기전 기반 시뮬레이터: 스트레스 강도·기간(포화형) × 개체 민감도 × 후보물질 보호효과(Hill 용량-반응 × 처리시점)
  → 지표 변화(% 무처리 대조 대비). 정답을 아는 학습·검증 데이터를 만든다.
- 예측 모델: 지표별 가우시안 프로세스 회귀(RBF, 주변우도 기반 하이퍼파라미터 선택) → 예측값 + 95% 예측구간
- 검증: K-fold 교차검증 R²·RMSE·예측구간 포함률, 순열 중요도
- 예방효과: 지표별 예방률 = (스트레스 단독 변화 − 병용 변화) / 스트레스 단독 변화 × 100, 기전 축 가중 종합점수
- 용량-반응: 농도 스윕 → 4PL 적합 → EC50, 최대 예방률의 80% 도달 최소 농도(권장 농도)
- 실데이터 병합: 표준 양식(long format) 업로드 → 시뮬레이션 데이터와 합쳐 재학습
numpy 만 사용.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- 정의
ORGANISMS = {
    "cell": {"label": "세포 (in vitro)", "sens": 1.0, "deliver": 1.0, "dur": (1, 7),
             "conc_unit": "µM", "conc_max": 100.0, "ec50_scale": 1.0},
    "mouse": {"label": "마우스 (HR-1, in vivo)", "sens": 0.75, "deliver": 0.75, "dur": (14, 84),
              "conc_unit": "% (도포)", "conc_max": 2.0, "ec50_scale": 0.02},
    "human": {"label": "인체 (적용시험)", "sens": 0.45, "deliver": 0.55, "dur": (28, 168),
              "conc_unit": "% (도포)", "conc_max": 2.0, "ec50_scale": 0.025},
}
# 질환(목표 상태): 질환 지수는 관련 기전 축의 손상 정도를 가중 합산해 0~100으로 환산
DISEASES = {
    "aging": {"label": "피부노화", "short": "노화 지수",
              "axes": {"ECM": 1.0, "COL": 1.2, "SEN": 1.0, "INF": 0.6, "OX": 0.8, "PIG": 0.6, "BAR": 0.5, "MOR": 1.2},
              "stresses": ["uv", "ox", "gly"]},
    "atopic": {"label": "아토피 피부염", "short": "아토피 중증도 지수",
               "axes": {"TH2": 1.2, "IGE": 0.8, "ITCH": 1.0, "BAR": 1.2, "INF": 0.9, "MOR": 0.5, "OX": 0.4},
               "stresses": ["allergen", "barrier", "ox"]},
}
STRESSES = {
    "uv": {"label": "광노화 (UVB/UVA)",
           "unit": {"cell": "UVB 10–50 mJ/cm²", "mouse": "UVB 60–180 mJ/cm², 주 3회", "human": "일광 누적노출 지수"},
           "w": {"ECM": 1.0, "COL": 0.8, "SEN": 0.5, "INF": 0.9, "OX": 0.8, "PIG": 1.0, "BAR": 0.7, "MOR": 0.9,
                 "TH2": 0.15, "IGE": 0.05, "ITCH": 0.2}},
    "ox": {"label": "산화스트레스 (H₂O₂ 등)",
           "unit": {"cell": "H₂O₂ 50–400 µM", "mouse": "산화 유도 모델 강도", "human": "산화스트레스 노출 지수"},
           "w": {"ECM": 0.6, "COL": 0.5, "SEN": 1.0, "INF": 0.6, "OX": 1.0, "PIG": 0.5, "BAR": 0.4, "MOR": 0.5,
                 "TH2": 0.25, "IGE": 0.1, "ITCH": 0.25}},
    "gly": {"label": "당화(AGEs)·내인성 노화",
            "unit": {"cell": "AGEs/고당 또는 반복계대", "mouse": "D-gal·고당 식이·연령", "human": "연령·당화 지수"},
            "w": {"ECM": 0.4, "COL": 1.0, "SEN": 0.9, "INF": 0.5, "OX": 0.6, "PIG": 0.3, "BAR": 0.5, "MOR": 0.8,
                  "TH2": 0.05, "IGE": 0.0, "ITCH": 0.1}},
    "allergen": {"label": "알레르겐·햅텐 감작 (DNCB·집먼지진드기·MC903)",
                 "unit": {"cell": "TNF-α/IFN-γ 또는 IL-4/IL-13 자극 (HaCaT)", "mouse": "DNCB 반복 도포 (BALB/c·NC/Nga)",
                          "human": "알레르겐 노출·증상 지수"},
                 "w": {"ECM": 0.3, "COL": 0.2, "SEN": 0.2, "INF": 0.8, "OX": 0.4, "PIG": 0.2, "BAR": 0.8, "MOR": 0.7,
                       "TH2": 1.0, "IGE": 1.0, "ITCH": 0.9}},
    "barrier": {"label": "피부장벽 손상 (계면활성제·테이프 스트리핑)",
                "unit": {"cell": "SDS 처리 (HaCaT)", "mouse": "SDS·아세톤 반복 처리", "human": "세정제·건조 환경 노출 지수"},
                "w": {"ECM": 0.2, "COL": 0.2, "SEN": 0.2, "INF": 0.6, "OX": 0.3, "PIG": 0.1, "BAR": 1.0, "MOR": 0.4,
                      "TH2": 0.5, "IGE": 0.3, "ITCH": 0.6}},
}
AXES = {"ECM": "ECM 분해", "COL": "콜라겐 합성", "SEN": "세포노화", "INF": "염증", "OX": "산화스트레스",
        "PIG": "색소침착", "BAR": "피부장벽", "MOR": "피부 형태(주름·두께)", "TH2": "Th2 면역반응",
        "IGE": "IgE 감작", "ITCH": "가려움"}
# 지표: (이름, 기전축, 스트레스 시 방향, 측정 가능한 개체, 최대 변화폭 %)
BIOMARKERS = [
    ("MMP-1", "ECM", +1, ("cell", "mouse"), 180),
    ("Type I collagen (COL1A1·PIP)", "COL", -1, ("cell", "mouse", "human"), 55),
    ("p16INK4a", "SEN", +1, ("cell", "mouse"), 150),
    ("SA-β-gal 양성률", "SEN", +1, ("cell",), 200),
    ("IL-6", "INF", +1, ("cell", "mouse"), 220),
    ("ROS", "OX", +1, ("cell", "mouse"), 160),
    ("멜라닌", "PIG", +1, ("cell", "mouse", "human"), 90),
    ("TEWL", "BAR", +1, ("mouse", "human"), 110),
    ("Filaggrin", "BAR", -1, ("cell", "mouse", "human"), 60),
    ("주름 깊이", "MOR", +1, ("mouse", "human"), 70),
    ("피부 탄력 (R2)", "MOR", -1, ("human",), 40),
    ("표피 두께", "MOR", +1, ("mouse",), 120),
    ("TSLP", "TH2", +1, ("cell", "mouse"), 250),
    ("IL-4·IL-13", "TH2", +1, ("mouse", "human"), 220),
    ("혈청 IgE", "IGE", +1, ("mouse", "human"), 300),
    ("긁기 횟수·가려움 점수", "ITCH", +1, ("mouse", "human"), 200),
]
BM = pd.DataFrame(BIOMARKERS, columns=["biomarker", "axis", "direction", "organisms", "max_change"])
# 시연용 가상 후보물질 (기전 축별 최대 보호효과 0~1). 실제 물질이 아니며 실데이터 입력 시 실제 물질로 대체된다.
COMPOUNDS = {
    "후보물질 A (항산화형)": {"ec50": 12, "hill": 1.2, "E": {"ECM": .5, "COL": .3, "SEN": .55, "INF": .4, "OX": .9, "PIG": .45, "BAR": .2, "MOR": .35, "TH2": .2, "IGE": .1, "ITCH": .2}},
    "후보물질 B (항염·장벽형)": {"ec50": 20, "hill": 1.0, "E": {"ECM": .45, "COL": .2, "SEN": .25, "INF": .9, "OX": .3, "PIG": .2, "BAR": .8, "MOR": .4, "TH2": .55, "IGE": .3, "ITCH": .6}},
    "후보물질 C (콜라겐·항노화형)": {"ec50": 8, "hill": 1.5, "E": {"ECM": .7, "COL": .85, "SEN": .6, "INF": .3, "OX": .25, "PIG": .1, "BAR": .3, "MOR": .7, "TH2": .05, "IGE": .0, "ITCH": .1}},
    "후보물질 D (미백형)": {"ec50": 15, "hill": 1.3, "E": {"ECM": .15, "COL": .1, "SEN": .15, "INF": .25, "OX": .45, "PIG": .9, "BAR": .1, "MOR": .1, "TH2": .05, "IGE": .0, "ITCH": .05}},
    "후보물질 E (면역조절·항알레르기형)": {"ec50": 10, "hill": 1.2, "E": {"ECM": .2, "COL": .1, "SEN": .1, "INF": .7, "OX": .2, "PIG": .05, "BAR": .45, "MOR": .45, "TH2": .9, "IGE": .8, "ITCH": .85}},
}
TIMING = {"pre": ("전처리 (예방)", 1.0), "co": ("동시 처리", 0.8), "post": ("후처리 (개선)", 0.5)}
AXIS_WEIGHT = DISEASES["aging"]["axes"]


def biomarkers_for(org):
    return [b for b, o in zip(BM.biomarker, BM.organisms) if org in o]


# ---------------------------------------------------------------- 기전 시뮬레이터 (정답 함수)
def true_effect(org, stress, level, days, compound, conc, timing):
    """지표별 % 변화(무처리 대조 대비)의 기대값 dict."""
    o, s = ORGANISMS[org], STRESSES[stress]
    d0, d1 = o["dur"]
    dn = np.clip((days - d0) / max(d1 - d0, 1), 0, 1.5)
    dur_f = 1 - np.exp(-2.2 * (0.25 + dn))
    lev_f = level / (level + 0.35) * 1.35  # 포화형 강도 반응 (level 0~1)
    out = {}
    for b, ax, dirn, orgs, mx in BIOMARKERS:
        if org not in orgs:
            continue
        damage = mx * o["sens"] * s["w"][ax] * lev_f * dur_f
        prot = 0.0
        if compound in COMPOUNDS and conc > 0:
            c = COMPOUNDS[compound]
            ec50 = c["ec50"] * o["ec50_scale"]
            hill = conc ** c["hill"] / (ec50 ** c["hill"] + conc ** c["hill"])
            prot = c["E"][ax] * hill * TIMING[timing][1] * o["deliver"]
        out[b] = dirn * damage * (1 - prot)
    return out


def simulate(n=1400, seed=7, noise_cv=0.12):
    """실험 조건을 무작위 생성해 측정값(노이즈 포함)과 정답값을 함께 반환 (long format)."""
    rng = np.random.default_rng(seed)
    rows = []
    comps = list(COMPOUNDS) + ["없음"]
    for i in range(n):
        org = rng.choice(list(ORGANISMS))
        stress = rng.choice(list(STRESSES), p=[0.18, 0.2, 0.17, 0.25, 0.2])
        level = float(rng.choice([0.25, 0.5, 0.75, 1.0]))
        d0, d1 = ORGANISMS[org]["dur"]
        days = float(rng.choice(np.linspace(d0, d1, 4)).round())
        comp = rng.choice(comps, p=[0.17] * 5 + [0.15])
        cmax = ORGANISMS[org]["conc_max"]
        conc = 0.0 if comp == "없음" else float(rng.choice(np.array([0.03, 0.1, 0.3, 1.0]) * cmax))
        timing = rng.choice(list(TIMING)) if comp != "없음" else "pre"
        truth = true_effect(org, stress, level, days, comp, conc, timing)
        for b, v in truth.items():
            obs = v + rng.normal(0, noise_cv * abs(v) + 2.0)
            rows.append(dict(source="시뮬레이션", condition_id=i, organism=org, stress=stress, level=level,
                             days=days, compound=comp, conc=conc, timing=timing, biomarker=b,
                             value=round(obs, 2), truth=round(v, 2)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- 특징 인코딩
def encode(df: pd.DataFrame, compounds: list[str]) -> np.ndarray:
    feats = []
    for org in ORGANISMS:
        feats.append((df.organism == org).astype(float))
    for st_ in STRESSES:
        feats.append((df.stress == st_).astype(float) * df.level)
    feats.append(df.level.astype(float))
    dn = []
    for org, d in zip(df.organism, df.days):
        d0, d1 = ORGANISMS[org]["dur"]
        dn.append((d - d0) / max(d1 - d0, 1))
    feats.append(pd.Series(dn, index=df.index))
    cn = np.array([c / ORGANISMS[o]["conc_max"] for o, c in zip(df.organism, df.conc)])
    lc = np.log10(cn * 100 + 1) / 2  # 0~1
    for comp in compounds:
        feats.append((df.compound == comp).astype(float) * lc)
    for t, (_, f) in TIMING.items():
        feats.append((df.timing == t).astype(float) * (df.conc > 0) * lc)
    return np.column_stack([np.asarray(f, float) for f in feats])


def feature_names(compounds):
    return ([f"개체:{ORGANISMS[o]['label']}" for o in ORGANISMS] + [f"스트레스:{STRESSES[s]['label']}" for s in STRESSES]
            + ["스트레스 강도", "노출 기간"] + [f"물질:{c}" for c in compounds]
            + [f"처리시점:{TIMING[t][0]}" for t in TIMING])


FEATURE_GROUPS = {"개체": "개체:", "스트레스 종류·강도": "스트레스", "노출 기간": "노출 기간",
                  "후보물질·농도": "물질:", "처리 시점": "처리시점:"}


# ---------------------------------------------------------------- 가우시안 프로세스 회귀
class GP:
    def __init__(self, ls=1.0, noise=0.1, sf=1.0):
        self.ls, self.noise, self.sf = ls, noise, sf

    def _k(self, A, B):
        d = (A[:, None, :] - B[None, :, :]) / self.ls
        return self.sf ** 2 * np.exp(-0.5 * np.sum(d * d, axis=2))

    def fit(self, X, y):
        self.X = X
        self.ym, self.ys = y.mean(), y.std() + 1e-9
        yt = (y - self.ym) / self.ys
        K = self._k(X, X) + (self.noise ** 2 + 1e-8) * np.eye(len(X))
        self.L = np.linalg.cholesky(K)
        self.alpha = np.linalg.solve(self.L.T, np.linalg.solve(self.L, yt))
        self.lml = (-0.5 * yt @ self.alpha - np.log(np.diag(self.L)).sum() - 0.5 * len(X) * np.log(2 * np.pi))
        return self

    def predict(self, Xs, return_std=True):
        Ks = self._k(Xs, self.X)
        mu = Ks @ self.alpha
        if not return_std:
            return mu * self.ys + self.ym
        v = np.linalg.solve(self.L, Ks.T)
        var = np.clip(self.sf ** 2 - np.sum(v * v, axis=0), 1e-12, None) + self.noise ** 2
        return mu * self.ys + self.ym, np.sqrt(var) * self.ys


def fit_gp(X, y, ls_grid=(0.8, 1.4, 2.2), noise_grid=(0.1, 0.25)):
    best = None
    for ls in ls_grid:
        for nz in noise_grid:
            try:
                g = GP(ls, nz).fit(X, y)
            except np.linalg.LinAlgError:
                continue
            if best is None or g.lml > best.lml:
                best = g
    return best


class SkinAgeAI:
    """지표별 GP 모델 묶음."""

    def __init__(self, data: pd.DataFrame, max_n=700, seed=0):
        self.compounds = sorted(c for c in data.compound.unique() if c != "없음")
        self.models, self.train = {}, {}
        rng = np.random.default_rng(seed)
        for b, d in data.groupby("biomarker"):
            if len(d) < 15:
                continue
            if len(d) > max_n:  # 실데이터는 모두 유지, 시뮬레이션은 표본추출
                real = d[d.source != "시뮬레이션"]
                sim = d[d.source == "시뮬레이션"]
                sim = sim.iloc[rng.choice(len(sim), max(max_n - len(real), 50), replace=False)]
                d = pd.concat([real, sim])
            X = encode(d, self.compounds)
            self.models[b] = fit_gp(X, d.value.values.astype(float))
            self.train[b] = d

    def predict(self, cond: pd.DataFrame, biomarker):
        m = self.models[biomarker]
        return m.predict(encode(cond, self.compounds))


def cross_validate(data: pd.DataFrame, k=5, max_n=600, seed=1):
    """지표별 K-fold CV: R², RMSE, 95% 예측구간 포함률, (시뮬레이션이면) 정답 대비 R²."""
    rng = np.random.default_rng(seed)
    comps = sorted(c for c in data.compound.unique() if c != "없음")
    rows, pairs = [], []
    for b, d in data.groupby("biomarker"):
        if len(d) < 30:
            continue
        if len(d) > max_n:
            d = d.iloc[rng.choice(len(d), max_n, replace=False)]
        d = d.reset_index(drop=True)
        folds = rng.permutation(len(d)) % k
        pred, sd = np.zeros(len(d)), np.zeros(len(d))
        for f in range(k):
            tr, te = folds != f, folds == f
            g = fit_gp(encode(d[tr], comps), d.value.values[tr])
            pred[te], sd[te] = g.predict(encode(d[te], comps))
        y = d.value.values
        r2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
        cover = np.mean(np.abs(y - pred) <= 1.96 * sd)
        r2t = np.nan
        if "truth" in d and d.truth.notna().all():
            t = d.truth.values
            r2t = 1 - np.sum((t - pred) ** 2) / np.sum((t - t.mean()) ** 2)
        rows.append(dict(biomarker=b, n=len(d), R2=r2, RMSE=float(np.sqrt(np.mean((y - pred) ** 2))),
                         PI95_coverage=cover, R2_vs_truth=r2t))
        pairs.append(pd.DataFrame({"biomarker": b, "observed": y, "predicted": pred}))
    return pd.DataFrame(rows), pd.concat(pairs, ignore_index=True)


def permutation_importance(ai: SkinAgeAI, biomarker, n_rep=3, seed=2):
    rng = np.random.default_rng(seed)
    d = ai.train[biomarker]
    X = encode(d, ai.compounds)
    y = d.value.values
    base = np.mean((ai.models[biomarker].predict(X, return_std=False) - y) ** 2)
    names = feature_names(ai.compounds)
    out = {}
    for grp, prefix in FEATURE_GROUPS.items():
        idx = [i for i, n in enumerate(names) if n.startswith(prefix)]
        inc = []
        for _ in range(n_rep):
            Xp = X.copy()
            perm = rng.permutation(len(X))
            Xp[:, idx] = X[perm][:, idx]
            inc.append(np.mean((ai.models[biomarker].predict(Xp, return_std=False) - y) ** 2) - base)
        out[grp] = max(np.mean(inc), 0)
    s = sum(out.values()) or 1
    return pd.Series({k: v / s * 100 for k, v in out.items()}).sort_values(ascending=False)


# ---------------------------------------------------------------- 예방효과
def make_cond(org, stress, level, days, compound, conc, timing):
    return pd.DataFrame([dict(organism=org, stress=stress, level=level, days=days,
                              compound=compound, conc=conc, timing=timing)])


def prevention_profile(ai: SkinAgeAI, org, stress, level, days, compound, conc, timing, disease=None):
    """지표별: 스트레스 단독 변화, 병용 변화(±95%), 예방률. disease 를 주면 그 질환의 기전 축 지표만."""
    rows = []
    bm_axis = BM.set_index("biomarker")["axis"]
    for b in biomarkers_for(org):
        if disease and bm_axis[b] not in DISEASES[disease]["axes"]:
            continue
        if b not in ai.models:
            continue
        s_mu, s_sd = ai.predict(make_cond(org, stress, level, days, "없음", 0.0, "pre"), b)
        c_mu, c_sd = ai.predict(make_cond(org, stress, level, days, compound, conc, timing), b)
        s_mu, c_mu = float(s_mu[0]), float(c_mu[0])
        pi = (s_mu - c_mu) / s_mu * 100 if abs(s_mu) > 10 else np.nan  # 스트레스 영향이 작은 지표는 예방률 미산출
        ax = BM.set_index("biomarker").loc[b, "axis"]
        rows.append(dict(biomarker=b, axis=AXES[ax], axis_key=ax, stress_only=s_mu, with_compound=c_mu,
                         ci95=1.96 * float(c_sd[0]), prevention=pi))
    df = pd.DataFrame(rows)
    return df


def composite_score(profile: pd.DataFrame, disease="aging"):
    """질환 관련 기전 축 가중 예방점수 (%)."""
    ax_w = DISEASES[disease]["axes"]
    p = profile.dropna(subset=["prevention"])
    p = p[p.axis_key.isin(ax_w)]
    if p.empty:
        return np.nan
    w = p.axis_key.map(ax_w)
    return float(np.clip((p.prevention.clip(-50, 120) * w).sum() / w.sum(), -50, 120))


def disease_index(values: dict, org, disease):
    """지표 변화(%) → 질환 지수 0~100 (관련 축 손상 정도의 가중 평균, 최대 변화폭으로 정규화)."""
    ax_w = DISEASES[disease]["axes"]
    bm = BM.set_index("biomarker")
    num = den = 0.0
    for b, v in values.items():
        ax = bm.loc[b, "axis"]
        if ax not in ax_w:
            continue
        sev = np.clip(v * bm.loc[b, "direction"] / (bm.loc[b, "max_change"] * ORGANISMS[org]["sens"]), 0, 1.2)
        num += ax_w[ax] * sev
        den += ax_w[ax]
    return float(np.clip(num / den * 100, 0, 100)) if den else np.nan


def disease_prediction(profile: pd.DataFrame, org, disease):
    s = disease_index(dict(zip(profile.biomarker, profile.stress_only)), org, disease)
    c = disease_index(dict(zip(profile.biomarker, profile.with_compound)), org, disease)
    return {"stress_only": s, "with_compound": c, "prevention": (s - c) / s * 100 if s and s > 1 else np.nan}


def dose_response(ai: SkinAgeAI, org, stress, level, days, compound, timing, n=14, disease="aging"):
    cmax = ORGANISMS[org]["conc_max"]
    concs = np.unique(np.concatenate([[0], np.geomspace(cmax * 0.005, cmax, n)]))
    scores = []
    for c in concs:
        prof = prevention_profile(ai, org, stress, level, days, compound, c, timing, disease) if c > 0 else None
        scores.append(0.0 if c == 0 else composite_score(prof, disease))
    return pd.DataFrame({"conc": concs, "score": scores})


def fit_4pl(conc, resp):
    """resp = top * c^h / (ec50^h + c^h) (bottom=0). 격자탐색 + 선형 top 추정."""
    c = np.asarray(conc, float)
    r = np.asarray(resp, float)
    m = c > 0
    c, r = c[m], r[m]
    best = None
    for le in np.linspace(np.log10(c.min() / 3), np.log10(c.max() * 3), 120):
        for h in np.linspace(0.5, 3, 26):
            g = c ** h / ((10 ** le) ** h + c ** h)
            top = float(g @ r / max(g @ g, 1e-12))
            sse = float(np.sum((r - top * g) ** 2))
            if best is None or sse < best[0]:
                best = (sse, 10 ** le, h, top)
    _, ec50, h, top = best
    rec = ec50 * (0.8 / 0.2) ** (1 / h)  # 최대효과의 80% 도달 농도
    return dict(ec50=ec50, hill=h, top=top, rec80=rec)


# ---------------------------------------------------------------- 실데이터 양식
TEMPLATE_COLUMNS = ["organism", "model", "stress", "level", "days", "compound", "conc", "timing", "biomarker", "value", "pmid", "note"]


def template_df():
    return pd.DataFrame([
        dict(organism="cell", model="HDF", stress="uv", level=0.5, days=2, compound="물질X", conc=10, timing="pre",
             biomarker="MMP-1", value=85.0, note="value = 무처리 대조 대비 % 변화 (스트레스로 증가하면 +)"),
        dict(organism="cell", model="HDF", stress="uv", level=0.5, days=2, compound="없음", conc=0, timing="pre",
             biomarker="MMP-1", value=160.0, note="스트레스 단독군도 함께 입력"),
        dict(organism="mouse", model="HR-1", stress="uv", level=0.75, days=56, compound="물질X", conc=0.5, timing="co",
             biomarker="주름 깊이", value=30.0, pmid="", note="level: 약 0.25 / 중 0.5 / 강 0.75 / 최강 1.0"),
        dict(organism="mouse", model="BALB/c", stress="allergen", level=0.75, days=28, compound="물질Y", conc=1.0,
             timing="co", biomarker="혈청 IgE", value=120.0, pmid="12345678", note="문헌값은 PMID 입력 (예시 번호)"),
    ], columns=TEMPLATE_COLUMNS)


def standardize_real(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    miss = [c for c in ["organism", "stress", "level", "days", "compound", "conc", "timing", "biomarker", "value"] if c not in df]
    if miss:
        raise ValueError("필수 컬럼이 없습니다: " + ", ".join(miss))
    df["organism"] = df.organism.str.lower().str.strip()
    df["stress"] = df.stress.str.lower().str.strip()
    df["timing"] = df.timing.str.lower().str.strip()
    bad = df[~df.organism.isin(ORGANISMS) | ~df.stress.isin(STRESSES) | ~df.timing.isin(TIMING)]
    if len(bad):
        raise ValueError("organism(cell/mouse/human), stress(uv/ox/gly/allergen/barrier), timing(pre/co/post) "
                         f"값을 확인하세요: {len(bad)}행")
    unknown = sorted(set(df.biomarker) - set(BM.biomarker))
    if unknown:
        raise ValueError("지표 이름은 다음 중 하나여야 합니다: " + ", ".join(BM.biomarker) + f" (확인 필요: {', '.join(unknown)})")
    for c in ["level", "days", "conc", "value"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["level", "days", "conc", "value"])
    pm = df["pmid"].astype(str).str.strip() if "pmid" in df else pd.Series("", index=df.index)
    df["source"] = np.where(pm.str.fullmatch(r"\d{5,9}"), "문헌", "실데이터")
    df["truth"] = np.nan
    return df


def reversal_score(deg: pd.DataFrame, panel_rows):
    """후보물질 전사체 변화가 피부노화 방향을 되돌리는 정도 (-1~1, 양수 = 노화 역전)."""
    d = deg.set_index("gene")
    s, n, used = 0.0, 0, []
    for r in panel_rows:
        if r["expected"] == "var" or r["gene"] not in d.index:
            continue
        lfc = float(d.loc[r["gene"], "log2fc"])
        sign = 1 if r["expected"] == "up" else -1
        s += -sign * np.tanh(lfc)
        n += 1
        used.append((r["gene"], r["category"], r["expected"], lfc, -sign * np.tanh(lfc)))
    return (s / n if n else np.nan), pd.DataFrame(used, columns=["gene", "category", "aging_dir", "log2fc", "reversal"])
