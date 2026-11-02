"""AI(머신러닝) 모듈: 전사체 기반 피부노화 점수.

- PCA (비지도 학습): 샘플 간 전사체 구조 확인
- L2-정규화 로지스틱 회귀 (지도 학습): 젊은/노화 분류기 → '전사체 피부노화 점수'(0~1)
- Leave-one-out 교차검증: 특징 선택(분산 상위 + 피부노화 패널)을 매 fold 안에서만 수행해 정보 누수 방지
- 특징 중요도: 표준화 계수 크기 순위
numpy 만 사용 (외부 ML 라이브러리 의존 없음).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .skin_aging_genes import all_panel_genes


def log_expr(norm: pd.DataFrame) -> pd.DataFrame:
    x = np.log2(norm + 1)
    x.index = x.index.astype(str).str.upper()
    return x[~x.index.duplicated()]


def pca(norm: pd.DataFrame, n_top=1000, n_comp=2):
    x = log_expr(norm)
    x = x.loc[x.var(axis=1).nlargest(min(n_top, len(x))).index]
    xc = (x.T - x.T.mean()).values
    u, s, vt = np.linalg.svd(xc, full_matrices=False)
    pcs = u[:, :n_comp] * s[:n_comp]
    var = s ** 2 / np.sum(s ** 2)
    return pd.DataFrame(pcs, index=x.columns, columns=[f"PC{i+1}" for i in range(n_comp)]), var[:n_comp]


def _fit_logreg(X, y, lam=1.0, iters=100):
    """Newton-Raphson L2 로지스틱 회귀 (절편 비정규화)."""
    n, p = X.shape
    Xb = np.column_stack([np.ones(n), X])
    w = np.zeros(p + 1)
    R = np.eye(p + 1) * lam
    R[0, 0] = 0
    for _ in range(iters):
        z = np.clip(Xb @ w, -30, 30)
        pr = 1 / (1 + np.exp(-z))
        g = Xb.T @ (pr - y) + R @ w
        H = Xb.T @ (Xb * (pr * (1 - pr))[:, None]) + R
        step = np.linalg.solve(H + 1e-9 * np.eye(p + 1), g)
        w -= step
        if np.max(np.abs(step)) < 1e-8:
            break
    return w


def _predict(w, X):
    z = np.clip(np.column_stack([np.ones(len(X)), X]) @ w, -30, 30)
    return 1 / (1 + np.exp(-z))


def _select(x: pd.DataFrame, n_var):
    panel = [g for g in all_panel_genes() if g in x.index]
    top = list(x.var(axis=1).nlargest(n_var).index)
    return list(dict.fromkeys(panel + top))


def _auc(y, s):
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    return float(np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg]))


def aging_classifier(norm: pd.DataFrame, groups: pd.Series, aged_label: str,
                     n_var=200, lam=10.0):
    """LOOCV 로 노화 점수와 성능을 계산하고, 전체 데이터로 최종 모델을 학습한다."""
    x = log_expr(norm)
    y = (groups.loc[x.columns].astype(str).values == str(aged_label)).astype(float)
    n = len(y)
    loo = np.zeros(n)
    for i in range(n):
        tr = np.arange(n) != i
        xt = x.iloc[:, tr]
        feats = _select(xt, n_var)
        mu, sd = xt.loc[feats].mean(axis=1), xt.loc[feats].std(axis=1).replace(0, 1)
        Xtr = ((xt.loc[feats].T - mu) / sd).values
        Xte = ((x.loc[feats].iloc[:, [i]].T - mu) / sd).values
        w = _fit_logreg(Xtr, y[tr], lam)
        loo[i] = _predict(w, Xte)[0]
    feats = _select(x, n_var)
    mu, sd = x.loc[feats].mean(axis=1), x.loc[feats].std(axis=1).replace(0, 1)
    X = ((x.loc[feats].T - mu) / sd).values
    w = _fit_logreg(X, y, lam)
    imp = pd.DataFrame({"gene": feats, "weight": w[1:]})
    imp["abs"] = imp["weight"].abs()
    imp["방향"] = np.where(imp["weight"] > 0, "노화 쪽", "젊음 쪽")
    imp["피부노화 패널"] = imp["gene"].isin(all_panel_genes())
    imp = imp.sort_values("abs", ascending=False).drop(columns="abs").reset_index(drop=True)
    pred = (loo >= 0.5).astype(float)
    other = [g for g in groups.astype(str).unique() if g != str(aged_label)]
    other = other[0] if len(other) == 1 else "기타"
    scores = pd.DataFrame({"sample": x.columns, "group": groups.loc[x.columns].values,
                           "노화점수_LOOCV": loo.round(3), "예측": np.where(pred == 1, aged_label, other)})
    metrics = {"n_samples": n, "n_features": len(feats), "loocv_accuracy": float((pred == y).mean()),
               "loocv_auc": _auc(y, loo)}
    model = {"features": feats, "mu": mu, "sd": sd, "w": w}
    return scores, imp, metrics, model


def score_new(model, norm_new: pd.DataFrame):
    """학습된 모델로 새 샘플의 노화점수 계산 (없는 유전자는 평균값으로 대체)."""
    x = log_expr(norm_new).reindex(model["features"])
    x = x.T.fillna(model["mu"]).T
    X = ((x.T - model["mu"]) / model["sd"]).values
    return pd.Series(_predict(model["w"], X), index=x.columns)
