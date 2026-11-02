"""피부노화 관련 큐레이션 유전자 패널.

각 유전자에는 기능 범주와 '노화(또는 광노화) 피부에서 문헌상 보고되는 일반적 변화 방향'을
붙였다. 방향은 해석을 돕기 위한 참고값이며, 조직/연령대/광노출 여부에 따라 달라질 수 있다.
  "up"   : 노화 피부에서 증가 경향
  "down" : 노화 피부에서 감소 경향
  "var"  : 연구마다 방향이 다르거나 상황 의존적
"""

SKIN_AGING_PANEL = {
    "ECM 분해 (MMP)": {
        "MMP1": "up", "MMP2": "up", "MMP3": "up", "MMP9": "up",
        "MMP12": "up", "MMP13": "up",
    },
    "콜라겐·ECM 합성": {
        "COL1A1": "down", "COL1A2": "down", "COL3A1": "down", "COL4A1": "down",
        "COL7A1": "down", "COL17A1": "down", "ELN": "var", "FBN1": "down",
        "DCN": "down", "LUM": "down", "HAS2": "down", "HAS3": "down",
    },
    "MMP 억제 (TIMP)": {
        "TIMP1": "var", "TIMP2": "var", "TIMP3": "down",
    },
    "세포노화·세포주기": {
        "CDKN2A": "up", "CDKN1A": "up", "TP53": "up", "LMNB1": "down",
        "GLB1": "up", "SERPINE1": "up",
    },
    "SASP·염증": {
        "IL6": "up", "CXCL8": "up", "IL1A": "up", "IL1B": "up", "TNF": "up",
        "CCL2": "up", "NFKB1": "up", "PTGS2": "up",
    },
    "산화스트레스·항산화": {
        "SOD1": "down", "SOD2": "var", "CAT": "down", "GPX1": "down",
        "NFE2L2": "down", "HMOX1": "var", "NQO1": "var",
    },
    "표피 장벽·분화": {
        "FLG": "down", "LOR": "down", "IVL": "var", "KRT10": "down",
        "KRT14": "down", "AQP3": "down", "CLDN1": "down", "TGM1": "var",
    },
    "장수·대사 신호": {
        "SIRT1": "down", "SIRT3": "down", "SIRT6": "down", "FOXO3": "down",
        "MTOR": "up", "IGF1": "down", "IGF1R": "var", "PPARGC1A": "down",
    },
    "TGF-β·AP-1 신호": {
        "TGFB1": "var", "TGFBR2": "down", "SMAD3": "down", "SMAD7": "up",
        "JUN": "up", "FOS": "up", "CCN2": "down", "CCN1": "up",
    },
    "색소·광노화": {
        "MITF": "var", "TYR": "var", "DCT": "var", "PMEL": "var",
    },
}

# 흔히 쓰이는 이전 유전자 기호 → 현재 HGNC 기호
ALIASES = {"IL8": "CXCL8", "CTGF": "CCN2", "CYR61": "CCN1", "P16": "CDKN2A",
           "P21": "CDKN1A", "NRF2": "NFE2L2", "PGC1A": "PPARGC1A"}


def panel_table():
    """패널을 (gene, category, expected) 행 목록으로 반환."""
    rows = []
    for cat, genes in SKIN_AGING_PANEL.items():
        for g, d in genes.items():
            rows.append({"gene": g, "category": cat, "expected": d})
    return rows


def all_panel_genes():
    return [r["gene"] for r in panel_table()]
