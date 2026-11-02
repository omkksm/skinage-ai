"""앱 화면 자동 캡처 (증빙용).

사용법:
  1) 터미널 1:  streamlit run app.py
  2) 터미널 2:  python scripts/capture_screenshots.py
     → screenshots/ 폴더에 탭별 전체 화면 PNG 저장
처음 한 번:  pip install playwright  &&  python -m playwright install chromium
"""
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

URL = sys.argv[1] if len(sys.argv) > 1 else "http://localhost:8501"
OUT = Path(__file__).resolve().parents[1] / "screenshots"
OUT.mkdir(exist_ok=True)
TABS = ["① DEG 시각화", "② 피부노화 유전자 탐색", "③ GO·KEGG 경로 분석",
        "④ 기존 분석과 비교", "⑤ AI 피부노화 점수", "ⓘ 정보"]

with sync_playwright() as p:
    b = p.chromium.launch()
    pg = b.new_page(viewport={"width": 1600, "height": 1000}, device_scale_factor=1.5)
    pg.goto(URL)
    pg.wait_for_selector("text=SkinAge-Omics", timeout=60000)
    time.sleep(8)
    for i, t in enumerate(TABS, 1):
        pg.get_by_role("tab", name=t).click()
        if t.startswith("④"):
            try:
                pg.get_by_text("앱 내 다른 방법").click()
            except Exception:
                pass
        time.sleep(6)
        pg.screenshot(path=str(OUT / f"{i:02d}_{t[2:].replace(' ', '_').replace('·', '')}.png"), full_page=True)
        print("saved", t)
    b.close()
print("완료:", OUT)
