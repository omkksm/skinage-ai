@echo off
chcp 65001 > nul
cd /d %~dp0
echo [SkinAge-Omics] 필요한 패키지를 설치합니다 (처음 한 번만 시간이 걸립니다)...
python -m pip install -r requirements.txt
echo [SkinAge-Omics] 앱을 실행합니다. 브라우저가 자동으로 열립니다.
python -m streamlit run streamlit_app.py
pause
