@echo off
REM 启动合金液相线/固相线计算GUI (Windows)

echo 启动合金相图计算GUI...
echo ================================

REM 检查Python是否安装
python --version >nul 2>&1
if errorlevel 1 (
    echo 错误: 未找到 Python
    echo 请先安装 Python 3.7+
    pause
    exit /b 1
)

REM 启动GUI
echo 启动GUI程序...
python phase_diagram_gui.py

echo ================================
echo 程序已退出
pause
