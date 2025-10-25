#!/bin/bash
# 启动合金液相线/固相线计算GUI

echo "启动合金相图计算GUI..."
echo "================================"

# 检查Python是否安装
if ! command -v python3 &> /dev/null; then
    echo "错误: 未找到 Python3"
    echo "请先安装 Python 3.7+"
    exit 1
fi

# 检查依赖
echo "检查依赖..."
python3 -c "import pycalphad, numpy, matplotlib" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "警告: 部分依赖未安装"
    echo "正在安装依赖..."
    pip3 install pycalphad numpy matplotlib
fi

# 启动GUI
echo "启动GUI程序..."
python3 phase_diagram_gui.py

echo "================================"
echo "程序已退出"
