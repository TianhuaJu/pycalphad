@echo off
echo ============================================================
echo 清理Python缓存
echo ============================================================
echo.

echo 正在删除__pycache__目录...
for /d /r pycalphad %%d in (__pycache__) do (
    if exist "%%d" (
        echo 删除: %%d
        rd /s /q "%%d"
    )
)

echo.
echo 正在删除.pyc文件...
del /s /q pycalphad\*.pyc 2>nul

echo.
echo ============================================================
echo 缓存清理完成！
echo ============================================================
echo.
echo 现在请运行: python test_uem_example.py
echo.
pause
