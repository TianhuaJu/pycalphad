"""
R-K-UEM 模型示例脚本

该脚本演示如何使用R-K-UEM模块计算多元合金系统的热力学性质。

示例包括:
1. 二元系统的混合焓和过剩Gibbs能计算
2. 三元系统使用不同外推模型的计算
3. 活度系数和活度计算
4. 不同外推模型的比较
"""

import sys
import os

# 添加pycalphad到路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from pycalphad.models.rk_uem import (
    RKBinaryPolynomial,
    ThermodynamicCalculator,
    UEM1, UEM2_N, Kohler, Muggianu
)


def example_binary_system(database_path):
    """
    示例1: 二元系统计算
    """
    print("=" * 70)
    print("示例 1: 二元系统 - R-K多项式计算")
    print("=" * 70)

    # 创建二元系统对象（需要实际的数据库文件）
    # 这里使用占位符，实际使用时需要提供真实的数据库路径
    try:
        binary = RKBinaryPolynomial(("Fe", "Cr"), database_path)

        T = 1800.0  # K
        x_Fe = 0.5

        # 计算混合焓
        H_mix = binary.mixing_enthalpy("Fe", x_Fe, T)
        print(f"Fe-Cr 二元系 @ T={T}K, x_Fe={x_Fe}")
        print(f"  混合焓: {H_mix:.2f} J/mol" if H_mix is not None else "  混合焓: 数据缺失")

        # 计算过剩Gibbs能
        G_ex = binary.excess_gibbs_energy("Fe", x_Fe, T)
        print(f"  过剩Gibbs能: {G_ex:.2f} J/mol" if G_ex is not None else "  过剩Gibbs能: 数据缺失")

        # 计算无限稀释性质
        W_Fe_inf = binary.infinite_dilution_property("Fe", T)
        print(f"  Fe在Cr中的无限稀释性质: {W_Fe_inf:.2f} J/mol")

    except Exception as e:
        print(f"二元系统计算失败: {e}")
        print("提示: 需要提供包含Fe-Cr参数的R-K数据库文件")

    print()


def example_ternary_system(database_path):
    """
    示例2: 三元系统计算 - 比较不同外推模型
    """
    print("=" * 70)
    print("示例 2: 三元系统 - UEM外推模型")
    print("=" * 70)

    try:
        # 创建热力学计算器
        calc = ThermodynamicCalculator(database_path)

        # 定义三元系组成
        composition = {
            "Fe": 0.5,
            "Cr": 0.3,
            "Ni": 0.2
        }

        T = 1800.0  # K

        print(f"系统: Fe-Cr-Ni")
        print(f"组成: {composition}")
        print(f"温度: {T} K")
        print()

        # 使用不同外推模型计算
        models = [
            ("UEM1", UEM1),
            ("Kohler", Kohler),
            ("Muggianu", Muggianu),
        ]

        for model_name, model_func in models:
            try:
                # 计算混合焓
                H_mix = calc.get_mixing_enthalpy(composition, T, model_func)

                # 计算过剩Gibbs能
                G_ex = calc.get_excess_gibbs(composition, T, model_func)

                print(f"{model_name:12s} 模型:")
                print(f"  混合焓:         {H_mix:>12.2f} J/mol")
                print(f"  过剩Gibbs能:    {G_ex:>12.2f} J/mol")

                if calc.warnings:
                    print(f"  警告: {calc.warnings}")
                    calc.warnings.clear()

            except Exception as e:
                print(f"{model_name:12s} 模型: 计算失败 - {e}")

            print()

    except Exception as e:
        print(f"三元系统计算失败: {e}")
        print("提示: 需要提供包含Fe-Cr, Fe-Ni, Cr-Ni参数的R-K数据库文件")

    print()


def example_activity_coefficient(database_path):
    """
    示例3: 活度系数计算
    """
    print("=" * 70)
    print("示例 3: 活度系数计算")
    print("=" * 70)

    try:
        calc = ThermodynamicCalculator(database_path)

        composition = {
            "Fe": 0.6,
            "Cr": 0.25,
            "Ni": 0.15
        }

        T = 1800.0
        solvent = "Fe"  # 参考组分

        print(f"系统: Fe-Cr-Ni")
        print(f"组成: {composition}")
        print(f"温度: {T} K")
        print(f"溶剂 (参考组分): {solvent}")
        print()

        # 计算所有组分的活度系数
        ln_gammas = calc.calculate_all_activity_coefficients(
            composition, solvent, T, UEM1
        )

        print("活度系数 (ln γ):")
        for comp, ln_gamma in ln_gammas.items():
            import math
            gamma = math.exp(ln_gamma)
            print(f"  {comp:3s}: ln(γ) = {ln_gamma:>8.4f}, γ = {gamma:>8.4f}")

        print()

        # 计算活度
        print("活度:")
        for comp in composition.keys():
            activity = calc.calculate_activity(composition, comp, solvent, T, UEM1)
            print(f"  {comp:3s}: a = {activity:>8.4f}")

    except Exception as e:
        print(f"活度系数计算失败: {e}")
        print("提示: 需要提供包含所有二元对参数的R-K数据库文件")

    print()


def example_composition_scan(database_path):
    """
    示例4: 组成扫描 - 绘制混合焓随组成变化
    """
    print("=" * 70)
    print("示例 4: 组成扫描")
    print("=" * 70)

    try:
        calc = ThermodynamicCalculator(database_path)
        T = 1800.0

        print(f"扫描 Fe-Cr 二元系的混合焓")
        print(f"温度: {T} K")
        print()
        print("x_Fe      H_mix (J/mol)")
        print("-" * 30)

        for x_fe in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
            composition = {
                "Fe": x_fe,
                "Cr": 1.0 - x_fe
            }

            try:
                H_mix = calc.get_mixing_enthalpy(composition, T, UEM1)
                print(f"{x_fe:.1f}     {H_mix:>12.2f}")
            except Exception as e:
                print(f"{x_fe:.1f}     计算失败")

    except Exception as e:
        print(f"组成扫描失败: {e}")

    print()


def print_usage():
    """打印使用说明"""
    print("=" * 70)
    print("R-K-UEM 模型使用示例")
    print("=" * 70)
    print()
    print("该脚本演示了R-K-UEM模块的基本用法。")
    print()
    print("注意:")
    print("  1. 需要提供包含R-K参数的SQLite数据库文件")
    print("  2. 数据库应包含表: R_K_gE 和 R_K_Hmix")
    print("  3. 表结构应包含: Symbol, A0, A1, A2, A3, A4, A5")
    print()
    print("使用方法:")
    print("  python rk_uem_example.py [database_path]")
    print()
    print("如果不提供数据库路径，将运行演示模式（仅展示API用法）")
    print("=" * 70)
    print()


def main():
    """主函数"""
    print_usage()

    # 检查命令行参数
    if len(sys.argv) > 1:
        database_path = sys.argv[1]
        if not os.path.exists(database_path):
            print(f"错误: 数据库文件不存在: {database_path}")
            sys.exit(1)
    else:
        database_path = None
        print("警告: 未提供数据库路径，运行演示模式")
        print("      实际计算需要R-K参数数据库")
        print()

    # 运行示例
    if database_path:
        example_binary_system(database_path)
        example_ternary_system(database_path)
        example_activity_coefficient(database_path)
        example_composition_scan(database_path)
    else:
        print("演示模式: API用法示例")
        print()
        print("# 1. 创建二元系统对象")
        print("binary = RKBinaryPolynomial(('Fe', 'Cr'), database_path)")
        print("H_mix = binary.mixing_enthalpy('Fe', 0.5, 1800.0)")
        print()
        print("# 2. 创建热力学计算器")
        print("calc = ThermodynamicCalculator(database_path)")
        print()
        print("# 3. 计算多元系性质")
        print("composition = {'Fe': 0.5, 'Cr': 0.3, 'Ni': 0.2}")
        print("H_mix = calc.get_mixing_enthalpy(composition, 1800.0, UEM1)")
        print()
        print("# 4. 计算活度系数")
        print("ln_gamma = calc.calculate_activity_coefficient(")
        print("    composition, 'Fe', 'Ni', 1800.0, UEM1")
        print(")")

    print()
    print("示例运行完成！")


if __name__ == "__main__":
    main()
