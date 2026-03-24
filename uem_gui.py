"""UEM 模型可视化测试界面"""
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pycalphad import Database, calculate, equilibrium, variables as v
from pycalphad.models.model_uem import ModelUEM, ModelMuggianu, ModelKohler
import time

dbf = Database('examples/alcrni.tdb')
comps = ['AL', 'NI', 'CR', 'VA']
phases = ['LIQUID']
T = 1800

print("正在计算...")

# 三个模型计算
t0 = time.perf_counter()
res_uem = calculate(dbf, comps, phases, model=ModelUEM, T=T, P=101325, output='GM')
t1 = time.perf_counter()
res_mug = calculate(dbf, comps, phases, model=ModelMuggianu, T=T, P=101325, output='GM')
t2 = time.perf_counter()
res_rkm = calculate(dbf, comps, phases, T=T, P=101325, output='GM')
t3 = time.perf_counter()

print(f"UEM1:     {t1-t0:.2f}s")
print(f"Muggianu: {t2-t1:.2f}s")
print(f"R-K-M:    {t3-t2:.2f}s")

# 提取数据: 沿 AL-NI 截面 (CR 固定)
# 取 X(CR) ≈ 0.1 的截面
fig = plt.figure(figsize=(14, 10))
fig.suptitle(f'Al-Cr-Ni LIQUID @ T={T}K — UEM vs Traditional Models', fontsize=14)
gs = GridSpec(2, 2, figure=fig)

# --- Plot 1: GM vs X(AL) at X(CR)≈0.1 ---
ax1 = fig.add_subplot(gs[0, 0])

for res, label, color in [(res_uem, 'UEM1', 'red'),
                           (res_mug, 'Muggianu', 'blue'),
                           (res_rkm, 'R-K-M', 'green')]:
    gm = res.GM.values.flatten()
    x_al = res.X.sel(component='AL').values.flatten()
    x_cr = res.X.sel(component='CR').values.flatten()

    mask = (x_cr > 0.08) & (x_cr < 0.12)
    if mask.sum() > 0:
        order = np.argsort(x_al[mask])
        ax1.plot(x_al[mask][order], gm[mask][order], '.', label=label, color=color, markersize=2, alpha=0.5)

ax1.set_xlabel('X(AL)')
ax1.set_ylabel('GM (J/mol)')
ax1.set_title('X(CR) ≈ 0.1 截面')
ax1.legend()
ax1.grid(True, alpha=0.3)

# --- Plot 2: GM vs X(AL) at X(CR)≈0.3 ---
ax2 = fig.add_subplot(gs[0, 1])

for res, label, color in [(res_uem, 'UEM1', 'red'),
                           (res_mug, 'Muggianu', 'blue'),
                           (res_rkm, 'R-K-M', 'green')]:
    gm = res.GM.values.flatten()
    x_al = res.X.sel(component='AL').values.flatten()
    x_cr = res.X.sel(component='CR').values.flatten()

    mask = (x_cr > 0.28) & (x_cr < 0.32)
    if mask.sum() > 0:
        order = np.argsort(x_al[mask])
        ax2.plot(x_al[mask][order], gm[mask][order], '.', label=label, color=color, markersize=2, alpha=0.5)

ax2.set_xlabel('X(AL)')
ax2.set_ylabel('GM (J/mol)')
ax2.set_title('X(CR) ≈ 0.3 截面')
ax2.legend()
ax2.grid(True, alpha=0.3)

# --- Plot 3: 差值 (UEM - RKM) 等值线 ---
ax3 = fig.add_subplot(gs[1, 0])

gm_uem = res_uem.GM.values.flatten()
gm_rkm = res_rkm.GM.values.flatten()
x_al = res_uem.X.sel(component='AL').values.flatten()
x_cr = res_uem.X.sel(component='CR').values.flatten()

diff = gm_uem - gm_rkm
valid = np.isfinite(diff) & (x_al + x_cr <= 1.01)

sc = ax3.scatter(x_al[valid], x_cr[valid], c=diff[valid], cmap='RdBu_r',
                  s=1, alpha=0.6, vmin=-200, vmax=200)
plt.colorbar(sc, ax=ax3, label='ΔGM (J/mol)')
ax3.set_xlabel('X(AL)')
ax3.set_ylabel('X(CR)')
ax3.set_title('GM(UEM1) - GM(R-K-M)')
ax3.set_xlim(0, 1)
ax3.set_ylim(0, 1)
ax3.plot([0, 1], [1, 0], 'k-', lw=0.5)
ax3.grid(True, alpha=0.3)

# --- Plot 4: 统计信息 ---
ax4 = fig.add_subplot(gs[1, 1])
ax4.axis('off')

info = f"""
  模型对比统计 (T = {T} K, LIQUID)
  ─────────────────────────────
  UEM1:
    GM range: [{gm_uem[np.isfinite(gm_uem)].min():.1f}, {gm_uem[np.isfinite(gm_uem)].max():.1f}]
    计算耗时: {t1-t0:.3f}s

  Muggianu:
    GM range: [{res_mug.GM.values.min():.1f}, {res_mug.GM.values.max():.1f}]
    计算耗时: {t2-t1:.3f}s

  R-K-M (原始):
    GM range: [{gm_rkm[np.isfinite(gm_rkm)].min():.1f}, {gm_rkm[np.isfinite(gm_rkm)].max():.1f}]
    计算耗时: {t3-t2:.3f}s

  ─────────────────────────────
  UEM1 vs R-K-M 差值:
    max|ΔGM| = {np.abs(diff[valid]).max():.1f} J/mol
    mean|ΔGM| = {np.abs(diff[valid]).mean():.1f} J/mol
"""
ax4.text(0.05, 0.95, info, transform=ax4.transAxes, fontsize=10,
         verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.show()
