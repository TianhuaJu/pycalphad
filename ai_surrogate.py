"""
AI 代理加速模块 — 基于 RBF 插值的液相面快速计算

使用稀疏锚点 + 径向基函数插值 + 自适应细化，将三元液相面计算
加速约 5-10 倍，同时通过 progress_callback 支持渐进式 GUI 显示。

依赖: numpy, scipy（均为项目已有依赖，无需额外安装）
"""

import hashlib
import os
import pickle
import time
import traceback
from collections import OrderedDict
from enum import Enum
from multiprocessing import Pool
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.interpolate import RBFInterpolator
from scipy.spatial.distance import cdist

from thermocal_core import (
	CalculationMode,
	CalculationResult,
	LiquidusSolidusCalculator,
	LiquidusSurfaceCalculator,
	_surface_calculate_point,
	_surface_worker_init,
	_worker_state,
)


# ============================================================================
# 计算阶段枚举
# ============================================================================

class SurrogateStage(Enum):
	"""代理计算的各阶段"""
	SPARSE   = 'sparse'     # 阶段1：稀疏网格计算
	PREVIEW  = 'preview'    # 阶段2：RBF插值预览
	REFINING = 'refining'   # 阶段3：自适应细化
	REFINED  = 'refined'    # 阶段4：细化后最终结果


# ============================================================================
# 缓存
# ============================================================================

class SurrogateCache:
	"""液相面代理模型缓存

	按系统签名（组分+相+模型+温度范围）缓存已计算的锚点数据，
	避免重复计算同一系统。FIFO 淘汰，上限 max_entries 条。
	"""

	def __init__(self, max_entries: int = 10):
		self._cache: OrderedDict = OrderedDict()
		self._max_entries = max_entries

	@staticmethod
	def make_key(components: List[str],
	             phases: Optional[List[str]],
	             model_spec,
	             T_range: Tuple[float, float]) -> str:
		"""生成缓存键（SHA256 前16位）"""
		key_parts = [
			','.join(sorted(c.upper() for c in components)),
			','.join(sorted(p.upper() for p in (phases or []))),
			repr(model_spec),
			f"{T_range[0]:.1f}-{T_range[1]:.1f}",
		]
		raw = '|'.join(key_parts)
		return hashlib.sha256(raw.encode()).hexdigest()[:16]

	def get(self, key: str) -> Optional[Dict]:
		"""获取缓存数据。命中时移到末尾（LRU）"""
		if key in self._cache:
			self._cache.move_to_end(key)
			return self._cache[key]
		return None

	def put(self, key: str, xb: np.ndarray, xc: np.ndarray, t: np.ndarray):
		"""存入缓存。超过容量时移除最早条目"""
		self._cache[key] = {
			'xb': xb.copy(),
			'xc': xc.copy(),
			't': t.copy(),
		}
		self._cache.move_to_end(key)
		while len(self._cache) > self._max_entries:
			self._cache.popitem(last=False)

	def clear(self):
		"""清空缓存"""
		self._cache.clear()


# ============================================================================
# RBF 代理模型
# ============================================================================

class RBFSurrogateModel:
	"""基于 RBF 的液相面代理插值模型

	使用 scipy.interpolate.RBFInterpolator 对稀疏锚点进行
	径向基函数插值，生成完整网格上的液相线温度预测。
	"""

	def __init__(self, kernel: str = 'thin_plate_spline', smoothing: float = 0.0):
		"""
		Parameters
		----------
		kernel : str
			RBF 核函数。'thin_plate_spline' 适用于光滑温度场。
		smoothing : float
			平滑参数。0=精确插值，>0=平滑拟合。
		"""
		self.kernel = kernel
		self.smoothing = smoothing
		self._interpolator: Optional[RBFInterpolator] = None
		self._anchor_xb: Optional[np.ndarray] = None
		self._anchor_xc: Optional[np.ndarray] = None
		self._anchor_t: Optional[np.ndarray] = None

	@property
	def n_anchors(self) -> int:
		"""当前锚点数"""
		return len(self._anchor_t) if self._anchor_t is not None else 0

	def fit(self, xb: np.ndarray, xc: np.ndarray, t: np.ndarray):
		"""用锚点数据拟合 RBF 插值器

		自动过滤 NaN 值。要求至少 6 个有效锚点。
		"""
		valid = ~np.isnan(t)
		xb_v, xc_v, t_v = xb[valid], xc[valid], t[valid]

		if len(t_v) < 6:
			raise ValueError(f"有效锚点不足: {len(t_v)} < 6")

		self._anchor_xb = xb_v.copy()
		self._anchor_xc = xc_v.copy()
		self._anchor_t = t_v.copy()

		coords = np.column_stack([xb_v, xc_v])
		self._interpolator = RBFInterpolator(
			coords, t_v,
			kernel=self.kernel,
			smoothing=self.smoothing,
		)

	def predict(self, xb: np.ndarray, xc: np.ndarray) -> np.ndarray:
		"""在目标网格点上预测液相线温度

		对 xb+xc > 1 的点返回 NaN。
		"""
		if self._interpolator is None:
			raise RuntimeError("模型未拟合，请先调用 fit()")

		coords = np.column_stack([xb, xc])
		t_pred = self._interpolator(coords)

		# 三元约束外的点标为 NaN
		invalid = (xb + xc) > 1.0 + 1e-9
		t_pred[invalid] = np.nan

		return t_pred

	def estimate_uncertainty(self, xb: np.ndarray, xc: np.ndarray) -> np.ndarray:
		"""估算插值不确定度（到最近锚点的距离，归一化）

		距离越大 → 不确定度越高 → 越需要细化。
		"""
		if self._anchor_xb is None:
			return np.zeros(len(xb))

		target = np.column_stack([xb, xc])
		anchor = np.column_stack([self._anchor_xb, self._anchor_xc])
		dists = cdist(target, anchor)
		min_dist = dists.min(axis=1)

		# 归一化到 [0, 1]
		d_max = min_dist.max()
		if d_max > 1e-12:
			return min_dist / d_max
		return np.zeros(len(xb))

	def compute_gradient(self, xb: np.ndarray, xc: np.ndarray,
	                     h: float = 0.005) -> np.ndarray:
		"""计算插值面的梯度幅值（有限差分）

		高梯度区域通常对应相边界，需要更多锚点。
		"""
		if self._interpolator is None:
			return np.zeros(len(xb))

		# ∂T/∂xb（中心差分）
		t_bp = self.predict(xb + h, xc)
		t_bm = self.predict(xb - h, xc)
		dtdb = (t_bp - t_bm) / (2.0 * h)

		# ∂T/∂xc（中心差分）
		t_cp = self.predict(xb, xc + h)
		t_cm = self.predict(xb, xc - h)
		dtdc = (t_cp - t_cm) / (2.0 * h)

		# 用 nan 安全方式计算幅值
		grad_mag = np.sqrt(np.nan_to_num(dtdb**2, 0) + np.nan_to_num(dtdc**2, 0))
		return grad_mag


# ============================================================================
# 自适应液相面加速计算器
# ============================================================================

class AdaptiveSurfaceCalculator:
	"""自适应液相面加速计算器

	算法流程:
	  阶段1: 在稀疏网格（6×6→~21个有效点）上精确计算液相线温度
	  阶段2: 用 RBF 插值生成完整网格预测 → 回调通知"预览"结果
	  阶段3: 识别高梯度/高不确定区域，计算额外锚点 (~50-80个)
	  阶段4: 用增强锚点集重新 RBF 插值 → 回调通知"细化"结果

	预期加速: ~5-10x（取决于网格密度）
	"""

	# ============== 算法参数 ==============
	SPARSE_GRID_N = 6           # 稀疏网格每边点数 → ~21个有效三元点
	MAX_REFINE_POINTS = 80      # 最大细化点数
	REFINE_TOP_FRACTION = 0.15  # 选择不确定度最高的前15%网格点
	GRADIENT_WEIGHT = 0.6       # 梯度在综合评分中的权重
	DISTANCE_WEIGHT = 0.4       # 距离在综合评分中的权重
	MIN_ANCHOR_POINTS = 6       # RBF 最少需要的锚点数
	DEDUP_RADIUS = 0.03         # 去重半径（距已有锚点此范围内的候选点被排除）

	def __init__(self, database, logger: Optional[Callable] = None,
	             cache: Optional[SurrogateCache] = None):
		"""
		Parameters
		----------
		database : pycalphad.Database
			PyCalphad 数据库对象
		logger : Callable, optional
			日志输出函数，签名 logger(message: str)
		cache : SurrogateCache, optional
			缓存实例。为 None 时创建新实例。
		"""
		self.db = database
		self.logger = logger or (lambda msg: None)
		self.cache = cache or SurrogateCache()
		self.surrogate = RBFSurrogateModel()

	def calculate(self,
	              comp_b_range: Tuple[float, float],
	              comp_c_range: Tuple[float, float],
	              grid_points: int,
	              T_range: Tuple[float, float],
	              components: List[str],
	              phases: Optional[List[str]] = None,
	              model_spec: Optional[Dict] = None,
	              pdens: int = 2000,
	              verbose: bool = False,
	              progress_callback: Optional[Callable] = None
	              ) -> CalculationResult:
		"""执行加速液相面计算

		Parameters
		----------
		comp_b_range, comp_c_range : Tuple[float, float]
			组分B/C的范围（通常 (0, 1)）
		grid_points : int
			目标网格每边点数
		T_range : Tuple[float, float]
			温度范围 (K)
		components : List[str]
			组分列表（含 VA）
		phases : List[str], optional
			参与计算的相
		model_spec : Dict, optional
			模型规格
		pdens : int
			pycalphad 点密度
		verbose : bool
			详细日志
		progress_callback : Callable, optional
			进度回调: callback(stage: str, result: CalculationResult, info: dict)

		Returns
		-------
		CalculationResult
			mode='surface'，结构与 LiquidusSurfaceCalculator 完全兼容
		"""
		t_start = time.time()
		callback = progress_callback or (lambda *a, **kw: None)

		try:
			# 提取三元组分名
			components = [c.upper() for c in components]
			non_va = [c for c in components if c != 'VA']
			if len(non_va) < 3:
				return self._make_error("需要至少3个非VA组分")
			comp_a, comp_b, comp_c = non_va[:3]

			# 生成完整目标网格
			full_xb, full_xc = self._generate_full_grid(grid_points)
			total_points = len(full_xb)

			if verbose:
				self.logger(f"三元系统: {comp_a}-{comp_b}-{comp_c}")
				self.logger(f"目标网格: {total_points} 个有效点")

			# 网格太小时直接全量计算（无加速意义）
			if grid_points < 8 or total_points < 30:
				self.logger("网格较小，直接使用全量计算")
				return self._fallback_calculate(
					comp_b_range, comp_c_range, grid_points,
					T_range, components, phases, model_spec, pdens, verbose)

			# ============== 阶段1: 稀疏锚点 ==============
			cache_key = self.cache.make_key(components, phases, model_spec, T_range)
			cached = self.cache.get(cache_key)

			if cached is not None:
				anchor_xb = cached['xb']
				anchor_xc = cached['xc']
				anchor_t = cached['t']
				self.logger(f"✓ 缓存命中，{len(anchor_t)} 个锚点")
			else:
				sparse_xb, sparse_xc = self._generate_sparse_grid(self.SPARSE_GRID_N)
				n_sparse = len(sparse_xb)
				self.logger(f"阶段1: 计算 {n_sparse} 个稀疏锚点...")

				anchor_t = self._compute_points_parallel(
					sparse_xb, sparse_xc, T_range, components,
					comp_a, comp_b, comp_c, phases, model_spec, pdens, verbose)
				anchor_xb = sparse_xb
				anchor_xc = sparse_xc

				# 存入缓存
				self.cache.put(cache_key, anchor_xb, anchor_xc, anchor_t)

			valid_count = int(np.sum(~np.isnan(anchor_t)))
			elapsed = time.time() - t_start

			callback(SurrogateStage.SPARSE.value, None, {
				'anchor_count': valid_count,
				'total_target': total_points,
				'elapsed_seconds': elapsed,
			})

			if valid_count < self.MIN_ANCHOR_POINTS:
				self.logger(f"⚠ 有效锚点不足 ({valid_count})，降级为全量计算")
				return self._fallback_calculate(
					comp_b_range, comp_c_range, grid_points,
					T_range, components, phases, model_spec, pdens, verbose)

			# ============== 阶段2: RBF 插值预览 ==============
			self.logger(f"阶段2: RBF 插值生成预览...")
			self.surrogate.fit(anchor_xb, anchor_xc, anchor_t)
			t_preview = self.surrogate.predict(full_xb, full_xc)

			preview_result = self._make_result(
				full_xb, full_xc, t_preview,
				comp_a, comp_b, comp_c, grid_points)

			elapsed = time.time() - t_start
			self.logger(f"✓ 预览生成 ({elapsed:.1f}s), {valid_count} 个锚点 → {total_points} 个网格点")

			callback(SurrogateStage.PREVIEW.value, preview_result, {
				'anchor_count': valid_count,
				'total_target': total_points,
				'elapsed_seconds': elapsed,
			})

			# ============== 阶段3: 自适应细化 ==============
			ref_xb, ref_xc = self._select_refinement_points(
				full_xb, full_xc, anchor_xb, anchor_xc)
			n_refine = len(ref_xb)

			if n_refine > 0:
				self.logger(f"阶段3: 计算 {n_refine} 个细化点...")

				ref_t = self._compute_points_parallel(
					ref_xb, ref_xc, T_range, components,
					comp_a, comp_b, comp_c, phases, model_spec, pdens, verbose)

				# 合并锚点集
				anchor_xb = np.concatenate([anchor_xb, ref_xb])
				anchor_xc = np.concatenate([anchor_xc, ref_xc])
				anchor_t = np.concatenate([anchor_t, ref_t])

				# 更新缓存
				self.cache.put(cache_key, anchor_xb, anchor_xc, anchor_t)

				total_anchors = int(np.sum(~np.isnan(anchor_t)))
				elapsed = time.time() - t_start

				callback(SurrogateStage.REFINING.value, None, {
					'anchor_count': total_anchors,
					'total_target': total_points,
					'elapsed_seconds': elapsed,
				})
			else:
				self.logger("阶段3: 无需细化点")
				total_anchors = valid_count

			# ============== 阶段4: 重新 RBF 插值 ==============
			self.logger(f"阶段4: 使用 {total_anchors} 个锚点重新插值...")
			self.surrogate.fit(anchor_xb, anchor_xc, anchor_t)
			t_refined = self.surrogate.predict(full_xb, full_xc)

			# 将已计算锚点的精确值覆写插值值
			t_refined = self._overlay_exact_values(
				full_xb, full_xc, t_refined,
				anchor_xb, anchor_xc, anchor_t)

			final_result = self._make_result(
				full_xb, full_xc, t_refined,
				comp_a, comp_b, comp_c, grid_points)

			elapsed = time.time() - t_start
			speedup = total_points / max(1, total_anchors)
			self.logger(f"✓ 代理计算完成: {total_anchors}/{total_points} 点, "
			            f"加速比~{speedup:.1f}x, 耗时 {elapsed:.1f}s")

			callback(SurrogateStage.REFINED.value, final_result, {
				'anchor_count': total_anchors,
				'total_target': total_points,
				'elapsed_seconds': elapsed,
				'estimated_speedup': speedup,
			})

			return final_result

		except Exception as e:
			self.logger(f"⚠ 代理加速失败 ({e})，降级为全量计算...")
			if verbose:
				self.logger(traceback.format_exc())
			return self._fallback_calculate(
				comp_b_range, comp_c_range, grid_points,
				T_range, components, phases, model_spec, pdens, verbose)

	# ====================================================================
	# 内部辅助方法
	# ====================================================================

	def _generate_sparse_grid(self, n: int = 6) -> Tuple[np.ndarray, np.ndarray]:
		"""生成稀疏三元网格 (n×n 等分，过滤 xb+xc<=1)"""
		b_vals = np.linspace(0, 1, n)
		c_vals = np.linspace(0, 1, n)
		bb, cc = np.meshgrid(b_vals, c_vals)
		mask = bb + cc <= 1.0 + 1e-9
		return bb[mask].ravel(), cc[mask].ravel()

	def _generate_full_grid(self, grid_points: int) -> Tuple[np.ndarray, np.ndarray]:
		"""生成完整三元网格（与 LiquidusSurfaceCalculator 兼容）"""
		step = 1.0 / (grid_points - 1)
		b_coords, c_coords = np.meshgrid(
			np.arange(0, 1 + step, step),
			np.arange(0, 1 + step, step)
		)
		mask = b_coords + c_coords <= 1.0001
		return b_coords[mask].ravel(), c_coords[mask].ravel()

	def _compute_points_parallel(self,
	                             xb_points: np.ndarray,
	                             xc_points: np.ndarray,
	                             T_range: Tuple[float, float],
	                             components: List[str],
	                             comp_a: str, comp_b: str, comp_c: str,
	                             phases: Optional[List[str]],
	                             model_spec,
	                             pdens: int,
	                             verbose: bool = False
	                             ) -> np.ndarray:
		"""并行计算一批网格点的液相线温度

		复用 thermocal_core 的 _surface_worker_init / _surface_calculate_point
		"""
		n_points = len(xb_points)
		t_result = np.full(n_points, np.nan)

		# 构建与原始格式一致的 task 列表
		tasks = []
		for i in range(n_points):
			xb = float(xb_points[i])
			xc = float(xc_points[i])
			xa = 1.0 - xb - xc
			composition = self._build_composition(comp_a, xa, comp_b, xb, comp_c, xc)
			tasks.append((i, composition, T_range, components, phases,
			              model_spec, pdens))

		# 多进程并行
		n_workers = max(1, os.cpu_count() - 1) if os.cpu_count() else 1
		use_parallel = n_points >= 4 and n_workers > 1

		if use_parallel:
			try:
				db_bytes = pickle.dumps(self.db)
				chunksize = max(1, n_points // (n_workers * 4))

				with Pool(n_workers, initializer=_surface_worker_init,
				          initargs=(db_bytes,)) as pool:
					results = pool.map(_surface_calculate_point, tasks,
					                   chunksize=chunksize)

				for idx, t_liq in results:
					if t_liq is not None and not np.isnan(t_liq):
						t_result[idx] = t_liq

			except Exception as e:
				self.logger(f"  并行失败: {e}, 降级为串行")
				use_parallel = False

		if not use_parallel:
			# 串行回退
			_worker_state['db'] = self.db
			_worker_state['ls_calc'] = LiquidusSolidusCalculator(self.db)
			for task in tasks:
				try:
					idx, t_liq = _surface_calculate_point(task)
					if t_liq is not None and not np.isnan(t_liq):
						t_result[idx] = t_liq
				except Exception:
					pass

		return t_result

	def _select_refinement_points(self,
	                               full_xb: np.ndarray,
	                               full_xc: np.ndarray,
	                               anchor_xb: np.ndarray,
	                               anchor_xc: np.ndarray
	                               ) -> Tuple[np.ndarray, np.ndarray]:
		"""选择需要细化的网格点

		综合考虑 RBF 梯度（相边界）+ 到锚点距离（覆盖不足区域）。
		"""
		# 1. 计算梯度和距离
		grad = self.surrogate.compute_gradient(full_xb, full_xc)
		dist = self.surrogate.estimate_uncertainty(full_xb, full_xc)

		# 2. 归一化后加权
		g_min, g_max = np.nanmin(grad), np.nanmax(grad)
		d_min, d_max = np.nanmin(dist), np.nanmax(dist)

		grad_norm = (grad - g_min) / (g_max - g_min + 1e-12)
		dist_norm = (dist - d_min) / (d_max - d_min + 1e-12)

		score = self.GRADIENT_WEIGHT * grad_norm + self.DISTANCE_WEIGHT * dist_norm

		# 3. 去除 NaN 评分的点
		valid_mask = ~np.isnan(score)
		indices = np.where(valid_mask)[0]
		if len(indices) == 0:
			return np.array([]), np.array([])

		# 4. 排除已有锚点附近的候选点
		target_coords = np.column_stack([full_xb[indices], full_xc[indices]])
		anchor_coords = np.column_stack([anchor_xb, anchor_xc])
		min_dists = cdist(target_coords, anchor_coords).min(axis=1)
		far_mask = min_dists > self.DEDUP_RADIUS
		indices = indices[far_mask]

		if len(indices) == 0:
			return np.array([]), np.array([])

		# 5. 取评分最高的 N 个
		max_pts = min(self.MAX_REFINE_POINTS,
		              int(self.REFINE_TOP_FRACTION * len(full_xb)))
		max_pts = max(max_pts, 10)  # 至少选10个

		scores_filtered = score[indices]
		top_k = min(max_pts, len(indices))
		top_idx = np.argsort(scores_filtered)[::-1][:top_k]

		selected = indices[top_idx]
		return full_xb[selected], full_xc[selected]

	def _overlay_exact_values(self,
	                          full_xb: np.ndarray, full_xc: np.ndarray,
	                          t_interp: np.ndarray,
	                          anchor_xb: np.ndarray, anchor_xc: np.ndarray,
	                          anchor_t: np.ndarray,
	                          tol: float = 1e-6) -> np.ndarray:
		"""将已精确计算的锚点值覆写到插值结果中"""
		result = t_interp.copy()
		for i in range(len(anchor_xb)):
			if np.isnan(anchor_t[i]):
				continue
			# 找到最近的网格点
			d = np.sqrt((full_xb - anchor_xb[i])**2 +
			            (full_xc - anchor_xc[i])**2)
			nearest = np.argmin(d)
			if d[nearest] < tol:
				result[nearest] = anchor_t[i]
		return result

	@staticmethod
	def _build_composition(comp_a: str, xa: float,
	                       comp_b: str, xb: float,
	                       comp_c: str, xc: float,
	                       min_fraction: float = 1e-9) -> Dict:
		"""构建三元成分字典（与 LiquidusSurfaceCalculator 一致）"""
		composition = {}
		if xa > min_fraction:
			composition[comp_a] = xa
		if xb > min_fraction:
			composition[comp_b] = xb
		if xc > min_fraction:
			composition[comp_c] = xc
		if not composition:
			composition = {comp_a: 1.0}
		return composition

	def _make_result(self, xb: np.ndarray, xc: np.ndarray, t_data: np.ndarray,
	                 comp_a: str, comp_b: str, comp_c: str,
	                 grid_points: int) -> CalculationResult:
		"""构造标准 CalculationResult（与 LiquidusSurfaceCalculator 格式一致）"""
		valid_count = int(np.sum(~np.isnan(t_data)))
		return CalculationResult(
			success=True,
			mode=CalculationMode.SURFACE.value,
			message=f"液相面代理计算完成: {valid_count}/{len(t_data)}有效点",
			x_axis=xb,
			y_axis=xc,
			z_axis=t_data,
			surface_data={
				'xb_data': xb,
				'xc_data': xc,
				't_data': t_data,
				'comp_a': comp_a,
				'comp_b': comp_b,
				'comp_c': comp_c,
				'grid_points': grid_points,
			}
		)

	@staticmethod
	def _make_error(message: str) -> CalculationResult:
		"""构造失败结果"""
		return CalculationResult(
			success=False,
			mode=CalculationMode.SURFACE.value,
			message=message,
		)

	def _fallback_calculate(self, *args, **kwargs) -> CalculationResult:
		"""降级为原始全量计算"""
		self.logger("降级: 使用 LiquidusSurfaceCalculator 全量计算")
		fallback = LiquidusSurfaceCalculator(self.db, self.logger)
		return fallback.calculate(*args, **kwargs)


# ============================================================================
# 导出
# ============================================================================

__all__ = [
	'SurrogateCache',
	'RBFSurrogateModel',
	'AdaptiveSurfaceCalculator',
	'SurrogateStage',
]
