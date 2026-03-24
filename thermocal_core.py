#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyCalphad 核心计算模块 (Refactored v2.0)
==========================================

统一的后端核心模块，包含所有计算功能：
- 数据库管理
- 液相线/固相线计算
- 热力学性质计算
- 溶解度计算
- 相图计算

设计理念：
1. 计算器初始化时只需database和logger
2. 具体计算参数在调用calculate方法时传入
3. 所有Calculator统一返回CalculationResult对象
4. 严禁抛出未捕获的异常
"""

import os
import re
import traceback
import pickle
import time
import numpy as np
from enum import Enum
from typing import Dict, List, Tuple, Optional, Callable, Any, Union
from dataclasses import dataclass, field
from pycalphad import Database, equilibrium, variables as v, Model
from pycalphad.plot.binary import binplot
import matplotlib
from collections import defaultdict
matplotlib.use('Agg')  # 线程安全的后端
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

# 导入高级模型（如果可用）
try:
	from pycalphad.models.model_uem import ModelMuggianu, ModelToop, ModelUEM as ModelUEM1, ModelWithUEM

	ADVANCED_MODELS_AVAILABLE = True
except ImportError:
	ADVANCED_MODELS_AVAILABLE = False
	ModelMuggianu = None
	ModelToop = None
	ModelUEM1 = None
	ModelWithUEM = None


# ============================================================================
# 统一数据契约 (Data Contract)
# ============================================================================

class CalculationMode(Enum):
	"""计算模式枚举"""
	BINPLOT = 'binplot'  # 真二元相图（PyCalphad binplot）
	PHASE_MAP = 'phase_map'  # 伪二元相区图
	LINE = 'line'  # 线性数据（液相线/固相线曲线）
	TERNARY = 'ternary'  # 三元相图
	SURFACE = 'surface'  # 液相面投影
	SOLUBILITY = 'solubility'  # 溶解度计算
	PROPERTIES = 'properties'  # 热力学性质
	SOLVUS = 'solvus'  # 溶解度曲线


@dataclass
class CalculationResult:
	"""
	统一的计算结果数据契约

	所有Calculator的calculate方法必须返回此类型
	"""
	success: bool
	mode: str  # 使用字符串以便JSON序列化，值应为CalculationMode的value
	message: str = ""
	
	# ========== 通用数据载荷 (Payload) ==========
	figure: Optional[Figure] = None  # 用于 binplot / ternary
	x_axis: Optional[np.ndarray] = None  # X轴数据（成分、温度等）
	y_axis: Optional[np.ndarray] = None  # Y轴数据
	z_axis: Optional[np.ndarray] = None  # Z轴数据（3D图）
	
	# ========== 相图专用 ==========
	phase_map: Optional[np.ndarray] = None  # 相区图矩阵 (T x X)，int32编码
	phase_legend: Optional[Dict[int, str]] = None  # 相区图编码映射: int -> phase_name
	t_grid: Optional[np.ndarray] = None  # 温度网格
	x_range: Optional[Tuple[float, float]] = None  # X轴范围
	t_range: Optional[Tuple[float, float]] = None  # 温度范围
	components: Optional[List[str]] = None  # 组分列表
	
	# ========== 液相线/固相线专用 ==========
	liquidus_data: Optional[Dict[str, Any]] = None
	# 结构: {
	#     'liquidus_K': float, 'solidus_K': float,
	#     'liquidus_C': float, 'solidus_C': float,
	#     'range_K': float, 'range_C': float,
	#     'temperatures': np.ndarray, 'liquid_fraction': np.ndarray,
	#     'composition': Dict
	# }
	
	# ========== 溶解度专用 ==========
	solubility_data: Optional[Dict[str, Any]] = None
	# 结构: {
	#     'temperature_K': float, 'solubility': float,
	#     'matrix_phase': str, 'composition': Dict,
	#     # 或用于solvus曲线:
	#     'temperatures': np.ndarray, 'solubilities': np.ndarray,
	#     'phases': List[str]
	# }
	
	# ========== 热力学性质专用 ==========
	properties_data: Optional[Dict[str, Any]] = None
	# 结构: {
	#     'compositions': List[Dict],
	#     'temperature': float,
	#     'gibbs_energy': np.ndarray,      # GM
	#     'enthalpy_mix': np.ndarray,      # HM_MIX
	#     'activity': Dict[str, np.ndarray],
	#     'x_values': np.ndarray           # 变化组分的摩尔分数
	# }
	
	# ========== 液相面专用 ==========
	surface_data: Optional[Dict[str, Any]] = None
	
	# 结构: {
	#     'xb_data': np.ndarray, 'xc_data': np.ndarray,
	#     't_data': np.ndarray,
	#     'comp_a': str, 'comp_b': str, 'comp_c': str
	# }
	
	def to_dict (self) -> Dict[str, Any]:
		"""转换为字典（用于JSON序列化）"""
		result = {
			'success': self.success,
			'mode': self.mode,
			'message': self.message
		}
		
		# 添加非None的数据字段
		if self.x_axis is not None:
			result['x_axis'] = self.x_axis.tolist() if isinstance(self.x_axis, np.ndarray) else self.x_axis
		if self.y_axis is not None:
			result['y_axis'] = self.y_axis.tolist() if isinstance(self.y_axis, np.ndarray) else self.y_axis
		if self.z_axis is not None:
			result['z_axis'] = self.z_axis.tolist() if isinstance(self.z_axis, np.ndarray) else self.z_axis
		if self.phase_map is not None:
			result['phase_map'] = self.phase_map.tolist() if isinstance(self.phase_map, np.ndarray) else self.phase_map
		if self.phase_legend is not None:
			result['phase_legend'] = {str(k): v for k, v in self.phase_legend.items()}
		if self.t_grid is not None:
			result['t_grid'] = self.t_grid.tolist() if isinstance(self.t_grid, np.ndarray) else self.t_grid
		if self.liquidus_data is not None:
			result['liquidus_data'] = self._serialize_dict(self.liquidus_data)
		if self.solubility_data is not None:
			result['solubility_data'] = self._serialize_dict(self.solubility_data)
		if self.properties_data is not None:
			result['properties_data'] = self._serialize_dict(self.properties_data)
		if self.surface_data is not None:
			result['surface_data'] = self._serialize_dict(self.surface_data)
		
		return result
	
	@staticmethod
	def _serialize_dict (d: Dict) -> Dict:
		"""递归序列化字典中的numpy数组"""
		result = {}
		for k, v in d.items():
			if isinstance(v, np.ndarray):
				result[k] = v.tolist()
			elif isinstance(v, dict):
				result[k] = CalculationResult._serialize_dict(v)
			elif isinstance(v, (np.floating, np.integer)):
				result[k] = float(v) if isinstance(v, np.floating) else int(v)
			else:
				result[k] = v
		return result


# ============================================================================
# 贡献系数打印混入类（消除重复代码）
# ============================================================================

class _ContributionCoefficientMixin:
	"""提供 _get_model_name_from_spec 和 _print_contribution_coefficients 方法的混入类"""

	@staticmethod
	def _get_model_name_from_spec(model_spec):
		"""从 model_spec 中提取模型名称字符串

		model_spec 可以是:
		1. None (RKM 默认模型)
		2. 模型类 (例如 ModelToop, ModelMuggianu)
		3. 字符串 (模型名称)
		4. 字典 (相到模型的映射)
		"""
		if model_spec is None:
			return None

		if isinstance(model_spec, str):
			return model_spec

		if isinstance(model_spec, dict):
			liquid_model = model_spec.get('LIQUID')
			if liquid_model and hasattr(liquid_model, '__name__'):
				class_name = liquid_model.__name__
			else:
				return None
		elif hasattr(model_spec, '__name__'):
			class_name = model_spec.__name__
		else:
			return None

		class_to_name = {
			'ModelUEM1': 'uem1',
			'ModelUEMAdv': 'uem_adv',
			'ModelUEM2N': 'uem2_n',
			'ModelToop': 'toop',
			'ModelMuggianu': 'muggianu',
			'ModelKohler': 'kohler',
			'Model': None
		}

		return class_to_name.get(class_name, None)

	def _print_contribution_coefficients(self, components: List[str],
	                                     temperature: float,
	                                     model_spec):
		"""打印贡献系数

		Parameters
		----------
		components : List[str]
			组分列表（不含 VA）
		temperature : float
			温度（K）
		model_spec :
			模型规格
		"""
		try:
			model_name = self._get_model_name_from_spec(model_spec)

			if not model_name or model_name.lower() not in ['uem1', 'uem_adv', 'uem2_n', 'toop', 'muggianu']:
				return

			from pycalphad.models.model_uem import ModelUEM as ModelUEM1, ModelToop, ModelMuggianu, ModelWithUEM

			components_with_va = sorted([c.upper() for c in components]) + ['VA']

			model_classes = {
				'uem1': ModelUEM1,
				'uem_adv': ModelUEM1,
				'uem2_n': ModelUEM1,
				'toop': ModelToop,
				'muggianu': ModelMuggianu
			}

			model_class = model_classes.get(model_name.lower())
			if model_class is None:
				return

			model = model_class(self.db, components_with_va, 'LIQUID')

			active_components = [c for c in components_with_va if c != 'VA']
			alpha_coeffs = model.get_contribution_coefficients(
				active_components, subl_idx=0, temperature=temperature,
				logger=self.logger
			)

			if alpha_coeffs:
				self.logger(f"  贡献系数 α @ {temperature:.0f} K ({model_name.upper()}模型):")
				sorted_keys = sorted(alpha_coeffs.keys())
				for k, i, j in sorted_keys:
					alpha_ki, alpha_kj = alpha_coeffs[(k, i, j)]
					self.logger(f"    α({k},{i}|{j})={alpha_ki:.4f}  α({k},{j}|{i})={alpha_kj:.4f}")
			else:
				self.logger(f"  贡献系数 @ {temperature:.0f}K: 无系数")

		except Exception as e:
			self.logger(f"  贡献系数计算失败: {e}")


# ============================================================================
# 数据库管理类
# ============================================================================

class DatabaseManager:
	"""数据库管理器 - 负责数据库加载和管理"""
	
	def __init__ (self):
		self.dbe = None
		self.available_comps = []
		self.available_phases = []
		self.pdens = 2000
		
		# 可用模型类字典
		self.available_models_cls = {'RKM': Model}
		if ADVANCED_MODELS_AVAILABLE:
			self.available_models_cls.update({
				'Muggianu': ModelMuggianu,
				'Toop': ModelToop,
				'UEM1': ModelUEM1
			})
	
	def set_pdens (self, val: int):
		"""设置点密度"""
		self.pdens = int(val)
	
	def load_database (self, file_paths: List[str]) -> List[str]:
		"""
		加载TDB数据库文件

		Parameters:
		-----------
		file_paths : list
			TDB文件路径列表

		Returns:
		--------
		list : 成功加载的文件名列表
		"""
		if not file_paths:
			return None
		
		all_tdb_content = ""
		loaded_files_list = []
		
		for path in file_paths:
			try:
				with open(path, 'r', encoding='latin-1') as f:
					all_tdb_content += f.read() + "\n"
				loaded_files_list.append(os.path.basename(path))
			except Exception:
				try:
					with open(path, 'r', encoding='utf-8') as f:
						all_tdb_content += f.read() + "\n"
					loaded_files_list.append(os.path.basename(path))
				except Exception as e:
					raise IOError(f"无法读取文件 {os.path.basename(path)}: {e}")
		
		self.dbe = Database(all_tdb_content)
		self._update_available_data()
		return loaded_files_list
	
	def _update_available_data (self):
		"""更新可用组分和相列表"""
		if self.dbe is None:
			return
		
		# 提取组分
		all_elements = []
		for element in self.dbe.elements:
			elem_str = str(element).strip().upper()
			if elem_str not in ['VA', '/-']:
				if elem_str:
					all_elements.append(elem_str)
		
		# 标准化组分名称
		self.available_comps = []
		for elem in all_elements:
			if len(elem) <= 2 and elem[0].isalpha():
				if len(elem) == 1:
					standardized = elem.upper()
				else:
					standardized = elem[0].upper() + elem[1].lower()
				self.available_comps.append(standardized)
		
		self.available_comps = sorted(list(set(self.available_comps)))
		self.available_phases = sorted(list(self.dbe.phases.keys()))
	
	def get_model_spec (self, model_key: str, uem1_liquid_only: bool = False):
		"""
		获取模型规格

		Parameters:
		-----------
		model_key : str
			模型键名 ('RKM', 'Muggianu', 'Toop', 'UEM1')
		uem1_liquid_only : bool
			UEM1模型是否仅应用于液相

		Returns:
		--------
		模型规格 (None, Model类, 或模型字典)
		"""
		if model_key == 'RKM':
			return None
		
		model_class = self.available_models_cls.get(model_key)
		if not model_class:
			return None
		
		# 溶液相列表
		solution_phases = ['LIQUID', 'FCC_A1', 'BCC_A2', 'HCP_A3', 'HCP_ZN', 'DIAMOND_A4']
		
		model_dict = {}
		for phase_name in self.available_phases:
			phase_upper = phase_name.upper()
			use_custom = False
			
			if model_key == 'UEM1' and uem1_liquid_only:
				if 'LIQUID' in phase_upper:
					use_custom = True
			else:
				for sol_ph in solution_phases:
					if (phase_upper == sol_ph or
							phase_upper.startswith(sol_ph + '#') or
							phase_upper.startswith(sol_ph + ':')):
						use_custom = True
						break
			
			model_dict[phase_name] = model_class if use_custom else Model
		
		return model_dict
	
	@staticmethod
	def parse_base_alloy (alloy_str: str, available_comps: List[str] = None) -> Dict[str, float]:
		"""
		解析基础合金字符串

		Parameters:
		-----------
		alloy_str : str
			合金字符串，如 "AL", "AL1CU1", "ALSI0.2"
		available_comps : list, optional
			可用组分列表，用于验证

		Returns:
		--------
		dict : 组分比例字典，如 {'AL': 1.0, 'CU': 1.0}
		"""
		alloy_str_upper = alloy_str.strip().upper().replace(" ", "")
		if not alloy_str_upper:
			raise ValueError("基础合金字符串不能为空")
		
		# 匹配 (AL)(1), (SI)(0.2)
		matches = re.findall(r'([A-Z]{1,2})(\d*\.?\d*)', alloy_str_upper)
		
		if not matches:
			raise ValueError(f"无法解析基础合金: '{alloy_str}'")
		
		# 验证解析完整性
		parsed_length = sum(len(m[0]) + len(m[1]) for m in matches)
		if parsed_length != len(alloy_str_upper):
			unparsed_part = alloy_str_upper[parsed_length:]
			raise ValueError(f"基础合金格式错误: '{alloy_str}'. 无法解析 '{unparsed_part}'")
		
		base_comps_dict = {}
		available_comps_upper = [c.upper() for c in (available_comps or [])]
		
		for comp, ratio_str in matches:
			# 验证组分是否在数据库中
			if available_comps and comp not in available_comps_upper:
				if comp == 'VA':
					continue
				raise ValueError(f"基础合金中的组分 '{comp}' 在数据库中不存在！")
			
			# 解析比例
			if ratio_str == "":
				ratio = 1.0
			else:
				try:
					ratio = float(ratio_str)
				except ValueError:
					raise ValueError(f"基础合金中 '{comp}' 的比例 '{ratio_str}' 无效")
			
			if ratio <= 0:
				raise ValueError(f"基础合金中 '{comp}' 的比例必须大于 0")
			
			base_comps_dict[comp] = base_comps_dict.get(comp, 0) + ratio
		
		if not base_comps_dict:
			raise ValueError("基础合金未包含任何有效组分。")
		
		return base_comps_dict


# ============================================================================
# 液相线/固相线计算器
# ============================================================================

class LiquidusSolidusCalculator:
	"""
	液相线/固相线温度计算器 - 自适应网格细化版本

	初始化时只需database和logger，计算时传入具体参数
	返回统一的CalculationResult对象

	改进点：
	- 使用自适应网格细化替代暴力线性扫描
	- 在相变区自动加密，平坦区自动稀疏
	- 计算效率提升数倍，同时保证绘图精度
	"""
	
	# ============== 自适应算法参数（类级别常量）==============
	AMR_MAX_DEPTH = 12  # 最大递归深度（对应最小区间 ~0.024%）
	AMR_CURVATURE_TOL = 0.02  # 曲率容差：中点偏离线性插值的阈值
	AMR_MIN_INTERVAL = 0.5  # 最小温度区间 (K)，防止过度细分
	AMR_TRANSITION_ZONE_TOL = 0.05  # 两相区强制细分阈值
	
	def __init__ (self, database, logger: Optional[Callable] = None):
		"""
		初始化液相线/固相线计算器

		Parameters:
		-----------
		database : Database
			PyCalphad数据库对象
		logger : Callable, optional
			日志函数
		"""
		self.db = database
		self.logger = logger or (lambda msg: None)
	
	def calculate (self,
	               composition: Dict,
	               T_range: Tuple[float, float],
	               components: List[str],
	               phases: Optional[List[str]] = None,
	               model_spec: Optional[Dict] = None,
	               pdens: int = 2000,
	               T_step: float = 5.0,
	               refine: bool = True,
	               use_bisection: bool = True,
	               bisection_tol: float = 1.0,
	               verbose: bool = False) -> CalculationResult:
		"""
		计算单个成分的液相线和固相线温度

		Parameters:
		-----------
		composition : dict
			成分字典，如 {'AL': 0.7, 'CU': 0.3}
		T_range : tuple
			温度范围 (T_min, T_max)，单位K
		components : List[str]
			组分列表（包含VA）
		phases : List[str], optional
			相列表。如果为None，使用数据库中所有相
		model_spec : Dict, optional
			模型规范（如UEM模型字典）
		pdens : int
			点密度参数
		T_step : float
			温度步长，单位K（仅用于兼容性，自适应模式下作为参考）
		refine : bool
			是否进行精修
		use_bisection : bool
			是否使用二分法（True=快速模式，False=自适应扫描模式）
		bisection_tol : float
			二分法精度
		verbose : bool
			是否输出详细信息

		Returns:
		--------
		CalculationResult with mode='line'
		"""
		try:
			# 参数预处理
			composition = {k.upper(): val for k, val in composition.items()}
			components = [c.upper() for c in components]
			T_min, T_max = T_range
			active_phases = phases if phases else list(self.db.phases.keys())
			
			if verbose:
				method = "二分法（快速模式）" if use_bisection else "自适应网格细化"
				self.logger(f"计算成分: {composition}, 温度范围: {T_min}-{T_max}K, 方法: {method}")
			
			# ============== 模式选择 ==============
			if use_bisection:
				# 快速二分法 + 绘图数据生成（向量化优化）
				T_sol, T_liq, temperatures, liquid_fraction = self._bisection_with_optional_scan(
						composition, T_range, components, active_phases,
						model_spec, pdens, bisection_tol, T_step, verbose
				)
			else:
				# 完整模式：自适应网格细化扫描
				temperatures, liquid_fraction, T_sol, T_liq = self._adaptive_scan(
						composition, T_range, components, active_phases,
						model_spec, pdens, T_step, verbose
				)
			
			# ============== 结果封装 ==============
			T_liq_C = T_liq - 273.15 if T_liq is not None else None
			T_sol_C = T_sol - 273.15 if T_sol is not None else None
			range_K = T_liq - T_sol if (T_liq and T_sol) else None
			range_C = range_K if range_K is not None else None
			
			liquidus_data = {
				'liquidus_K': T_liq,
				'solidus_K': T_sol,
				'liquidus_C': T_liq_C,
				'solidus_C': T_sol_C,
				'range_K': range_K,
				'range_C': range_C,
				'temperatures': temperatures,
				'liquid_fraction': liquid_fraction,
				'composition': composition,
				'calculation_points': len(temperatures)  # 诊断信息
			}
			
			msg = ""
			if T_liq and T_sol:
				msg = f"液相线: {T_liq:.1f}K, 固相线: {T_sol:.1f}K (计算点数: {len(temperatures)})"
			else:
				msg = f"计算完成 (计算点数: {len(temperatures)})"
			
			return CalculationResult(
					success=True,
					mode=CalculationMode.LINE.value,
					message=msg,
					x_axis=temperatures,
					y_axis=liquid_fraction if len(liquid_fraction) > 0 else None,
					liquidus_data=liquidus_data
			)
		
		except Exception as e:
			if verbose:
				self.logger(traceback.format_exc())
			return CalculationResult(
					success=False,
					mode=CalculationMode.LINE.value,
					message=f"计算失败: {str(e)}",
					liquidus_data={'composition': composition, 'error': str(e)}
			)
	
	def calculate_batch (self,
	                     compositions: List[Dict],
	                     T_range: Tuple[float, float],
	                     components: List[str],
	                     phases: Optional[List[str]] = None,
	                     model_spec: Optional[Dict] = None,
	                     pdens: int = 2000,
	                     T_step: float = 5.0,
	                     refine: bool = True,
	                     use_bisection: bool = True,
	                     bisection_tol: float = 1.0,
	                     verbose: bool = False) -> CalculationResult:
		"""批量计算多个成分的液相线和固相线，返回统一结果"""
		try:
			results = []
			total = len(compositions)
			total_points = 0
			
			for i, comp in enumerate(compositions):
				if verbose and i % max(1, total // 10) == 0:
					self.logger(f"进度: {i}/{total}")
				
				result = self.calculate(
						comp, T_range, components, phases, model_spec,
						pdens, T_step, refine, use_bisection, bisection_tol,
						verbose=False
				)
				
				if result.liquidus_data:
					results.append(result.liquidus_data)
					total_points += result.liquidus_data.get('calculation_points', 0)
				
				if verbose and result.success:
					data = result.liquidus_data
					comp_str = ', '.join([f'{k}={v:.3f}' for k, v in comp.items()])
					if data.get('liquidus_K') and data.get('solidus_K'):
						self.logger(
								f"  [{i + 1}/{total}] {comp_str}: "
								f"液相线={data['liquidus_K']:.1f}K, "
								f"固相线={data['solidus_K']:.1f}K, "
								f"点数={data.get('calculation_points', 'N/A')}"
						)
			
			return CalculationResult(
					success=True,
					mode=CalculationMode.LINE.value,
					message=f"批量计算完成: {len(results)}/{total} 成功, 总计算点数: {total_points}",
					liquidus_data={'batch_results': results, 'count': len(results), 'total_points': total_points}
			)
		
		except Exception as e:
			return CalculationResult(
					success=False,
					mode=CalculationMode.LINE.value,
					message=f"批量计算失败: {str(e)}"
			)
	
	# ============== 核心算法：自适应网格细化 ==============
	
	def _adaptive_scan (self,
	                    composition: Dict,
	                    T_range: Tuple[float, float],
	                    components: List[str],
	                    phases: List[str],
	                    model_spec: Optional[Dict],
	                    pdens: int,
	                    T_step: float,
	                    verbose: bool) -> Tuple[np.ndarray, np.ndarray, Optional[float], Optional[float]]:
		"""
		自适应网格细化扫描算法

		核心思想：
		1. 从粗网格开始，计算端点和中点
		2. 根据"细分判据"决定是否需要进一步细化区间
		3. 细分判据包括：
		   - 曲率判据：中点fL偏离线性插值过大
		   - 相变判据：区间跨越fL=0（固相线）或fL=1（液相线）
		   - 两相区判据：在0<fL<1区域保持足够分辨率

		Parameters:
		-----------
		composition : dict
			成分字典
		T_range : tuple
			温度范围 (T_min, T_max)
		components : list
			组分列表
		phases : list
			相列表
		model_spec : dict, optional
			模型规范
		pdens : int
			点密度参数
		T_step : float
			参考温度步长（用于确定初始粗网格）
		verbose : bool
			详细输出

		Returns:
		--------
		tuple: (temperatures, liquid_fractions, T_sol, T_liq)
		"""
		T_min, T_max = T_range
		
		# 计算参数打包，便于传递
		calc_params = (composition, components, phases, model_spec, pdens)
		
		# 结果存储：使用字典避免重复计算
		computed_points: Dict[float, float] = {}
		
		def get_fL (T: float) -> float:
			"""获取某温度的液相分数（带缓存）"""
			# 四舍五入到0.1K精度，避免浮点数问题
			T_key = round(T, 1)
			if T_key not in computed_points:
				fL = self._calculate_liquid_fraction_at_T(T, *calc_params)
				computed_points[T_key] = fL if fL >= 0 else 0.0  # 处理计算失败情况
			return computed_points[T_key]
		
		# 初始化：计算端点
		fL_min = get_fL(T_min)
		fL_max = get_fL(T_max)
		
		if verbose:
			self.logger(f"  初始端点: T={T_min:.1f}K, fL={fL_min:.4f}; T={T_max:.1f}K, fL={fL_max:.4f}")
		
		# 使用栈实现非递归的自适应细分（避免Python递归深度限制）
		# 栈元素: (T_low, T_high, fL_low, fL_high, depth)
		stack = [(T_min, T_max, fL_min, fL_max, 0)]
		
		while stack:
			T_low, T_high, fL_low, fL_high, depth = stack.pop()
			
			# 达到最大深度，停止细分
			if depth >= self.AMR_MAX_DEPTH:
				continue
			
			# 区间太小，停止细分
			interval = T_high - T_low
			if interval < self.AMR_MIN_INTERVAL:
				continue
			
			# 计算中点
			T_mid = (T_low + T_high) / 2.0
			fL_mid = get_fL(T_mid)
			
			# ============== 细分判据 ==============
			need_refine = self._should_refine(
					T_low, T_high, T_mid,
					fL_low, fL_high, fL_mid,
					interval, depth
			)
			
			if need_refine:
				# 将两个子区间压入栈（注意顺序，先压右区间，后压左区间，保证左区间先处理）
				stack.append((T_mid, T_high, fL_mid, fL_high, depth + 1))
				stack.append((T_low, T_mid, fL_low, fL_mid, depth + 1))
		
		# 提取并排序结果
		sorted_temps = sorted(computed_points.keys())
		temperatures = np.array(sorted_temps)
		liquid_fractions = np.array([computed_points[T] for T in sorted_temps])
		
		# 从扫描结果中提取相变温度
		T_sol, T_liq = self._extract_boundaries_from_scan(temperatures, liquid_fractions)
		
		if verbose:
			self.logger(
				f"  自适应扫描完成: {len(temperatures)} 点 (vs 均匀网格 ~{int((T_max - T_min) / T_step) + 1} 点)")
			if T_sol and T_liq:
				self.logger(f"  相变温度: T_sol={T_sol:.1f}K, T_liq={T_liq:.1f}K")
		
		return temperatures, liquid_fractions, T_sol, T_liq
	
	def _should_refine (self,
	                    T_low: float, T_high: float, T_mid: float,
	                    fL_low: float, fL_high: float, fL_mid: float,
	                    interval: float, depth: int) -> bool:
		"""
		判断是否需要继续细分区间

		细分判据（满足任一即细分）：
		1. 曲率判据：中点fL与线性插值的偏差超过阈值
		2. 固相线判据：区间跨越fL≈0的边界
		3. 液相线判据：区间跨越fL≈1的边界
		4. 两相区判据：在两相共存区保持足够分辨率

		Returns:
		--------
		bool: True表示需要继续细分
		"""
		# 阈值定义
		SOLIDUS_THRESHOLD = 0.005  # 固相线判定阈值
		LIQUIDUS_THRESHOLD = 0.995  # 液相线判定阈值
		
		# 判据1：曲率判据
		# 计算中点的线性插值预测值
		fL_interp = (fL_low + fL_high) / 2.0
		curvature_error = abs(fL_mid - fL_interp)
		if curvature_error > self.AMR_CURVATURE_TOL:
			return True
		
		# 判据2：固相线判据（跨越fL≈0边界）
		# 检查是否有一侧是纯固相，另一侧出现液相
		crosses_solidus = (
				(fL_low <= SOLIDUS_THRESHOLD < fL_mid) or
				(fL_mid <= SOLIDUS_THRESHOLD < fL_high) or
				(fL_low > SOLIDUS_THRESHOLD >= fL_mid) or
				(fL_mid > SOLIDUS_THRESHOLD >= fL_high)
		)
		if crosses_solidus:
			return True
		
		# 判据3：液相线判据（跨越fL≈1边界）
		crosses_liquidus = (
				(fL_low < LIQUIDUS_THRESHOLD <= fL_mid) or
				(fL_mid < LIQUIDUS_THRESHOLD <= fL_high) or
				(fL_low >= LIQUIDUS_THRESHOLD > fL_mid) or
				(fL_mid >= LIQUIDUS_THRESHOLD > fL_high)
		)
		if crosses_liquidus:
			return True
		
		# 判据4：两相区判据
		# 在两相共存区（0 < fL < 1），如果区间仍然较大，继续细分以确保绘图分辨率
		in_transition = (
				(SOLIDUS_THRESHOLD < fL_low < LIQUIDUS_THRESHOLD) or
				(SOLIDUS_THRESHOLD < fL_mid < LIQUIDUS_THRESHOLD) or
				(SOLIDUS_THRESHOLD < fL_high < LIQUIDUS_THRESHOLD)
		)
		if in_transition:
			# 动态调整两相区的最小分辨率（基于深度）
			# 深度越深，允许的最小区间越小
			min_transition_interval = max(2.0, 50.0 / (2 ** (depth // 2)))
			if interval > min_transition_interval:
				return True
		
		return False
	
	def _extract_boundaries_from_scan (self,
	                                   temperatures: np.ndarray,
	                                   liquid_fractions: np.ndarray) -> Tuple[Optional[float], Optional[float]]:
		"""
		从自适应扫描结果中提取液相线和固相线温度

		使用线性插值精确定位fL=0.001和fL=0.999的位置

		Returns:
		--------
		tuple: (T_solidus, T_liquidus)
		"""
		T_sol = None
		T_liq = None
		
		n = len(temperatures)
		if n < 2:
			return None, None
		
		# 固相线：找到fL从≤0.001变为>0.001的位置
		for i in range(n - 1):
			if liquid_fractions[i] <= 0.001 < liquid_fractions[i + 1]:
				# 线性插值
				T1, T2 = temperatures[i], temperatures[i + 1]
				f1, f2 = liquid_fractions[i], liquid_fractions[i + 1]
				if abs(f2 - f1) > 1e-9:
					T_sol = T1 + (0.001 - f1) / (f2 - f1) * (T2 - T1)
				else:
					T_sol = (T1 + T2) / 2
				break
		
		# 液相线：找到fL从<0.999变为≥0.999的位置
		for i in range(n - 1):
			if liquid_fractions[i] < 0.999 <= liquid_fractions[i + 1]:
				T1, T2 = temperatures[i], temperatures[i + 1]
				f1, f2 = liquid_fractions[i], liquid_fractions[i + 1]
				if abs(f2 - f1) > 1e-9:
					T_liq = T1 + (0.999 - f1) / (f2 - f1) * (T2 - T1)
				else:
					T_liq = (T1 + T2) / 2
				break
		
		# 如果没找到清晰的液相线边界，使用最大fL点
		if T_liq is None:
			max_idx = np.argmax(liquid_fractions)
			if liquid_fractions[max_idx] > 0.5:
				T_liq = temperatures[max_idx]
		
		return T_sol, T_liq
	
	# ============== 批量液相分数计算 ==============

	def _calculate_liquid_fractions_batch (self, temperatures: np.ndarray,
	                                       composition: Dict, components: List[str],
	                                       phases: List[str], model_spec: Optional[Dict],
	                                       pdens: int) -> np.ndarray:
		"""批量计算多个温度点的液相分数（单次equilibrium调用）

		Parameters
		----------
		temperatures : np.ndarray
			温度数组
		composition, components, phases, model_spec, pdens :
			与 _calculate_liquid_fraction_at_T 相同

		Returns
		-------
		np.ndarray
			每个温度点的液相分数，失败的点为0.0
		"""
		try:
			conds = self._build_conditions(composition, temperatures)
			eq = self._call_equilibrium(conds, components, phases, model_spec, pdens)
			if eq is None:
				return np.zeros(len(temperatures))
			phase_info_list = self._extract_phase_info(eq, temperatures)
			fractions = np.array([
				max(0.0, info['liquid_fraction']) for info in phase_info_list
			])
			return fractions
		except Exception:
			return np.zeros(len(temperatures))

	# ============== 优化的二分法（带可选扫描）==============

	def _bisection_with_optional_scan (self,
	                                   composition: Dict,
	                                   T_range: Tuple[float, float],
	                                   components: List[str],
	                                   phases: List[str],
	                                   model_spec: Optional[Dict],
	                                   pdens: int,
	                                   tol: float,
	                                   T_step: float,
	                                   verbose: bool) -> Tuple[
		Optional[float], Optional[float], np.ndarray, np.ndarray]:
		"""
		优化的二分法：快速定位相变温度，同时生成基本的绘图数据

		策略：
		1. 先用低精度粗网格向量化扫描快速定位大致边界
		2. 在边界附近用二分法精确定位
		3. 生成用于绘图的关键点数据（不是完整扫描，但足够绘制主要特征）

		Returns:
		--------
		tuple: (T_solidus, T_liquidus, temperatures, liquid_fractions)
		"""
		T_min, T_max = T_range
		pdens_search = min(500, pdens)  # 搜索阶段用低精度
		calc_params = (composition, components, phases, model_spec, pdens)
		calc_params_search = (composition, components, phases, model_spec, pdens_search)

		# 粗网格扫描（向量化：单次equilibrium调用）
		coarse_step = max(50.0, (T_max - T_min) / 20.0)
		coarse_temps = np.arange(T_min, T_max + coarse_step, coarse_step)

		if verbose:
			self.logger(f"  粗网格扫描: {len(coarse_temps)}个温度点 (步长{coarse_step:.0f}K, pdens={pdens_search})")

		coarse_fractions = self._calculate_liquid_fractions_batch(
			coarse_temps, composition, components, phases, model_spec, pdens_search
		)

		# 二分法精确定位（使用低精度pdens）
		T_sol, _ = self._bisect_boundary(
				coarse_temps, coarse_fractions, calc_params_search,
				target_fL=0.001, search_rising=True, tol=tol
		)

		T_liq, _ = self._bisect_boundary(
				coarse_temps, coarse_fractions, calc_params_search,
				target_fL=0.999, search_rising=True, tol=tol
		)
		
		# 生成绘图数据：合并粗网格点 + 相变边界点 + 少量中间点
		plot_temps = set(coarse_temps)
		
		# 添加精确的边界点
		if T_sol is not None:
			plot_temps.add(round(T_sol, 1))
			# 在固相线附近添加几个点
			for delta in [-5, -2, 2, 5, 10, 20]:
				T_add = T_sol + delta
				if T_min <= T_add <= T_max:
					plot_temps.add(round(T_add, 1))
		
		if T_liq is not None:
			plot_temps.add(round(T_liq, 1))
			# 在液相线附近添加几个点
			for delta in [-20, -10, -5, -2, 2, 5]:
				T_add = T_liq + delta
				if T_min <= T_add <= T_max:
					plot_temps.add(round(T_add, 1))
		
		# 在两相区添加均匀点
		if T_sol is not None and T_liq is not None:
			transition_range = T_liq - T_sol
			if transition_range > 10:
				n_transition_points = min(10, int(transition_range / 5))
				for i in range(1, n_transition_points):
					T_add = T_sol + transition_range * i / n_transition_points
					plot_temps.add(round(T_add, 1))
		
		# 计算所有点的液相分数（使用完整pdens）
		sorted_temps = sorted(plot_temps)
		temperatures = np.array(sorted_temps)
		liquid_fractions = np.zeros(len(sorted_temps))

		# 分离：已有粗网格数据的点 vs 需要新计算的点
		new_temps_indices = []
		for idx_t, T in enumerate(sorted_temps):
			match = np.where(np.abs(coarse_temps - T) < 0.5)[0]
			if len(match) > 0:
				liquid_fractions[idx_t] = coarse_fractions[match[0]]
			else:
				new_temps_indices.append(idx_t)

		# 批量计算新点（使用完整pdens，向量化调用）
		if new_temps_indices:
			new_temps = np.array([sorted_temps[i] for i in new_temps_indices])
			new_fractions = self._calculate_liquid_fractions_batch(
				new_temps, composition, components, phases, model_spec, pdens
			)
			for k, idx_t in enumerate(new_temps_indices):
				liquid_fractions[idx_t] = new_fractions[k]
		
		if verbose:
			self.logger(f"  二分法完成: T_sol={T_sol:.1f}K, T_liq={T_liq:.1f}K" if T_sol and T_liq else "  二分法完成")
		
		return T_sol, T_liq, temperatures, liquid_fractions

	def _pure_bisection (self,
	                     composition: Dict,
	                     T_range: Tuple[float, float],
	                     components: List[str],
	                     phases: List[str],
	                     model_spec: Optional[Dict],
	                     pdens: int,
	                     tol: float,
	                     verbose: bool) -> Tuple[
		Optional[float], Optional[float], np.ndarray, np.ndarray]:
		"""
		优化的纯二分法：智能温度范围扩展 + 鲁棒边界搜索

		改进：
		1. 自适应温度范围扩展（当边界超出范围时）
		2. 更密集的初始网格（20个点而非11个）
		3. 失败回退机制（使用粗网格边界估计）
		4. 详细的诊断日志

		性能：每个成分点约15-25次计算（准确性显著提升）

		Returns:
		--------
		tuple: (T_solidus, T_liquidus, temperatures, liquid_fractions)
		"""
		T_min, T_max = T_range
		calc_params = (composition, components, phases, model_spec, pdens)

		# ========== 第1步：更密集的初始粗网格扫描 ==========
		# 从11个点增加到20个点，减少步长从120K到60K
		coarse_step = max(50.0, (T_max - T_min) / 20.0)
		coarse_temps = np.arange(T_min, T_max + coarse_step, coarse_step)

		if verbose:
			self.logger(f"  优化二分法: {len(coarse_temps)}个初始网格点 (步长{coarse_step:.0f}K)")

		coarse_fractions = []
		calc_count = len(coarse_temps)
		for T in coarse_temps:
			fL = self._calculate_liquid_fraction_at_T(T, *calc_params)
			coarse_fractions.append(max(0.0, fL))

		coarse_fractions = np.array(coarse_fractions)

		# ========== 第2步：智能液相线搜索（带范围扩展） ==========
		T_liq, bisect_count_liq = self._smart_boundary_search(
				coarse_temps, coarse_fractions, calc_params,
				target_fL=0.999, search_rising=True, tol=tol,
				T_range=(T_min, T_max), boundary_name="液相线", verbose=verbose
		)
		calc_count += bisect_count_liq

		# ========== 第3步：智能固相线搜索（带范围扩展） ==========
		T_sol, bisect_count_sol = self._smart_boundary_search(
				coarse_temps, coarse_fractions, calc_params,
				target_fL=0.001, search_rising=True, tol=tol,
				T_range=(T_min, T_max), boundary_name="固相线", verbose=verbose
		)
		calc_count += bisect_count_sol

		# ========== 第4步：保存计算统计 ==========
		temperatures = np.array([calc_count])
		liquid_fractions = np.array([])

		if verbose:
			if T_sol and T_liq:
				self.logger(f"  ✓ 二分法完成: T_sol={T_sol:.1f}K, T_liq={T_liq:.1f}K, 总点数={calc_count}")
			elif T_liq:
				self.logger(f"  ⚠ 仅找到液相线: T_liq={T_liq:.1f}K (固相线搜索失败), 总点数={calc_count}")
			elif T_sol:
				self.logger(f"  ⚠ 仅找到固相线: T_sol={T_sol:.1f}K (液相线搜索失败), 总点数={calc_count}")
			else:
				self.logger(f"  ✗ 液固相线搜索均失败, 总点数={calc_count}")

		return T_sol, T_liq, temperatures, liquid_fractions

	def _smart_boundary_search (self,
	                            coarse_temps: np.ndarray,
	                            coarse_fractions: np.ndarray,
	                            calc_params: tuple,
	                            target_fL: float,
	                            search_rising: bool,
	                            tol: float,
	                            T_range: Tuple[float, float],
	                            boundary_name: str,
	                            verbose: bool) -> Tuple[Optional[float], int]:
		"""
		智能边界搜索：自动扩展温度范围直到找到边界

		策略：
		1. 先在当前粗网格中搜索
		2. 如果没找到跨越点，根据趋势扩展温度范围
		3. 最多扩展3次，每次扩展500K
		4. 如果仍失败，返回粗网格的最佳估计值

		Returns:
		--------
		Tuple[float or None, int]: (边界温度, 额外计算次数)
		"""
		T_min, T_max = T_range
		calc_count = 0

		# 尝试在当前网格搜索
		T_boundary, count = self._bisect_boundary(
				coarse_temps, coarse_fractions, calc_params,
				target_fL, search_rising, tol
		)
		calc_count += count

		if T_boundary is not None:
			return T_boundary, calc_count

		# ========== 搜索失败，分析原因并扩展范围 ==========
		if verbose:
			self.logger(f"    {boundary_name}在初始范围({T_min:.0f}-{T_max:.0f}K)内未找到")

		# 分析液相分数趋势，决定扩展方向
		if target_fL > 0.5:  # 液相线（target_fL=0.999）
			# 如果最高温度处液相分数仍小于目标，需要向上扩展
			if coarse_fractions[-1] < target_fL:
				direction = "up"
				extend_reason = f"最高温{T_max:.0f}K处fL={coarse_fractions[-1]:.3f} < {target_fL}"
			# 如果最低温度处液相分数已大于目标，需要向下扩展
			elif coarse_fractions[0] > target_fL:
				direction = "down"
				extend_reason = f"最低温{T_min:.0f}K处fL={coarse_fractions[0]:.3f} > {target_fL}"
			else:
				# 数据有问题，无法确定方向
				if verbose:
					self.logger(f"    ⚠ {boundary_name}数据异常，无法扩展")
				return self._fallback_estimate(coarse_temps, coarse_fractions, target_fL), calc_count
		else:  # 固相线（target_fL=0.001）
			# 如果最低温度处液相分数仍大于目标，需要向下扩展
			if coarse_fractions[0] > target_fL:
				direction = "down"
				extend_reason = f"最低温{T_min:.0f}K处fL={coarse_fractions[0]:.3f} > {target_fL}"
			# 如果最高温度处液相分数已小于目标，需要向上扩展
			elif coarse_fractions[-1] < target_fL:
				direction = "up"
				extend_reason = f"最高温{T_max:.0f}K处fL={coarse_fractions[-1]:.3f} < {target_fL}"
			else:
				if verbose:
					self.logger(f"    ⚠ {boundary_name}数据异常，无法扩展")
				return self._fallback_estimate(coarse_temps, coarse_fractions, target_fL), calc_count

		# ========== 扩展温度范围并重试 ==========
		if verbose:
			self.logger(f"    正在向{direction}扩展温度范围: {extend_reason}")

		max_extensions = 3
		extend_step = 500.0  # 每次扩展500K

		for i in range(max_extensions):
			# 扩展范围
			if direction == "up":
				T_new_min = T_max + 50.0
				T_new_max = T_max + extend_step
			else:  # down
				T_new_max = T_min - 50.0
				T_new_min = T_min - extend_step

			# 防止温度低于绝对零度
			T_new_min = max(T_new_min, 100.0)
			T_new_max = max(T_new_max, 200.0)

			# 生成新的网格点
			new_temps = np.arange(T_new_min, T_new_max + 50.0, 50.0)
			new_fractions = []
			for T in new_temps:
				fL = self._calculate_liquid_fraction_at_T(T, *calc_params)
				new_fractions.append(max(0.0, fL))
			calc_count += len(new_temps)

			# 合并到原网格
			if direction == "up":
				extended_temps = np.concatenate([coarse_temps, new_temps])
				extended_fractions = np.concatenate([coarse_fractions, new_fractions])
			else:
				extended_temps = np.concatenate([new_temps, coarse_temps])
				extended_fractions = np.concatenate([new_fractions, coarse_fractions])

			# 重新搜索
			T_boundary, count = self._bisect_boundary(
					extended_temps, extended_fractions, calc_params,
					target_fL, search_rising, tol
			)
			calc_count += count

			if T_boundary is not None:
				if verbose:
					self.logger(f"    ✓ 扩展后找到{boundary_name}: {T_boundary:.1f}K (第{i+1}次扩展, 新增{len(new_temps)}点)")
				return T_boundary, calc_count

			# 更新范围继续扩展
			if direction == "up":
				T_max = T_new_max
			else:
				T_min = T_new_min

		# ========== 所有扩展都失败，使用回退估计 ==========
		if verbose:
			self.logger(f"    ⚠ {boundary_name}扩展3次后仍未找到，使用粗网格估计")

		T_fallback = self._fallback_estimate(coarse_temps, coarse_fractions, target_fL)
		return T_fallback, calc_count

	def _fallback_estimate (self, temps: np.ndarray, fractions: np.ndarray, target_fL: float) -> Optional[float]:
		"""
		当二分法失败时，使用粗网格的最接近值作为回退估计

		策略：找到液相分数最接近target_fL的温度点
		"""
		if len(fractions) == 0:
			return None

		# 找到最接近目标的点
		diff = np.abs(fractions - target_fL)
		closest_idx = np.argmin(diff)

		# 如果最接近的误差仍然很大（>0.1），返回None
		if diff[closest_idx] > 0.1:
			return None

		return float(temps[closest_idx])

	def _bisect_boundary (self,
	                      coarse_temps: np.ndarray,
	                      coarse_fractions: np.ndarray,
	                      calc_params: tuple,
	                      target_fL: float,
	                      search_rising: bool,
	                      tol: float) -> Tuple[Optional[float], int]:
		"""
		使用二分法精确定位特定液相分数对应的温度

		Parameters:
		-----------
		coarse_temps : array
			粗网格温度点
		coarse_fractions : array
			对应的液相分数
		calc_params : tuple
			计算参数
		target_fL : float
			目标液相分数（0.001表示固相线，0.999表示液相线）
		search_rising : bool
			True表示搜索fL上升方向的边界
		tol : float
			精度要求

		Returns:
		--------
		Tuple[float or None, int]: (精确的边界温度, 计算次数)
		"""
		# 找到粗网格中跨越目标值的区间
		T_low, T_high = None, None

		for i in range(len(coarse_temps) - 1):
			f1, f2 = coarse_fractions[i], coarse_fractions[i + 1]

			if search_rising:
				if f1 <= target_fL < f2:
					T_low, T_high = coarse_temps[i], coarse_temps[i + 1]
					break
			else:
				if f1 >= target_fL > f2:
					T_low, T_high = coarse_temps[i], coarse_temps[i + 1]
					break

		if T_low is None:
			return None, 0

		# 二分法精确定位
		bisect_count = 0
		for _ in range(20):  # 最多20次迭代
			if T_high - T_low < tol:
				break

			T_mid = (T_low + T_high) / 2.0
			fL_mid = self._calculate_liquid_fraction_at_T(T_mid, *calc_params)
			bisect_count += 1

			if search_rising:
				if fL_mid < target_fL:
					T_low = T_mid
				else:
					T_high = T_mid
			else:
				if fL_mid > target_fL:
					T_low = T_mid
				else:
					T_high = T_mid

		return (T_low + T_high) / 2.0, bisect_count
	
	# ============== 保留的辅助函数 ==============
	
	def _calculate_liquid_fraction_at_T (self, T: float, composition: Dict,
	                                     components: List[str], phases: List[str],
	                                     model_spec: Optional[Dict], pdens: int) -> float:
		"""计算单个温度点的液相分数"""
		try:
			conds = self._build_conditions(composition, np.array([T]))
			eq = self._call_equilibrium(conds, components, phases, model_spec, pdens)
			if eq is None:
				return -1
			phase_info = self._extract_phase_info(eq, np.array([T]))
			if phase_info:
				return phase_info[0]['liquid_fraction']
			return -1
		except Exception:
			return -1
	
	def _call_equilibrium (self, conditions: Dict, components: List[str],
	                       phases: List[str], model_spec: Optional[Dict],
	                       pdens: int) -> Optional[object]:
		"""调用equilibrium进行平衡计算"""
		try:
			# 注意：这里需要导入 pycalphad
			from pycalphad import equilibrium, variables as v
			
			composition = conditions.get('__composition__', {})
			
			if composition:
				filtered_components = []
				for comp_name, comp_value in composition.items():
					if comp_value > 1e-6:
						for comp in components:
							if comp.upper() == comp_name.upper():
								filtered_components.append(comp)
								break
				for comp in components:
					if comp.upper() == 'VA' and comp not in filtered_components:
						filtered_components.append(comp)
						break
			else:
				filtered_components = components
			
			clean_conds = {k: v for k, v in conditions.items() if k != '__composition__'}
			calc_opts = {'pdens': pdens}
			
			if model_spec:
				eq = equilibrium(self.db, filtered_components, phases,
				                 conditions=clean_conds, model=model_spec, calc_opts=calc_opts)
			else:
				eq = equilibrium(self.db, filtered_components, phases,
				                 conditions=clean_conds, calc_opts=calc_opts)
			return eq
		except Exception as e:
			self.logger(f"平衡计算失败: {e}")
			return None
	
	def _build_conditions (self, composition: Dict, temperatures: np.ndarray) -> Dict:
		"""构建平衡计算条件"""
		from pycalphad import variables as v
		
		conds = {v.T: temperatures, v.P: 101325, v.N: 1}
		
		non_zero_comps = {k: val for k, val in composition.items() if val > 1e-6}
		comp_list = list(non_zero_comps.keys())
		n_comps = len(comp_list)
		
		if n_comps > 1:
			total = sum(non_zero_comps.values())
			for comp in comp_list[:-1]:
				conds[v.X(comp)] = non_zero_comps[comp] / total
		
		conds['__composition__'] = composition
		return conds
	
	def _extract_phase_info (self, eq, temperatures):
		"""提取每个温度点的相信息"""
		phases = eq.Phase.values
		amounts = eq.NP.values
		
		phases = np.squeeze(phases)
		amounts = np.squeeze(amounts)
		
		if phases.ndim == 1:
			phases = phases.reshape(1, -1)
			amounts = amounts.reshape(1, -1)
		
		phase_info_list = []
		
		for i in range(len(temperatures)):
			phase_row = phases[i]
			amount_row = amounts[i]
			
			stable_phases = []
			liquid_fraction = 0.0
			
			for ph, amt in zip(phase_row, amount_row):
				ph_str = ph.decode('utf-8') if isinstance(ph, bytes) else str(ph)
				amt_val = float(amt) if hasattr(amt, '__float__') else amt
				
				if amt_val > 0.001:
					phase_name = ph_str.strip()
					if phase_name and phase_name != '':
						stable_phases.append({'name': phase_name, 'fraction': amt_val})
						if 'LIQUID' in phase_name.upper():
							liquid_fraction += amt_val
			
			phase_info_list.append({
				'temperature': temperatures[i],
				'stable_phases': stable_phases,
				'liquid_fraction': liquid_fraction,
				'has_liquid': liquid_fraction > 0.001,
				'is_all_liquid': liquid_fraction >= 0.9999
			})
		
		return phase_info_list



# ============================================================================
# 溶解度计算器
# ============================================================================

class SolubilityCalculator(_ContributionCoefficientMixin):
	"""溶解度计算器 - 返回统一的CalculationResult"""
	
	STANDARD_SOLUTION_PHASES = {
		'LIQUID', 'FCC_A1', 'BCC_A2', 'HCP_A3', 'HCP_ZN',
		'DIAMOND_A4', 'CBCC_A12', 'CUB_A13', 'BCT_A5', 'DHCP'
	}
	EXCLUDED_SUFFIXES = ['_INT_', '_ORD', '_L12', '_L10', '_D0', '_THETA']
	
	def __init__ (self, database: Database, model_spec=None, pdens: int = 2000, logger=None):
		self.db = database
		self.model_spec = model_spec
		self.pdens = pdens
		self.logger = logger if logger else print
		self._model_info_printed = False
	
	def calculate (self, base_composition: Dict[str, float],
	               solute: str, temperature: float,
	               x_max: float = 0.95,
	               tolerance: float = 1e-6,
	               max_iterations: int = 30,
	               phases: Optional[List[str]] = None,
	               verbose: bool = False,
	               bracket_hint: Optional[Tuple[float, float]] = None) -> CalculationResult:
		"""
		计算溶质在基础合金中的溶解度

		Returns:
		--------
		CalculationResult with mode='solubility'
		"""
		try:
			base_composition = {k.upper(): val for k, val in base_composition.items()}
			solute = solute.upper()
			if verbose and not self._model_info_printed:
				base_str = '+'.join([f'{k.upper()}{v:.2f}' for k, v in base_composition.items()])
				self.logger(f"基础组分: {base_str}, 溶质: {solute.upper()}")
				self._model_info_printed = True
			
			# 构建组分列表
			components = set(k.upper() for k in base_composition.keys())
			components.add(solute.upper())
			components.add('VA')
			components = sorted(list(components))

			# ★ 打印贡献系数（如果使用非默认模型）
			if verbose and self.model_spec:
				components_no_va = [c for c in components if c != 'VA']
				self._print_contribution_coefficients(components_no_va, temperature, self.model_spec)

			# 获取相列表（仅在调用方未指定时自动检测）
			if phases is None:
				phases = self._get_all_phases_for_system(components)
			
			if not phases:
				return CalculationResult(
						success=False,
						mode=CalculationMode.SOLUBILITY.value,
						message="没有找到兼容的相",
						solubility_data={'temperature_K': temperature, 'error': '没有找到兼容的相'}
				)
			
			# 识别基体相
			matrix_phase, error = self._identify_matrix_phase(base_composition, temperature, verbose=False)
			
			if error:
				return CalculationResult(
						success=False,
						mode=CalculationMode.SOLUBILITY.value,
						message=error,
						solubility_data={'temperature_K': temperature, 'error': error}
				)
			
			if matrix_phase not in phases:
				phases.append(matrix_phase)
			
			# 检查微量溶质
			trace = 1e-8
			is_single, phase_info = self._check_single_phase(
					base_composition, solute.upper(), trace, temperature,
					matrix_phase, components, phases, verbose=False
			)
			
			if not is_single:
				return CalculationResult(
						success=True,
						mode=CalculationMode.SOLUBILITY.value,
						message=f"溶解度: 0.0 at% @ {temperature}K (基体不稳定)",
						solubility_data={
							'temperature_K': temperature,
							'solubility': 0.0,
							'matrix_phase': matrix_phase,
							'composition': base_composition
						}
				)
			
			# 二分搜索（使用括号提示缩窄搜索范围）
			if bracket_hint is not None:
				hint_low, hint_high = bracket_hint
				margin = max((hint_high - hint_low) * 0.5, 0.01)
				x_low = max(0.0, hint_low - margin)
				x_high = min(x_max, hint_high + margin)
				# 验证括号有效性
				is_low_ok, _ = self._check_single_phase(
						base_composition, solute.upper(), x_low, temperature,
						matrix_phase, components, phases, verbose=False
				)
				is_high_ok, _ = self._check_single_phase(
						base_composition, solute.upper(), x_high, temperature,
						matrix_phase, components, phases, verbose=False
				)
				if not (is_low_ok and not is_high_ok):
					# 括号无效，回退到完整范围
					x_low, x_high = 0.0, x_max
			else:
				x_low, x_high = 0.0, x_max
			
			for i in range(max_iterations):
				x_mid = (x_low + x_high) / 2
				is_single, phase_info = self._check_single_phase(
						base_composition, solute.upper(), x_mid, temperature,
						matrix_phase, components, phases, verbose=False
				)
				if is_single:
					x_low = x_mid
				else:
					x_high = x_mid
				if x_high - x_low < tolerance:
					break
			
			if verbose:
				self.logger(f"✓ 溶解度: {x_low * 100:.4f} at%, 稳定相: {matrix_phase}")
			
			return CalculationResult(
					success=True,
					mode=CalculationMode.SOLUBILITY.value,
					message=f"溶解度: {x_low * 100:.4f} at% @ {temperature}K",
					solubility_data={
						'temperature_K': temperature,
						'solubility': x_low,
						'solubility_percent': x_low * 100,
						'matrix_phase': matrix_phase,
						'composition': base_composition,
						'solute': solute.upper()
					}
			)
		
		except Exception as e:
			if verbose:
				self.logger(traceback.format_exc())
			return CalculationResult(
					success=False,
					mode=CalculationMode.SOLUBILITY.value,
					message=f"计算失败: {str(e)}",
					solubility_data={'temperature_K': temperature, 'error': str(e)}
			)
	
	def calculate_solvus (self, base_composition: Dict[str, float],
	                      solute: str, T_range: Tuple[float, float],
	                      T_points: int = 20, phases: Optional[List[str]] = None, verbose: bool = False) -> CalculationResult:
		"""
		计算溶解度曲线

		Returns:
		--------
		CalculationResult with mode='solvus'
		"""
		try:
			temperatures = np.linspace(T_range[0], T_range[1], T_points)
			solubilities = []
			result_phases = []
			errors = []

			prev_bracket = None
			for T in temperatures:
				result = self.calculate(base_composition, solute, T, phases=phases,
				                        verbose=False, bracket_hint=prev_bracket)

				if result.success and result.solubility_data:
					sol = result.solubility_data.get('solubility', np.nan)
					solubilities.append(sol)
					result_phases.append(result.solubility_data.get('matrix_phase'))
					errors.append(None)
					# 用当前结果缩窄下一个温度点的搜索范围
					if not np.isnan(sol) and sol > 0:
						prev_bracket = (0.0, sol * 1.5)
					else:
						prev_bracket = None
				else:
					solubilities.append(np.nan)
					result_phases.append(None)
					errors.append(result.message)
					prev_bracket = None  # 失败时重置

			return CalculationResult(
					success=True,
					mode=CalculationMode.SOLVUS.value,
					message=f"溶解度曲线计算完成: {T_points}个温度点",
					x_axis=temperatures,
					y_axis=np.array(solubilities),
					solubility_data={
						'temperatures': temperatures,
						'solubilities': np.array(solubilities),
						'phases': result_phases,
						'errors': errors,
						'solute': solute.upper(),
						'base_composition': base_composition
					}
			)
		
		except Exception as e:
			return CalculationResult(
					success=False,
					mode=CalculationMode.SOLVUS.value,
					message=f"溶解度曲线计算失败: {str(e)}"
			)
	
	def _is_standard_solution_phase (self, phase_name: str) -> bool:
		"""判断是否是标准溶液相"""
		phase_upper = phase_name.upper()
		for suffix in self.EXCLUDED_SUFFIXES:
			if suffix in phase_upper:
				return False
		if phase_upper in self.STANDARD_SOLUTION_PHASES:
			return True
		if phase_upper.startswith('LIQUID'):
			return True
		return False
	
	def _get_phase_species (self, phase_name: str) -> set:
		"""获取相的物种集合"""
		phase = self.db.phases[phase_name]
		species_set = set()
		for sublattice in phase.constituents:
			for species in sublattice:
				name = species.name if hasattr(species, 'name') else str(species)
				if name.upper() not in ['VA', '*', 'VACANCY']:
					species_set.add(name.upper())
		return species_set
	
	def _get_all_phases_for_system (self, components: List[str]) -> List[str]:
		"""获取指定组分系统的所有相关相"""
		comp_set = set(c.upper() for c in components if c.upper() != 'VA')
		solution_phases = []
		compound_phases = []
		
		for phase_name in self.db.phases:
			try:
				phase_species = self._get_phase_species(phase_name)
				if not phase_species:
					continue
				if phase_species.issubset(comp_set):
					compound_phases.append(phase_name)
				elif comp_set.issubset(phase_species):
					if self._is_standard_solution_phase(phase_name):
						solution_phases.append(phase_name)
			except Exception:
				continue
		
		return list(set(solution_phases + compound_phases))
	
	def _get_solution_phases_for_components (self, components: List[str]) -> List[str]:
		"""获取能容纳指定组分的溶液相"""
		comp_set = set(c.upper() for c in components if c.upper() != 'VA')
		solution_phases = []
		
		for phase_name in self.db.phases:
			try:
				if not self._is_standard_solution_phase(phase_name):
					continue
				phase_species = self._get_phase_species(phase_name)
				if comp_set.issubset(phase_species):
					solution_phases.append(phase_name)
			except Exception:
				continue

		return solution_phases
	
	def _call_equilibrium (self, components: List[str], phases: List[str],
	                       conditions: Dict) -> Optional[object]:
		"""调用平衡计算"""
		try:
			if self.model_spec:
				eq = equilibrium(self.db, components, phases, conditions,
				                 model=self.model_spec, calc_opts={'pdens': self.pdens})
			else:
				eq = equilibrium(self.db, components, phases, conditions,
				                 calc_opts={'pdens': self.pdens})
			return eq
		except Exception as e:
			self.logger(f"平衡计算失败: {e}")
			return None
	
	def _extract_stable_phases (self, eq) -> List[Tuple[str, float]]:
		"""提取稳定相"""
		if eq is None:
			return []
		
		phase_names = eq.Phase.values.flatten()
		phase_amounts = eq.NP.values.flatten()
		
		valid_phases = []
		for ph, amt in zip(phase_names, phase_amounts):
			if isinstance(ph, bytes):
				ph = ph.decode('utf-8')
			ph_str = str(ph).strip()
			if ph_str and ph_str != '' and not np.isnan(amt) and amt > 0.001:
				valid_phases.append((ph_str, float(amt)))
		
		valid_phases.sort(key=lambda x: x[1], reverse=True)
		return valid_phases
	
	def _identify_matrix_phase (self, base_composition: Dict[str, float],
	                            temperature: float, verbose: bool = False) -> Tuple[Optional[str], Optional[str]]:
		"""使用平衡计算确定基体相"""
		total = sum(base_composition.values())
		base_norm = {k.upper(): val / total for k, val in base_composition.items()}
		
		base_elements = sorted(base_norm.keys())
		base_components = base_elements + ['VA']
		
		solution_phases = self._get_solution_phases_for_components(base_components)
		
		if not solution_phases:
			return None, f"没有找到适用的溶液相"
		
		conds = {v.T: temperature, v.P: 101325, v.N: 1}
		if len(base_elements) > 1:
			for elem in base_elements[:-1]:
				conds[v.X(elem)] = base_norm[elem]
		
		eq = self._call_equilibrium(base_components, solution_phases, conds)
		
		if eq is None:
			return None, "基础合金平衡计算失败"
		
		valid_phases = self._extract_stable_phases(eq)
		
		if not valid_phases:
			return None, "无法确定稳定相"
		
		return valid_phases[0][0], None

	def _check_single_phase (self, base_composition: Dict[str, float],
	                         solute: str, x_solute: float, temperature: float,
	                         matrix_phase: str, components: List[str],
	                         phases: List[str], verbose: bool = False) -> Tuple[bool, str]:
		"""检查给定溶质含量下是否仍为单一基体相"""
		total = sum(base_composition.values())
		base_norm = {k.upper(): val / total for k, val in base_composition.items()}

		remaining = 1.0 - x_solute

		conds = {v.T: temperature, v.P: 101325, v.N: 1}
		conds[v.X(solute)] = x_solute

		base_elements = sorted(base_norm.keys())
		for elem in base_elements[:-1]:
			conds[v.X(elem)] = remaining * base_norm[elem]

		eq = self._call_equilibrium(components, phases, conds)

		if eq is None:
			return False, "计算失败"

		valid_phases = self._extract_stable_phases(eq)

		if not valid_phases:
			return False, "Unknown"

		if len(valid_phases) == 1 and valid_phases[0][0] == matrix_phase:
			return True, matrix_phase
		else:
			return False, "+".join([p for p, _ in valid_phases])


# ============================================================================
# 伪二元相图计算器
# 主类：PseudoBinaryCalculator
# ============================================================================

class PseudoBinaryCalculator(_ContributionCoefficientMixin):
    """
    伪二元相图计算器 - 返回统一的 CalculationResult

    功能：
    - 真二元系统: mode='binplot', 返回 matplotlib Figure
    - 伪二元系统: mode='phase_map', 返回相区图矩阵

    改进：
    - 自适应网格细化 (AMR) 优化计算效率
    - 智能相区标注（自动避让、边界细化）
    - 统一的相提取逻辑
    """

    def __init__(self, database: Database, logger: Optional[Callable] = None):
        """
        初始化计算器

        Parameters
        ----------
        database : Database
            PyCalphad 数据库对象
        logger : Callable, optional
            日志函数，默认为空操作
        """
        self.db = database
        self.logger = logger or (lambda msg: None)

    # ========================================================================
    # 主计算方法
    # ========================================================================

    def calculate(self,
                  base_composition: Dict,
                  varying_comp: str,
                  x_range: Tuple[float, float],
                  x_points: int,
                  T_range: Tuple[float, float],
                  components: List[str],
                  phases: Optional[List[str]] = None,
                  model_spec: Optional[Dict] = None,
                  pdens: int = 2000,
                  T_points: int = 51,
                  verbose: bool = False) -> CalculationResult:
        """
        计算伪二元相图 - 生成 T-X 相区图

        Parameters
        ----------
        base_composition : dict
            基础成分字典，例如 {'AL': 0.7, 'MG': 0.3}
        varying_comp : str
            变化的组分名称，例如 'CU'
        x_range : tuple
            变化组分的摩尔分数范围 (x_min, x_max)
        x_points : int
            成分轴上的网格点数
        T_range : tuple
            温度范围 (T_min, T_max)，单位 K
        components : List[str]
            组分列表（包含 VA）
        phases : List[str], optional
            相列表。如果为 None，使用数据库中所有相
        model_spec : Dict, optional
            模型规范（如 UEM 模型字典）
        pdens : int
            点密度参数
        T_points : int
            温度轴上的网格点数（默认 51）
        verbose : bool
            详细输出

        Returns
        -------
        CalculationResult
            mode='phase_map' 或 'binplot'
        """
        try:
            # 参数验证
            x_min, x_max = x_range
            T_min, T_max = T_range

            if not (0 <= x_min < x_max <= 1):
                return CalculationResult(
                    success=False,
                    mode=CalculationMode.PHASE_MAP.value,
                    message=f"x_range 必须在 [0, 1] 范围内: {x_range}"
                )

            # 标准化组分名称为全大写
            standardized_base = {k.upper(): val for k, val in base_composition.items()}
            standardized_varying = varying_comp.upper()

            # 过滤组分列表：只保留实际参与计算的组分
            active_components = set(standardized_base.keys())
            active_components.add(standardized_varying)
            active_components.add('VA')
            filtered_components = sorted(list(active_components))

            if verbose:
                self.logger(f"计算伪二元相图: 基础={list(standardized_base.keys())}, 变化={standardized_varying}")
                self.logger(f"活性组分: {filtered_components}")

            # 相过滤
            if phases is None:
                phases = list(self.db.phases.keys())

            # 使用严格的相过滤函数
            components_no_va = [c for c in filtered_components if c != 'VA']
            filtered_phases = self._filter_phases_strict(components_no_va)

            # 如果用户提供了 phases 列表，取交集
            if phases != list(self.db.phases.keys()):
                filtered_phases = [p for p in filtered_phases if p in phases]

            if verbose:
                self.logger(f"相过滤: {len(filtered_phases)}/{len(phases)} 个相")

            # 打印贡献系数（如果使用非默认模型）
            if verbose and model_spec:
                mid_temp = (T_min + T_max) / 2
                self._print_contribution_coefficients(components_no_va, mid_temp, model_spec)

            # ================================================================
            # 智能路由：真二元 vs 伪二元
            # ================================================================
            is_true_binary = self._is_true_binary(standardized_base, standardized_varying)

            if is_true_binary:
                # 真二元系统：使用 PyCalphad 内置 binplot
                if verbose:
                    self.logger(f"✓ 检测到真二元系统，使用 PyCalphad binplot 绘制相边界")

                return self._calculate_true_binary_with_binplot(
                    filtered_components, x_range, x_points,
                    T_range, T_points, filtered_phases,
                    model_spec, pdens, verbose
                )
            else:
                # 伪二元系统：使用多进程并行计算相组合图
                if verbose:
                    self.logger(f"✓ 检测到伪二元系统，使用多进程并行计算")

                total_points = T_points * x_points
                if verbose:
                    self.logger(f"网格: {x_points} 成分点 × {T_points} 温度点 = {total_points} 计算点")

                # 使用 AMR 优化的并行计算（返回整数编码的phase_map + 映射表）
                phase_map, phase_legend = self._calculate_via_parallel(
                    standardized_base, standardized_varying,
                    x_range, x_points, T_range, T_points,
                    filtered_components, filtered_phases,
                    model_spec, pdens, verbose
                )

                x_values = np.linspace(x_min, x_max, x_points)
                T_values = np.linspace(T_min, T_max, T_points)

                if verbose:
                    self.logger(f"✓ 相图计算完成")

                return CalculationResult(
                    success=True,
                    mode=CalculationMode.PHASE_MAP.value,
                    message=f"伪二元相图计算完成: {x_points}×{T_points} 网格点",
                    x_axis=x_values,
                    t_grid=T_values,
                    phase_map=phase_map,
                    phase_legend=phase_legend,
                    x_range=x_range,
                    t_range=T_range,
                    components=[c for c in filtered_components if c != 'VA']
                )

        except Exception as e:
            if verbose:
                self.logger(traceback.format_exc())
            return CalculationResult(
                success=False,
                mode=CalculationMode.PHASE_MAP.value,
                message=f"相图计算失败: {str(e)}"
            )

    # ========================================================================
    # 核心方法：相提取（静态方法）
    # ========================================================================

    @staticmethod
    def _extract_stable_phases(eq_result, threshold: float = 1e-3) -> str:
        """
        从平衡计算结果中提取稳定相组合字符串

        Parameters
        ----------
        eq_result : xarray.Dataset
            PyCalphad 平衡计算结果
        threshold : float
            相分数阈值，低于此值的相被忽略（默认 0.001）

        Returns
        -------
        str
            稳定相组合字符串，格式如 "FCC_A1 + LIQUID"
            - 如果无稳定相返回 "No Phase"
            - 如果计算出错返回 "Error"

        Notes
        -----
        ★ 关键改进：相名按字母排序，确保 "FCC + LIQUID" 和 "LIQUID + FCC"
          被识别为同一相区
        """
        try:
            phase_values = eq_result.Phase.values
            np_values = eq_result.NP.values

            # 压缩维度
            phase_values = np.squeeze(phase_values)
            np_values = np.squeeze(np_values)

            # 处理标量情况
            if phase_values.ndim == 0:
                if phase_values.size == 0:
                    return "No Phase"
                phase_values = np.array([phase_values.item()])
                np_values = np.array([np_values.item()])

            # 收集稳定相
            stable_phases = []
            for phase, amount in zip(phase_values, np_values):
                # 解码 bytes
                if isinstance(phase, bytes):
                    phase = phase.decode('utf-8')

                phase_str = str(phase).strip()

                # 转换 amount 为 float
                try:
                    amt_float = float(amount)
                except (TypeError, ValueError):
                    continue

                # 过滤无效相和微量相
                if phase_str and phase_str != '' and amt_float > threshold:
                    clean_name = phase_str.strip()
                    if clean_name:
                        stable_phases.append(clean_name)

            if not stable_phases:
                return "No Phase"

            # ★ 关键：排序 + 去重后连接，确保相名一致性
            return " + ".join(sorted(set(stable_phases)))

        except Exception:
            return "Error"

    # ========================================================================
    # 真二元相图：binplot + 相区标注
    # ========================================================================

    def _calculate_true_binary_with_binplot(self,
                                            binary_comps: List[str],
                                            x_range: Tuple[float, float],
                                            x_points: int,
                                            T_range: Tuple[float, float],
                                            T_points: int,
                                            phases: List[str],
                                            model_spec: Optional[Dict],
                                            pdens: int,
                                            verbose: bool) -> CalculationResult:
        """
        使用 PyCalphad 的 binplot 绘制真二元相图
        """
        try:
            from pycalphad.mapping.compat_api import binplot

            x_min, x_max = x_range
            T_min, T_max = T_range

            # 确保只有两个非 VA 组分
            binary_comps_filtered = [c for c in binary_comps if c != 'VA']
            if len(binary_comps_filtered) != 2:
                return CalculationResult(
                    success=False,
                    mode=CalculationMode.BINPLOT.value,
                    message=f"真二元系统需要恰好 2 个非 VA 组分，得到: {binary_comps_filtered}"
                )

            # 添加 VA 组分
            comps_with_va = sorted(binary_comps_filtered) + ['VA']

            if verbose:
                self.logger(f"使用 PyCalphad binplot: {binary_comps_filtered[0]}-{binary_comps_filtered[1]}")
                self.logger(f"温度范围: {T_min:.0f}-{T_max:.0f}K, 成分范围: {x_min:.3f}-{x_max:.3f}")

            # 设置计算条件
            conditions = {
                v.T: (T_min, T_max, T_points),
                v.X(comps_with_va[0]): (x_min, x_max, x_points),
                v.P: 101325,
                v.N: 1
            }

            calc_opts = {'pdens': pdens}

            # 创建图形
            fig, ax = plt.subplots(figsize=(10, 8))

            plot_kwargs = {
                'ax': ax,
                'label_nodes': False
            }

            # 调用 binplot
            if model_spec:
                ax_result, strategy = binplot(
                    self.db, comps_with_va, phases, conditions,
                    plot_kwargs=plot_kwargs,
                    calc_opts=calc_opts,
                    model=model_spec,
                    return_strategy=True
                )
            else:
                ax_result, strategy = binplot(
                    self.db, comps_with_va, phases, conditions,
                    plot_kwargs=plot_kwargs,
                    calc_opts=calc_opts,
                    return_strategy=True
                )

            # 添加相区标注
            self._add_phase_labels_to_binplot(
                ax_result, strategy, binary_comps_filtered,
                x_range, T_range, phases, model_spec, pdens, verbose
            )

            if verbose:
                self.logger(f"✓ binplot 相图计算完成")

            return CalculationResult(
                success=True,
                mode=CalculationMode.BINPLOT.value,
                message=f"真二元相图 {binary_comps_filtered[0]}-{binary_comps_filtered[1]} 计算完成",
                figure=fig,
                x_range=x_range,
                t_range=T_range,
                components=binary_comps_filtered
            )

        except Exception as e:
            if verbose:
                self.logger(f"✗ binplot 计算失败: {e}")
                self.logger(traceback.format_exc())
            return CalculationResult(
                success=False,
                mode=CalculationMode.BINPLOT.value,
                message=f"binplot 计算失败: {str(e)}"
            )

    def _add_phase_labels_to_binplot(self,
                                      ax,
                                      strategy,
                                      components: List[str],
                                      x_range: Tuple[float, float],
                                      T_range: Tuple[float, float],
                                      phases: List[str],
                                      model_spec,
                                      pdens: int,
                                      verbose: bool,
                                      # ===== 可调参数 =====
                                      n_samples_x: int = 35,
                                      n_samples_T: int = 40,
                                      min_points_for_label: int = 3,
                                      use_adaptive_refinement: bool = True):
        """
        在 binplot 相图上添加相区名称标注

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Matplotlib 坐标轴对象
        strategy : BinaryStrategy
            PyCalphad 二元策略对象
        components : List[str]
            组分列表（不含 VA）
        x_range : Tuple[float, float]
            成分范围 (x_min, x_max)
        T_range : Tuple[float, float]
            温度范围 (T_min, T_max) [K]
        phases : List[str]
            相列表
        model_spec : Optional
            模型规格
        pdens : int
            点密度参数
        verbose : bool
            是否输出详细日志
        n_samples_x : int
            X 方向采样点数（默认 35）
        n_samples_T : int
            T 方向采样点数（默认 40）
        min_points_for_label : int
            标注相区所需的最少采样点数（默认 3）
        use_adaptive_refinement : bool
            是否使用自适应边界细化（默认 True）
        """
        try:
            start_time = time.time()

            # =================================================================
            # 1. 准备工作
            # =================================================================

            # 使用严格的相过滤
            filtered_phases = self._filter_phases_strict(components)

            # 如果用户指定了 phases，取交集
            if phases != list(self.db.phases.keys()):
                filtered_phases = [p for p in filtered_phases if p in phases]

            if not filtered_phases:
                if verbose:
                    self.logger("  [标注] 警告: 没有兼容的相，跳过标注")
                return

            if verbose:
                self.logger(f"  [标注] 相过滤: {len(phases)} → {len(filtered_phases)} 个兼容相")

            # 构建组分列表（含 VA）
            comps = sorted(components) + ['VA']

            x_min, x_max = x_range
            T_min, T_max = T_range

            # =================================================================
            # 2. 粗网格采样
            # =================================================================

            x_samples = np.linspace(x_min, x_max, n_samples_x)
            T_samples = np.linspace(T_min, T_max, n_samples_T)

            # 存储采样结果
            sample_results: Dict[Tuple[int, int], str] = {}
            phase_regions: Dict[str, List[Tuple[float, float]]] = defaultdict(list)

            total_samples = n_samples_x * n_samples_T
            success_count = 0

            if verbose:
                self.logger(f"  [标注] 粗网格采样: {n_samples_x}×{n_samples_T} = {total_samples} 点")

            # 降低 pdens 加速标注采样
            label_pdens = min(500, pdens)
            calc_opts = {'pdens': label_pdens}

            # 向量化采样：对每个x值，用T_samples数组做单次equilibrium调用
            for i, x_val in enumerate(x_samples):
                try:
                    conds = {
                        v.T: T_samples,  # 传入完整T数组（向量化）
                        v.P: 101325,
                        v.N: 1,
                        v.X(comps[0]): x_val
                    }

                    if model_spec:
                        eq = equilibrium(self.db, comps, filtered_phases, conds,
                                       model=model_spec, calc_opts=calc_opts)
                    else:
                        eq = equilibrium(self.db, comps, filtered_phases, conds,
                                       calc_opts=calc_opts)

                    # 从向量化结果中逐T提取相信息
                    phase_values = np.squeeze(eq.Phase.values)
                    np_values = np.squeeze(eq.NP.values)

                    if phase_values.ndim == 1:
                        phase_values = phase_values.reshape(1, -1)
                        np_values = np_values.reshape(1, -1)

                    for j in range(len(T_samples)):
                        try:
                            stable = []
                            for ph, amt in zip(phase_values[j], np_values[j]):
                                if isinstance(ph, bytes):
                                    ph = ph.decode('utf-8')
                                ph_str = str(ph).strip()
                                try:
                                    amt_f = float(amt)
                                except (TypeError, ValueError):
                                    continue
                                if ph_str and ph_str != '' and amt_f > 1e-3:
                                    stable.append(ph_str)
                            if stable:
                                phase_str = " + ".join(sorted(set(stable)))
                                sample_results[(i, j)] = phase_str
                                phase_regions[phase_str].append((x_val, T_samples[j]))
                                success_count += 1
                        except Exception:
                            continue

                except Exception:
                    # 回退：逐点计算此x值
                    for j, T_val in enumerate(T_samples):
                        try:
                            conds_single = {
                                v.T: T_val, v.P: 101325, v.N: 1,
                                v.X(comps[0]): x_val
                            }
                            if model_spec:
                                eq = equilibrium(self.db, comps, filtered_phases, conds_single,
                                               model=model_spec, calc_opts=calc_opts)
                            else:
                                eq = equilibrium(self.db, comps, filtered_phases, conds_single,
                                               calc_opts=calc_opts)
                            phase_str = self._extract_stable_phases(eq)
                            if phase_str and phase_str not in ("Error", "No Phase"):
                                sample_results[(i, j)] = phase_str
                                phase_regions[phase_str].append((x_val, T_val))
                                success_count += 1
                        except Exception:
                            continue

            coarse_time = time.time() - start_time

            if verbose:
                self.logger(f"  [标注] 粗网格完成: {success_count}/{total_samples} 成功, "
                          f"耗时 {coarse_time:.1f}s (向量化)")
                self.logger(f"  [标注] 识别到 {len(phase_regions)} 个相区")

            # =================================================================
            # 3. 自适应边界细化（可选）
            # =================================================================

            if use_adaptive_refinement and len(phase_regions) > 1:
                refine_start = time.time()

                # 检测相边界
                boundary_regions = self._detect_phase_boundaries(
                    sample_results, x_samples, T_samples, n_samples_x, n_samples_T
                )

                # 在边界区域细化采样
                refined_count = 0
                for (x_lo, x_hi, T_lo, T_hi) in boundary_regions:
                    x_fine = np.linspace(x_lo, x_hi, 5)
                    T_fine = np.linspace(T_lo, T_hi, 5)

                    for x_val in x_fine:
                        for T_val in T_fine:
                            try:
                                conds = {
                                    v.T: T_val,
                                    v.P: 101325,
                                    v.N: 1,
                                    v.X(comps[0]): x_val
                                }

                                if model_spec:
                                    eq = equilibrium(self.db, comps, filtered_phases, conds,
                                                   model=model_spec, calc_opts=calc_opts)
                                else:
                                    eq = equilibrium(self.db, comps, filtered_phases, conds,
                                                   calc_opts=calc_opts)

                                phase_str = self._extract_stable_phases(eq)

                                if phase_str and phase_str not in ("Error", "No Phase"):
                                    phase_regions[phase_str].append((x_val, T_val))
                                    refined_count += 1

                            except Exception:
                                continue

                refine_time = time.time() - refine_start

                if verbose and refined_count > 0:
                    self.logger(f"  [标注] 边界细化: +{refined_count} 点, 耗时 {refine_time:.1f}s")

            # =================================================================
            # 4. 计算相区中心
            # =================================================================

            # 按面积排序（大相区先标注）
            sorted_regions = sorted(
                phase_regions.items(),
                key=lambda x: -len(x[1])
            )

            # 已占用位置列表
            occupied_positions: List[Tuple[float, float, float, float]] = []

            # 估算标签尺寸
            label_width = (x_max - x_min) * 0.12
            label_height = (T_max - T_min) * 0.04

            labeled_count = 0

            for phase_str, points in sorted_regions:
                n_points = len(points)

                # 过滤点数不足的相区
                if n_points < min_points_for_label:
                    if verbose:
                        self.logger(f"  [标注] 跳过 '{phase_str}': 仅 {n_points} 点")
                    continue

                points_array = np.array(points)

                # 计算中心点
                center_x = np.mean(points_array[:, 0])
                center_y = np.mean(points_array[:, 1])

                # 对窄带状区域使用中位数
                x_spread = np.max(points_array[:, 0]) - np.min(points_array[:, 0])
                T_spread = np.max(points_array[:, 1]) - np.min(points_array[:, 1])

                if x_spread < (x_max - x_min) * 0.1:
                    center_x = np.median(points_array[:, 0])
                if T_spread < (T_max - T_min) * 0.1:
                    center_y = np.median(points_array[:, 1])

                # 寻找不重叠的位置
                best_x, best_y = self._find_non_overlapping_position(
                    center_x, center_y,
                    label_width, label_height,
                    occupied_positions,
                    x_range, T_range,
                    points
                )

                # 根据相区大小调整字体
                if n_points < 10:
                    fontsize = 8
                elif n_points < 30:
                    fontsize = 9
                else:
                    fontsize = 10

                # 长名称换行
                display_str = phase_str
                if len(display_str) > 15 and '+' in display_str:
                    display_str = display_str.replace(' + ', '\n+ ')

                # 添加标注
                ax.text(
                    best_x, best_y, display_str,
                    fontsize=fontsize,
                    fontweight='bold',
                    ha='center',
                    va='center',
                    bbox=dict(
                        boxstyle='round,pad=0.3',
                        facecolor='white',
                        alpha=0.85,
                        edgecolor='gray',
                        linewidth=0.8
                    ),
                    zorder=100
                )

                # 记录已占用位置
                occupied_positions.append((
                    best_x - label_width / 2,
                    best_y - label_height / 2,
                    label_width,
                    label_height
                ))

                labeled_count += 1

                if verbose:
                    self.logger(f"  [标注] ✓ '{phase_str}' @ ({best_x:.3f}, {best_y:.0f}K), "
                              f"{n_points} 点")

            # =================================================================
            # 5. 完成
            # =================================================================

            total_time = time.time() - start_time

            if verbose:
                self.logger(f"  [标注] 完成: {labeled_count} 个标注, 总耗时 {total_time:.1f}s")

        except Exception as e:
            if verbose:
                self.logger(f"  [标注] 错误: {e}")
                self.logger(f"  {traceback.format_exc()}")

    def _detect_phase_boundaries(self,
                                  sample_results: Dict[Tuple[int, int], str],
                                  x_samples: np.ndarray,
                                  T_samples: np.ndarray,
                                  n_x: int,
                                  n_T: int) -> List[Tuple[float, float, float, float]]:
        """
        检测相边界区域

        通过比较相邻采样点的相是否相同来识别边界

        Returns
        -------
        List[Tuple[float, float, float, float]]
            边界区域列表 [(x_lo, x_hi, T_lo, T_hi), ...]
        """
        boundary_regions = []

        for i in range(n_x - 1):
            for j in range(n_T - 1):
                # 获取 2×2 区块的四个角
                corners = [
                    sample_results.get((i, j)),
                    sample_results.get((i + 1, j)),
                    sample_results.get((i, j + 1)),
                    sample_results.get((i + 1, j + 1))
                ]

                # 过滤无效值
                valid_corners = [c for c in corners if c and c not in ("Error", "No Phase")]

                # 如果有效角点的相不完全相同 → 边界区域
                if len(set(valid_corners)) > 1:
                    x_lo, x_hi = x_samples[i], x_samples[i + 1]
                    T_lo, T_hi = T_samples[j], T_samples[j + 1]
                    boundary_regions.append((x_lo, x_hi, T_lo, T_hi))

        return boundary_regions

    def _find_non_overlapping_position(self,
                                         init_x: float,
                                         init_y: float,
                                         label_w: float,
                                         label_h: float,
                                         occupied: List[Tuple[float, float, float, float]],
                                         x_range: Tuple[float, float],
                                         T_range: Tuple[float, float],
                                         region_points: List[Tuple[float, float]]) -> Tuple[float, float]:
        """
        寻找不与已有标签重叠的位置

        策略：
        1. 首先尝试原始中心位置
        2. 如果重叠，尝试在相区内的其他位置
        3. 如果仍然重叠，使用轻微偏移

        Returns
        -------
        Tuple[float, float]
            最佳标注位置 (x, y)
        """
        x_min, x_max = x_range
        T_min, T_max = T_range

        def is_overlapping(x: float, y: float) -> bool:
            for ox, oy, ow, oh in occupied:
                if (abs(x - (ox + ow / 2)) < (label_w + ow) / 2 and
                    abs(y - (oy + oh / 2)) < (label_h + oh) / 2):
                    return True
            return False

        def is_in_bounds(x: float, y: float) -> bool:
            margin_x = label_w / 2
            margin_y = label_h / 2
            return (x_min + margin_x <= x <= x_max - margin_x and
                    T_min + margin_y <= y <= T_max - margin_y)

        # 尝试 1：原始位置
        if not is_overlapping(init_x, init_y):
            return init_x, init_y

        # 尝试 2：相区内的四分位点
        if len(region_points) > 5:
            points_array = np.array(region_points)
            for q in [0.25, 0.5, 0.75]:
                sorted_by_x = points_array[points_array[:, 0].argsort()]
                idx = int(len(sorted_by_x) * q)
                candidate_x, candidate_y = sorted_by_x[idx]

                if is_in_bounds(candidate_x, candidate_y) and not is_overlapping(candidate_x, candidate_y):
                    return candidate_x, candidate_y

        # 尝试 3：八方向偏移
        offsets = [
            (0, label_h * 1.5),
            (0, -label_h * 1.5),
            (label_w * 1.2, 0),
            (-label_w * 1.2, 0),
            (label_w, label_h),
            (-label_w, label_h),
            (label_w, -label_h),
            (-label_w, -label_h),
        ]

        for dx, dy in offsets:
            new_x, new_y = init_x + dx, init_y + dy
            if is_in_bounds(new_x, new_y) and not is_overlapping(new_x, new_y):
                return new_x, new_y

        # 所有尝试失败，返回原始位置
        return init_x, init_y

    # ========================================================================
    # 伪二元相图：AMR 优化的并行计算
    # ========================================================================

    def _calculate_via_parallel(self,
                                standardized_base: Dict,
                                standardized_varying: str,
                                x_range: Tuple[float, float],
                                x_points: int,
                                T_range: Tuple[float, float],
                                T_points: int,
                                filtered_components: List[str],
                                filtered_phases: List[str],
                                model_spec: Optional[Dict],
                                pdens: int,
                                verbose: bool) -> np.ndarray:
        """
        [AMR 优化版] 自适应网格细化 + 多进程并行计算伪二元相图

        算法策略 (Coarse-to-Fine):
        ─────────────────────────
        1. 定义全分辨率网格 phase_map[T_points, x_points]
        2. 第一阶段：计算粗网格节点 (stride 间隔)
        3. 第二阶段：检测区块类型
           - 同质区块 (4角相同): 直接填充，跳过计算
           - 边界区块 (4角不同): 收集内部点待细化
        4. 第三阶段：批量计算边界区块内部点
        5. 返回完整填充的 phase_map
        """
        from multiprocessing import Pool, cpu_count

        start_time = time.time()

        # =====================================================================
        # 1. 初始化全尺寸网格（整数编码）
        # =====================================================================
        _PHASE_UNSET = -1
        _PHASE_ERROR = -2
        _PHASE_NO_PHASE = -3

        phase_map = np.full((T_points, x_points), _PHASE_UNSET, dtype=np.int32)
        phase_to_int = {}  # str -> int
        int_to_phase = {}  # int -> str
        next_phase_id = [0]

        def encode_phase(phase_str: str) -> int:
            """将相字符串编码为整数"""
            if phase_str == "Error":
                return _PHASE_ERROR
            if phase_str == "No Phase":
                return _PHASE_NO_PHASE
            if phase_str not in phase_to_int:
                phase_to_int[phase_str] = next_phase_id[0]
                int_to_phase[next_phase_id[0]] = phase_str
                next_phase_id[0] += 1
            return phase_to_int[phase_str]

        x_values = np.linspace(x_range[0], x_range[1], x_points)
        T_values = np.linspace(T_range[0], T_range[1], T_points)

        # 序列化数据库 (用于多进程传递)
        db_bytes = pickle.dumps(self.db)

        # =====================================================================
        # 2. AMR 参数配置
        # =====================================================================
        if max(T_points, x_points) >= 150:
            stride = 6
        elif max(T_points, x_points) >= 100:
            stride = 5
        else:
            stride = 4

        # 构建粗网格索引 (确保包含边界)
        t_indices = self._build_grid_indices(T_points, stride)
        x_indices = self._build_grid_indices(x_points, stride)

        if verbose:
            self.logger(f"  [AMR] 网格: {T_points}×{x_points} = {T_points * x_points} 总点数")
            self.logger(f"  [AMR] 步长: {stride}, 粗网格: {len(t_indices)}×{len(x_indices)}")

        # =====================================================================
        # 3. 第一阶段：粗网格计算
        # =====================================================================
        coarse_tasks = []
        for i in t_indices:
            for j in x_indices:
                task = (
                    i, j, T_values[i], x_values[j],
                    filtered_components, standardized_base, standardized_varying,
                    filtered_phases, model_spec, pdens
                )
                coarse_tasks.append(task)

        if verbose:
            self.logger(f"  [AMR] 第一阶段：粗网格计算 ({len(coarse_tasks)} 个点)...")

        # 执行粗网格计算（使用Pool initializer避免逐任务反序列化）
        n_workers = min(cpu_count(), max(4, len(coarse_tasks) // 10))
        coarse_results = self._execute_parallel_batch(coarse_tasks, n_workers, verbose, db_bytes=db_bytes)

        # 填入 phase_map（整数编码）
        for i, j, phase_str in coarse_results:
            phase_map[i, j] = encode_phase(phase_str)

        coarse_time = time.time() - start_time
        if verbose:
            self.logger(f"  [AMR] 粗网格完成，耗时 {coarse_time:.1f}s")

        # =====================================================================
        # 4. 第二阶段：区块检测 + 填充/标记（整数比较）
        # =====================================================================
        refine_coords = []
        filled_count = 0

        for i_idx in range(len(t_indices) - 1):
            for j_idx in range(len(x_indices) - 1):
                r0, r1 = t_indices[i_idx], t_indices[i_idx + 1]
                c0, c1 = x_indices[j_idx], x_indices[j_idx + 1]

                corners = [
                    phase_map[r0, c0],
                    phase_map[r0, c1],
                    phase_map[r1, c0],
                    phase_map[r1, c1],
                ]

                is_homogeneous = self._is_block_homogeneous(corners)

                if is_homogeneous:
                    fill_phase = corners[0]
                    for r in range(r0, r1 + 1):
                        for c in range(c0, c1 + 1):
                            if phase_map[r, c] == _PHASE_UNSET:
                                phase_map[r, c] = fill_phase
                                filled_count += 1
                else:
                    for r in range(r0, r1 + 1):
                        for c in range(c0, c1 + 1):
                            if phase_map[r, c] == _PHASE_UNSET:
                                refine_coords.append((r, c))

        if verbose:
            self.logger(f"  [AMR] 区块分析完成: {filled_count} 点直接填充, {len(refine_coords)} 点待细化")

        # =====================================================================
        # 5. 第三阶段：边界区块细化计算
        # =====================================================================
        if refine_coords:
            refine_tasks = []
            for r, c in refine_coords:
                task = (
                    r, c, T_values[r], x_values[c],
                    filtered_components, standardized_base, standardized_varying,
                    filtered_phases, model_spec, pdens
                )
                refine_tasks.append(task)

            if verbose:
                self.logger(f"  [AMR] 第三阶段：边界细化计算 ({len(refine_tasks)} 个点)...")

            refine_results = self._execute_parallel_batch(refine_tasks, n_workers, verbose, db_bytes=db_bytes)

            for i, j, phase_str in refine_results:
                phase_map[i, j] = encode_phase(phase_str)

        # =====================================================================
        # 6. 最终清理：处理遗留的未计算值
        # =====================================================================
        none_count = 0
        for i in range(T_points):
            for j in range(x_points):
                if phase_map[i, j] == _PHASE_UNSET:
                    phase_map[i, j] = self._fill_from_neighbors(phase_map, i, j, T_points, x_points)
                    none_count += 1

        if none_count > 0 and verbose:
            self.logger(f"  [AMR] 警告: {none_count} 个遗留点已用邻近值填充")

        # =====================================================================
        # 7. 统计与日志
        # =====================================================================
        total_time = time.time() - start_time
        total_points = T_points * x_points
        calculated_points = len(coarse_tasks) + len(refine_coords)
        skip_ratio = (total_points - calculated_points) / total_points * 100

        if verbose:
            self.logger(f"  [AMR] ═══════════════════════════════════════")
            self.logger(f"  [AMR] 总点数: {total_points}, 实际计算: {calculated_points}")
            self.logger(f"  [AMR] 跳过率: {skip_ratio:.1f}%, 总耗时: {total_time:.1f}s")
            self.logger(f"  [AMR] ═══════════════════════════════════════")

        return phase_map, int_to_phase

    # ========================================================================
    # 辅助方法：AMR 相关
    # ========================================================================

    def _build_grid_indices(self, n_points: int, stride: int) -> List[int]:
        """构建包含边界的网格索引列表"""
        indices = list(range(0, n_points, stride))
        if indices[-1] != n_points - 1:
            indices.append(n_points - 1)
        return indices

    def _is_block_homogeneous(self, corners) -> bool:
        """判断区块是否为同质区块（整数编码版本）"""
        _PHASE_UNSET = -1
        _PHASE_ERROR = -2
        _PHASE_NO_PHASE = -3
        if _PHASE_UNSET in corners:
            return False
        if _PHASE_ERROR in corners:
            return False
        if _PHASE_NO_PHASE in corners:
            return all(c == _PHASE_NO_PHASE for c in corners)
        return corners[0] == corners[1] == corners[2] == corners[3]

    def _execute_parallel_batch(self, tasks: List[Tuple], n_workers: int,
                                verbose: bool = False, db_bytes: bytes = None) -> List[Tuple]:
        """执行一批并行计算任务（使用 Pool initializer 优化）"""
        from multiprocessing import Pool

        if not tasks:
            return []

        try:
            chunksize = max(1, len(tasks) // (n_workers * 4))

            with Pool(n_workers, initializer=_worker_init, initargs=(db_bytes,)) as pool:
                results = pool.map(
                    _calculate_single_point_static,
                    tasks,
                    chunksize=chunksize
                )

            return list(results)

        except Exception as e:
            if verbose:
                self.logger(f"  [AMR] 并行计算失败: {e}, 降级为串行模式")

            # 串行回退：设置 _worker_state 后调用静态函数
            if db_bytes:
                _worker_state['db'] = pickle.loads(db_bytes)
            else:
                _worker_state['db'] = self.db

            results = []
            for task in tasks:
                try:
                    result = _calculate_single_point_static(task)
                    results.append(result)
                except Exception:
                    i, j = task[0], task[1]
                    results.append((i, j, "Error"))

            return results

    def _fill_from_neighbors(self, phase_map: np.ndarray, i: int, j: int,
                              n_rows: int, n_cols: int) -> int:
        """从邻近点填充缺失值（整数编码版本）"""
        neighbors = [
            (-1, 0), (1, 0), (0, -1), (0, 1),
            (-1, -1), (-1, 1), (1, -1), (1, 1),
        ]

        for di, dj in neighbors:
            ni, nj = i + di, j + dj
            if 0 <= ni < n_rows and 0 <= nj < n_cols:
                neighbor_val = phase_map[ni, nj]
                if neighbor_val >= 0:  # 有效相（非UNSET/ERROR/NO_PHASE）
                    return neighbor_val

        return -2  # PHASE_ERROR

    # ========================================================================
    # 辅助方法：相过滤与模型信息
    # ========================================================================

    def _filter_phases_strict(self, components: List[str]) -> List[str]:
        """严格的相过滤函数"""
        active_elements = set([c.upper() for c in components])
        if 'VA' not in active_elements:
            active_elements.add('VA')

        valid_phases = []

        for phase_name, phase_obj in self.db.phases.items():
            is_phase_possible = True

            for sublattice in phase_obj.constituents:
                sublattice_has_valid_species = False
                for species in sublattice:
                    species_elements = set(species.constituents.keys())
                    if species_elements.issubset(active_elements):
                        sublattice_has_valid_species = True
                        break

                if not sublattice_has_valid_species:
                    is_phase_possible = False
                    break

            if is_phase_possible:
                valid_phases.append(phase_name)

        return sorted(valid_phases)

    def _is_true_binary(self, base_composition: Dict, varying_comp: str) -> bool:
        """检测是否为真二元系统"""
        if not base_composition:
            return True

        non_zero_base = {k: val for k, val in base_composition.items() if val > 1e-6}

        if len(non_zero_base) == 0:
            return True
        elif len(non_zero_base) == 1:
            return True
        else:
            return False



# ============================================================================
# 液相面投影计算器
# ============================================================================

class LiquidusSurfaceCalculator(_ContributionCoefficientMixin):
	"""
	液相面投影计算器 - 返回统一的CalculationResult
	"""
	
	def __init__ (self, database: Database, logger: Optional[Callable] = None):
		self.db = database
		self.logger = logger or (lambda msg: None)
		self.ls_calculator = LiquidusSolidusCalculator(database, logger)
	
	def calculate (self,
	               comp_b_range: Tuple[float, float],
	               comp_c_range: Tuple[float, float],
	               grid_points: int,
	               T_range: Tuple[float, float],
	               components: List[str],
	               phases: Optional[List[str]] = None,
	               model_spec: Optional[Dict] = None,
	               pdens: int = 2000,
	               verbose: bool = False) -> CalculationResult:
		"""
		计算液相面投影

		Returns:
		--------
		CalculationResult with mode='surface'
		"""
		try:
			components = [c.upper() for c in components]
			# 提取三元组分
			non_va_comps = [c for c in components if c != 'VA']
			if len(non_va_comps) < 3:
				return CalculationResult(
						success=False,
						mode=CalculationMode.SURFACE.value,
						message="需要至少3个非VA组分"
				)
			
			comp_a, comp_b, comp_c = non_va_comps[:3]

			if verbose:
				self.logger(f"三元系统: {comp_a}-{comp_b}-{comp_c}")

			# ★ 打印贡献系数（如果使用非默认模型）
			if verbose and model_spec:
				# 使用温度范围的中点作为代表温度
				mid_temp = (T_range[0] + T_range[1]) / 2
				self._print_contribution_coefficients(non_va_comps[:3], mid_temp, model_spec)

			# 生成三元网格
			b_step = 1.0 / (grid_points - 1)
			c_step = 1.0 / (grid_points - 1)
			b_coords, c_coords = np.meshgrid(
					np.arange(0, 1 + b_step, b_step),
					np.arange(0, 1 + c_step, c_step)
			)
			
			# 过滤B+C > 1的点
			valid_mask = b_coords + c_coords <= 1.0001
			b_valid = b_coords[valid_mask]
			c_valid = c_coords[valid_mask]
			
			total_points = len(b_valid)
			if verbose:
				self.logger(f"液相面计算: {total_points}个有效网格点")
			
			# 构建任务列表
			t_data = np.full(total_points, np.nan)
			tasks = []
			for i in range(total_points):
				xb = b_valid[i]
				xc = c_valid[i]
				xa = 1.0 - xb - xc
				composition = self._build_ternary_composition(comp_a, xa, comp_b, xb, comp_c, xc)
				tasks.append((i, composition, T_range, components, phases,
				              model_spec, pdens))

			# 多进程并行计算
			n_workers = max(1, os.cpu_count() - 1) if os.cpu_count() else 1
			use_parallel = total_points >= 4 and n_workers > 1

			if use_parallel:
				try:
					from multiprocessing import Pool
					db_bytes = pickle.dumps(self.db)
					chunksize = max(1, total_points // (n_workers * 4))

					if verbose:
						self.logger(f"液相面并行计算: {n_workers}进程, {total_points}任务")

					with Pool(n_workers, initializer=_surface_worker_init,
					          initargs=(db_bytes,)) as pool:
						results = pool.map(_surface_calculate_point, tasks,
						                   chunksize=chunksize)

					for idx, t_liq in results:
						if t_liq is not None and not np.isnan(t_liq):
							t_data[idx] = t_liq

				except Exception as e:
					if verbose:
						self.logger(f"  并行失败: {e}, 降级为串行模式")
					use_parallel = False

			if not use_parallel:
				# 串行回退
				_worker_state['db'] = self.db
				_worker_state['ls_calc'] = LiquidusSolidusCalculator(self.db, self.logger)
				for task in tasks:
					if verbose and task[0] % max(1, total_points // 10) == 0:
						self.logger(f"  进度: {task[0]}/{total_points}")
					try:
						idx, t_liq = _surface_calculate_point(task)
						if t_liq is not None and not np.isnan(t_liq):
							t_data[idx] = t_liq
					except Exception:
						pass

			if verbose:
				valid_count = np.sum(~np.isnan(t_data))
				self.logger(f"✓ 完成，有效数据点: {valid_count}/{total_points}")
			
			return CalculationResult(
					success=True,
					mode=CalculationMode.SURFACE.value,
					message=f"液相面投影计算完成: {np.sum(~np.isnan(t_data))}/{total_points}有效点",
					x_axis=b_valid,
					y_axis=c_valid,
					z_axis=t_data,
					surface_data={
						'xb_data': b_valid,
						'xc_data': c_valid,
						't_data': t_data,
						'comp_a': comp_a,
						'comp_b': comp_b,
						'comp_c': comp_c,
						'grid_points': grid_points
					}
			)
		
		except Exception as e:
			self.logger(f"❌ 液相面计算失败: {e}")
			if verbose:
				self.logger(traceback.format_exc())
			return CalculationResult(
					success=False,
					mode=CalculationMode.SURFACE.value,
					message=f"液相面计算失败: {str(e)}"
			)
	
	def _build_ternary_composition (self, comp_a: str, xa: float,
	                                comp_b: str, xb: float,
	                                comp_c: str, xc: float,
	                                min_fraction: float = 1e-9) -> Dict:
		"""构建三元成分字典"""
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


# ============================================================================
# 热力学性质计算器
# ============================================================================

class ThermodynamicPropertiesCalculator(_ContributionCoefficientMixin):
	"""
	热力学性质计算器 - 返回统一的CalculationResult

	功能：计算液相中的Gibbs自由能、混合焓和活度
	"""
	
	def __init__ (self, database: Database, logger: Optional[Callable] = None):
		self.db = database
		self.logger = logger or (lambda msg: None)
	
	def calculate (self,
	               base_composition: Dict,
	               varying_comp: str,
	               x_range: Tuple[float, float],
	               x_points: int,
	               temperature: float,
	               components: List[str],
	               phases: Optional[List[str]] = None,
	               model_spec: Optional[Dict] = None,
	               pdens: int = 2000,
	               verbose: bool = False) -> CalculationResult:
		"""
		计算热力学性质

		Returns:
		--------
		CalculationResult with mode='properties'
		"""
		try:
			base_composition = {k.upper(): val for k, val in base_composition.items()}
			varying_comp = varying_comp.upper()
			components = [c.upper() for c in components]
			x_values = np.linspace(x_range[0], x_range[1], x_points)
			active_phases = phases if phases else list(self.db.phases.keys())
			
			# 找到液相
			liquid_phase = self._find_liquid_phase(active_phases)
			if liquid_phase is None:
				return CalculationResult(
						success=False,
						mode=CalculationMode.PROPERTIES.value,
						message="未找到LIQUID相"
				)
			
			# 构建成分列表
			compositions = []
			for x in x_values:
				comp = self._build_composition(base_composition, varying_comp, x)
				compositions.append(comp)
			
			# 确定所有组分
			all_comps = set()
			for comp in compositions:
				all_comps.update(comp.keys())

			# ★ 打印贡献系数（如果使用非默认模型）
			if verbose and model_spec:
				components_no_va = [c for c in all_comps if c != 'VA']
				self._print_contribution_coefficients(components_no_va, temperature, model_spec)

			# 计算参考化学势
			ref_mu = self._calculate_reference_potentials(
					all_comps, liquid_phase, temperature, model_spec, pdens, verbose
			)
			
			if verbose:
				self.logger(f"参考化学势: {ref_mu}")
				self.logger(f"计算热力学性质: {len(compositions)}个成分点")
			
			# 批量计算
			gibbs_list = []
			enthalpy_list = []
			activity_dict = {c: [] for c in all_comps}

			RT = 8.314 * temperature

			# 检测二元体系以使用向量化路径
			non_va_comps = sorted([c for c in all_comps if c != 'VA'])
			vectorized_ok = False

			if len(non_va_comps) == 2:
				try:
					active_comps = sorted(list(set(list(all_comps) + ['VA'])))
					conds = {v.T: temperature, v.P: 101325, v.N: 1,
					         v.X(varying_comp): x_values}

					eq = equilibrium(
							self.db, active_comps, [liquid_phase], conds,
							model=model_spec,
							output=['GM', 'HM_MIX', 'MU'],
							calc_opts={'pdens': pdens}
					)

					gm_vals = np.squeeze(eq.GM.values)
					hm_vals = np.squeeze(eq.HM_MIX.values)

					if gm_vals.ndim == 0:
						gm_vals = np.array([float(gm_vals)])
						hm_vals = np.array([float(hm_vals)])

					gibbs_list = [float(g) for g in gm_vals.flatten()[:len(x_values)]]
					enthalpy_list = [float(h) for h in hm_vals.flatten()[:len(x_values)]]

					for c in activity_dict.keys():
						if c in ref_mu and not np.isnan(ref_mu[c]) and c in active_comps:
							try:
								mu_vals = np.squeeze(eq.MU.sel(component=c).values)
								if mu_vals.ndim == 0:
									mu_vals = np.array([float(mu_vals)])
								mu_flat = mu_vals.flatten()[:len(x_values)]
								activity_dict[c] = [np.exp((float(mu) - ref_mu[c]) / RT) for mu in mu_flat]
							except Exception:
								activity_dict[c] = [np.nan] * len(x_values)
						else:
							activity_dict[c] = [np.nan] * len(x_values)

					vectorized_ok = True
					if verbose:
						self.logger(f"  二元向量化计算完成: {len(x_values)}点")
				except Exception as e:
					if verbose:
						self.logger(f"  向量化失败: {e}, 回退到逐点模式")
					gibbs_list = []
					enthalpy_list = []
					activity_dict = {c: [] for c in all_comps}

			if not vectorized_ok:
				for i, comp in enumerate(compositions):
					if verbose and i % max(1, len(compositions) // 10) == 0:
						self.logger(f"  进度: {i}/{len(compositions)}")

					conds = self._build_conditions(comp, temperature)
					active_comps = sorted(list(comp.keys()) + ['VA'])

					try:
						eq = equilibrium(
								self.db, active_comps, [liquid_phase], conds,
								model=model_spec,
								output=['GM', 'HM_MIX', 'MU'],
								calc_opts={'pdens': pdens}
						)

						gibbs = float(eq.GM.values.flatten()[0])
						enthalpy = float(eq.HM_MIX.values.flatten()[0])
						gibbs_list.append(gibbs)
						enthalpy_list.append(enthalpy)

						if verbose:
							comp_str = ', '.join([f'{k}={val:.3f}' for k, val in comp.items()])
							self.logger(f"  [{i + 1}/{len(compositions)}] {comp_str}: G={gibbs:.2e}, H_mix={enthalpy:.2e}")

						for c in activity_dict.keys():
							if c in ref_mu and not np.isnan(ref_mu[c]) and c in active_comps:
								try:
									mu = float(eq.MU.sel(component=c).values.flatten()[0])
									activity = np.exp((mu - ref_mu[c]) / RT)
									activity_dict[c].append(activity)
								except Exception:
									activity_dict[c].append(np.nan)
							else:
								activity_dict[c].append(np.nan)

					except Exception as e:
						if verbose:
							self.logger(f"  计算失败: {e}")
						gibbs_list.append(np.nan)
						enthalpy_list.append(np.nan)
						for c in activity_dict.keys():
							activity_dict[c].append(np.nan)

			activity_arrays = {c: np.array(vals) for c, vals in activity_dict.items()}
			
			return CalculationResult(
					success=True,
					mode=CalculationMode.PROPERTIES.value,
					message=f"热力学性质计算完成: {x_points}个成分点 @ {temperature}K",
					x_axis=x_values,
					y_axis=np.array(gibbs_list),  # GM as primary y_axis
					properties_data={
						'compositions': compositions,
						'temperature': temperature,
						'x_values': x_values,
						'varying_comp': varying_comp,
						'gibbs_energy': np.array(gibbs_list),
						'enthalpy_mix': np.array(enthalpy_list),
						'activity': activity_arrays,
						'reference_mu': ref_mu
					}
			)
		
		except Exception as e:
			if verbose:
				self.logger(traceback.format_exc())
			return CalculationResult(
					success=False,
					mode=CalculationMode.PROPERTIES.value,
					message=f"热力学性质计算失败: {str(e)}"
			)
	
	def _find_liquid_phase (self, phases: List[str]) -> Optional[str]:
		"""查找液相名称"""
		for ph in phases:
			if 'LIQUID' in ph.upper():
				return ph
		return None
	
	def _calculate_reference_potentials (self, components: set, liquid_phase: str,
	                                     temperature: float, model_spec: Optional[Dict],
	                                     pdens: int, verbose: bool) -> Dict:
		"""计算纯组分的参考化学势"""
		ref_mu = {}
		for comp in components:
			try:
				conds = {v.T: temperature, v.P: 101325, v.N: 1}
				eq = equilibrium(
						self.db, [comp, 'VA'], [liquid_phase], conds,
						model=model_spec, output='MU', calc_opts={'pdens': pdens}
				)
				ref_mu[comp] = float(eq.MU.sel(component=comp).values.flatten()[0])
			except Exception as e:
				if verbose:
					self.logger(f"  参考化学势计算失败 ({comp}): {e}")
				ref_mu[comp] = np.nan
		return ref_mu
	
	def _build_conditions (self, composition: Dict, temperature: float) -> Dict:
		"""构建平衡计算条件"""
		conds = {v.T: temperature, v.P: 101325, v.N: 1}
		
		comp_list = list(composition.keys())
		n_comps = len(comp_list)
		
		if n_comps > 1:
			total = sum(composition.values())
			for comp in comp_list[:-1]:
				conds[v.X(comp)] = composition[comp] / total
		
		return conds
	
	def _build_composition (self, base_dict: Dict, varying_comp: str, x_value: float) -> Dict:
		"""构建特定x值的成分"""
		base_total = sum(base_dict.values()) if base_dict else 1.0

		composition = {varying_comp: x_value}
		remaining = 1.0 - x_value

		for elem, frac in base_dict.items():
			base_fraction = frac / base_total
			contribution = base_fraction * remaining
			composition[elem] = composition.get(elem, 0.0) + contribution

		return composition


# ============================================================================
# 三元相图计算器（可选）
# ============================================================================

class TernaryCalculator(_ContributionCoefficientMixin):
	"""
	三元相图计算器 - 返回统一的CalculationResult

	使用PyCalphad的ternplot绘制等温截面
	"""
	
	def __init__ (self, database: Database, logger: Optional[Callable] = None):
		self.db = database
		self.logger = logger or (lambda msg: None)
	
	def calculate (self,
	               components: List[str],
	               temperature: float,
	               phases: Optional[List[str]] = None,
	               model_spec: Optional[Dict] = None,
	               pdens: int = 2000,
	               grid_points: int = 30,
	               verbose: bool = False) -> CalculationResult:
		"""
		计算三元等温截面相图
 
		Returns:
		--------
		CalculationResult with mode='ternary', figure=matplotlib.Figure
		"""
		try:
			from pycalphad.mapping.compat_api import ternplot
			
			components = [c.upper() for c in components]
			
			# 确保只有3个非VA组分
			non_va_comps = [c for c in components if c.upper() != 'VA']
			if len(non_va_comps) != 3:
				return CalculationResult(
						success=False,
						mode=CalculationMode.TERNARY.value,
						message=f"三元相图需要恰好3个非VA组分，得到: {non_va_comps}"
				)
			
			comps_with_va = sorted(non_va_comps) + ['VA']
			active_phases = phases if phases else list(self.db.phases.keys())

			if verbose:
				self.logger(f"计算三元相图: {non_va_comps[0]}-{non_va_comps[1]}-{non_va_comps[2]} @ {temperature}K")

			# ★ 打印贡献系数（如果使用非默认模型）
			if verbose and model_spec:
				self._print_contribution_coefficients(non_va_comps, temperature, model_spec)

			# 设置条件（第三个值是步长，不是网格点数）
			step_size = 1.0 / max(grid_points - 1, 1)
			conditions = {
				v.T: temperature,
				v.P: 101325,
				v.N: 1,
				v.X(comps_with_va[0]): (0, 1, step_size),
				v.X(comps_with_va[1]): (0, 1, step_size)
			}
			
			# 构建 map_kwargs（传给 TernaryStrategy，包含 model）
			map_kwargs = {}
			if model_spec:
				map_kwargs['model'] = model_spec

			# 调用 ternplot：让它自行创建三角投影坐标轴，返回策略对象
			ax_result, strategy = ternplot(
					self.db, comps_with_va, active_phases, conditions,
					map_kwargs=map_kwargs if map_kwargs else None,
					return_strategy=True
			)

			# 从返回的 axes 获取 figure
			fig = ax_result.get_figure()
			fig.set_size_inches(10, 8)

			if verbose:
				self.logger(f"✓ 三元相图计算完成")

			# 添加相区标注
			self._add_phase_labels_to_ternary(
				ax_result, strategy, non_va_comps, temperature,
				active_phases, model_spec, pdens, verbose
			)

			return CalculationResult(
					success=True,
					mode=CalculationMode.TERNARY.value,
					message=f"三元相图 {non_va_comps[0]}-{non_va_comps[1]}-{non_va_comps[2]} @ {temperature}K 计算完成",
					figure=fig,
					components=non_va_comps
			)
		
		except Exception as e:
			if verbose:
				self.logger(traceback.format_exc())
			return CalculationResult(
					success=False,
					mode=CalculationMode.TERNARY.value,
					message=f"三元相图计算失败: {str(e)}"
			)

	def _add_phase_labels_to_ternary(self, ax, strategy, non_va_comps,
	                                  temperature, phases, model_spec,
	                                  pdens, verbose,
	                                  n_grid=13, min_points=2):
		"""在三元相图上添加相区标注

		通过粗网格平衡采样识别各区域稳定相组合，在质心处放置标签。

		Parameters
		----------
		ax : matplotlib Axes (triangular projection)
		strategy : TernaryStrategy
		non_va_comps : list of 3 component names
		temperature : float
		phases : list of phase names
		model_spec : model specification or None
		pdens : int
		verbose : bool
		n_grid : int
			采样网格密度（每条边的分割数）
		min_points : int
			标注所需最少采样点数
		"""
		try:
			start_time = time.time()
			comps = sorted(non_va_comps) + ['VA']
			comp_x = comps[0]  # 第一个组分 = x 轴
			comp_y = comps[1]  # 第二个组分 = y 轴

			calc_opts = {'pdens': min(500, pdens)}
			phase_regions = defaultdict(list)
			margin = 0.02
			success_count = 0

			if verbose:
				self.logger(f"  [标注] 三元相区采样: 网格 {n_grid}×{n_grid}")

			# 粗网格采样（三角形内部）
			x1_vals = np.linspace(margin, 1.0 - 2 * margin, n_grid)

			for x1_val in x1_vals:
				x2_max = 1.0 - x1_val - margin
				if x2_max < margin:
					continue

				n_x2 = max(2, int(n_grid * x2_max))
				x2_vals = np.linspace(margin, x2_max, n_x2)

				try:
					# 向量化：一次 equilibrium 调用多个 x2 值
					conds = {
						v.T: temperature,
						v.P: 101325,
						v.N: 1,
						v.X(comp_x): x1_val,
						v.X(comp_y): x2_vals,
					}

					if model_spec:
						eq = equilibrium(self.db, comps, phases, conds,
						                 model=model_spec, calc_opts=calc_opts)
					else:
						eq = equilibrium(self.db, comps, phases, conds,
						                 calc_opts=calc_opts)

					phase_values = np.squeeze(eq.Phase.values)
					np_values = np.squeeze(eq.NP.values)

					if phase_values.ndim == 1:
						phase_values = phase_values.reshape(1, -1)
						np_values = np_values.reshape(1, -1)

					for j, x2_val in enumerate(x2_vals):
						try:
							stable = []
							for ph, amt in zip(phase_values[j], np_values[j]):
								if isinstance(ph, bytes):
									ph = ph.decode('utf-8')
								ph_str = str(ph).strip()
								try:
									amt_f = float(amt)
								except (TypeError, ValueError):
									continue
								if ph_str and ph_str != '' and amt_f > 1e-3:
									# 去除混溶间隙后缀 #N
									base_name = ph_str.split('#')[0]
									stable.append(base_name)
							if stable:
								phase_str = " + ".join(sorted(set(stable)))
								phase_regions[phase_str].append((x1_val, x2_val))
								success_count += 1
						except Exception:
							continue

				except Exception:
					# 向量化失败，逐点回退
					for x2_val in x2_vals:
						try:
							conds = {
								v.T: temperature, v.P: 101325, v.N: 1,
								v.X(comp_x): x1_val, v.X(comp_y): float(x2_val),
							}
							if model_spec:
								eq = equilibrium(self.db, comps, phases, conds,
								                 model=model_spec, calc_opts=calc_opts)
							else:
								eq = equilibrium(self.db, comps, phases, conds,
								                 calc_opts=calc_opts)
							phase_str = PseudoBinaryCalculator._extract_stable_phases(eq)
							if phase_str and phase_str not in ("Error", "No Phase"):
								phase_regions[phase_str].append((x1_val, x2_val))
								success_count += 1
						except Exception:
							continue

			sample_time = time.time() - start_time
			if verbose:
				self.logger(f"  [标注] 采样完成: {success_count} 点, "
				            f"{len(phase_regions)} 个相区, {sample_time:.1f}s")

			# 按区域大小排序（大区域先标注）
			sorted_regions = sorted(phase_regions.items(), key=lambda x: -len(x[1]))

			# 标签尺寸估算（数据坐标空间）
			label_w = 0.10
			label_h = 0.06
			occupied = []
			labeled_count = 0

			for phase_str, points in sorted_regions:
				if len(points) < min_points:
					continue

				pts = np.array(points)
				cx = np.mean(pts[:, 0])
				cy = np.mean(pts[:, 1])

				# 窄区域用中位数
				x_spread = pts[:, 0].max() - pts[:, 0].min()
				y_spread = pts[:, 1].max() - pts[:, 1].min()
				if x_spread < 0.08:
					cx = np.median(pts[:, 0])
				if y_spread < 0.08:
					cy = np.median(pts[:, 1])

				# 确保标签在三角形内
				cx = np.clip(cx, margin + label_w / 2, 1.0 - margin - label_w / 2)
				cy = np.clip(cy, margin + label_h / 2, 1.0 - cx - margin - label_h / 2)

				# 避免重叠
				best_x, best_y = cx, cy
				overlap = False
				for ox, oy, ow, oh in occupied:
					if (abs(best_x - ox) < (label_w + ow) / 2 and
					    abs(best_y - oy) < (label_h + oh) / 2):
						overlap = True
						break

				if overlap:
					# 尝试偏移
					offsets = [
						(0, label_h * 1.5), (0, -label_h * 1.5),
						(label_w * 1.2, 0), (-label_w * 1.2, 0),
						(label_w, label_h), (-label_w, -label_h),
					]
					for dx, dy in offsets:
						nx, ny = cx + dx, cy + dy
						if nx < margin or ny < margin or nx + ny > 1.0 - margin:
							continue
						no_overlap = True
						for ox, oy, ow, oh in occupied:
							if (abs(nx - ox) < (label_w + ow) / 2 and
							    abs(ny - oy) < (label_h + oh) / 2):
								no_overlap = False
								break
						if no_overlap:
							best_x, best_y = nx, ny
							break

				# 长名称换行
				display_name = phase_str
				if len(display_name) > 15 and '+' in display_name:
					display_name = display_name.replace(' + ', '\n+ ')

				# 根据区域大小调整字体
				fontsize = 8 if len(points) < 8 else (9 if len(points) < 20 else 10)

				ax.text(best_x, best_y, display_name,
				        fontsize=fontsize, fontweight='bold',
				        ha='center', va='center',
				        bbox=dict(boxstyle='round,pad=0.3',
				                  facecolor='white', alpha=0.85,
				                  edgecolor='gray', linewidth=0.8),
				        zorder=100)

				occupied.append((best_x, best_y, label_w, label_h))
				labeled_count += 1

				if verbose:
					self.logger(f"  [标注] ✓ '{phase_str}' @ ({best_x:.2f}, {best_y:.2f}), "
					            f"{len(points)} 点")

			total_time = time.time() - start_time
			if verbose:
				self.logger(f"  [标注] 完成: {labeled_count} 个标注, 总耗时 {total_time:.1f}s")

		except Exception as e:
			if verbose:
				self.logger(f"  [标注] 错误: {e}")
				self.logger(f"  {traceback.format_exc()}")


# ============================================================================
# 多进程 Worker 函数（模块级别，供 Pool 使用）
# ============================================================================

_worker_state = {}


def _worker_init(db_bytes):
	"""Pool initializer: 每个worker进程只反序列化一次数据库"""
	_worker_state['db'] = pickle.loads(db_bytes)


def _calculate_single_point_static(args):
	"""计算单个网格点的稳定相组合（使用预加载的db，减少序列化开销）"""
	(i, j, T, x, filtered_components, standardized_base,
	 standardized_varying, filtered_phases, model_spec, pdens) = args

	from pycalphad import equilibrium, variables as v

	try:
		db = _worker_state['db']

		# 构建成分：varying组分设为x，剩余(1-x)按基础合金中其他组分比例分配
		composition = {standardized_varying: x}
		remaining = 1.0 - x

		base_elements = {k: val for k, val in standardized_base.items()
		                 if k != standardized_varying}
		base_total = sum(base_elements.values()) if base_elements else 1.0

		for elem, frac in base_elements.items():
			composition[elem] = (frac / base_total) * remaining

		# 过滤有效组分
		non_zero_comps = [c for c in filtered_components
		                  if c == 'VA' or composition.get(c, 0) > 1e-6]

		non_zero_comp = {k: val for k, val in composition.items() if val > 1e-6}
		comp_list = list(non_zero_comp.keys())

		# 设置条件
		conditions = {v.T: T, v.P: 101325, v.N: 1}

		if len(comp_list) > 1:
			total = sum(non_zero_comp.values())
			for comp in comp_list[:-1]:
				conditions[v.X(comp)] = non_zero_comp[comp] / total

		# 平衡计算
		try:
			calc_opts = {'pdens': pdens}

			if model_spec:
				eq = equilibrium(db, non_zero_comps, filtered_phases,
				                 conditions=conditions, model=model_spec,
				                 calc_opts=calc_opts)
			else:
				eq = equilibrium(db, non_zero_comps, filtered_phases,
				                 conditions=conditions, calc_opts=calc_opts)

			phase_str = PseudoBinaryCalculator._extract_stable_phases(eq)
			return (i, j, phase_str)

		except Exception:
			return (i, j, "Error")

	except Exception:
		return (i, j, "Error")


def _surface_worker_init(db_bytes):
	"""LiquidusSurfaceCalculator Pool initializer"""
	_worker_state['db'] = pickle.loads(db_bytes)
	_worker_state['ls_calc'] = LiquidusSolidusCalculator(_worker_state['db'])


def _surface_calculate_point(args):
	"""计算单个液相面网格点的液相线温度"""
	(i, composition, T_range, components, phases,
	 model_spec, pdens) = args

	try:
		ls_calc = _worker_state['ls_calc']
		result = ls_calc.calculate(
			composition=composition,
			T_range=T_range,
			components=components,
			phases=phases,
			model_spec=model_spec,
			pdens=pdens,
			T_step=10.0,
			refine=False,
			use_bisection=True,
			bisection_tol=1.0,
			verbose=False
		)
		if result.success and result.liquidus_data:
			t_liq = result.liquidus_data.get('liquidus_K')
			return (i, t_liq if t_liq is not None else np.nan)
		return (i, np.nan)
	except Exception:
		return (i, np.nan)


# ============================================================================
# 导出
# ============================================================================

__all__ = [
	# 数据契约
	'CalculationResult',
	'CalculationMode',

	# 数据管理
	'DatabaseManager',

	# 计算器
	'LiquidusSolidusCalculator',
	'SolubilityCalculator',
	'PseudoBinaryCalculator',
	'LiquidusSurfaceCalculator',
	'ThermodynamicPropertiesCalculator',
	'TernaryCalculator',
]