#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
LLM 提示词模板与响应解析器
============================

为合金热力学计算提供5组领域特定的提示词模板：
1. 相预测 (predict_phases)
2. 成分解析 (parse_composition)
3. 参数建议 (suggest_parameters)
4. 结果解读 (interpret_result)
5. 模型推荐 (recommend_model)

每个模板配套一个防御性解析器，提取JSON并验证合法性。
"""

import json
import re
from typing import Dict, List, Optional, Any


# ============================================================================
# 1. 相预测 (Phase Prediction)
# ============================================================================

PHASE_PREDICTION_SYSTEM = """你是一位计算热力学（CALPHAD）专家。给定合金组分和热力学数据库中的候选相列表，
预测哪些相最可能参与平衡计算。

规则：
- 始终包含 LIQUID 相（如果候选列表中有）
- 包含与组分相关的标准溶液相（FCC_A1, BCC_A2, HCP_A3 等）
- 包含已知的金属间化合物相（如 AL3NI, MG2SI 等）
- 排除有序亚晶格变体（除非特别相关）
- 排除明显不相关的相

返回格式：仅返回JSON数组，如 ["LIQUID", "FCC_A1", "AL3NI"]
不要添加任何解释文字。"""


def build_phase_prediction_prompt(components: List[str], all_phases: List[str],
                                   calc_type: str = "general") -> str:
	"""构建相预测用户提示词"""
	calc_type_map = {
		"liquidus_solidus": "液相线/固相线计算",
		"solubility": "溶解度计算",
		"pseudo_binary": "伪二元相图计算",
		"ternary": "三元等温截面计算",
		"surface": "液相面投影计算",
		"properties": "热力学性质计算",
		"general": "通用平衡计算",
	}
	calc_desc = calc_type_map.get(calc_type, calc_type)

	return f"""合金组分: {', '.join(components)}
计算类型: {calc_desc}
数据库中的候选相列表: {', '.join(all_phases)}

请从候选列表中选择最可能参与此计算的相。仅返回JSON数组。"""


def parse_phase_prediction(raw_text: str, valid_phases: List[str]) -> Optional[List[str]]:
	"""解析相预测响应，验证每个相名存在于候选列表中"""
	try:
		# 提取JSON数组
		json_match = re.search(r'\[.*?\]', raw_text, re.DOTALL)
		if not json_match:
			return None

		phases = json.loads(json_match.group())
		if not isinstance(phases, list):
			return None

		# 构建大小写不敏感的映射
		valid_map = {p.upper(): p for p in valid_phases}

		result = []
		for phase in phases:
			if not isinstance(phase, str):
				continue
			p_upper = phase.strip().upper()
			if p_upper in valid_map:
				result.append(valid_map[p_upper])

		return result if result else None
	except (json.JSONDecodeError, Exception):
		return None


# ============================================================================
# 2. 成分解析 (Composition Parsing)
# ============================================================================

COMPOSITION_PARSE_SYSTEM = """你是一位材料科学专家。将用户输入的合金成分描述解析为结构化的元素比例格式。

支持的输入格式：
- 标准记法: "Al-7Si-0.3Mg" (重量百分比，基元素补余)
- 合金牌号: "6061", "Ti-6Al-4V", "Inconel 718", "AA356"
- 化学式: "Al0.7Cu0.3", "AL1CU1"
- 原子百分比: "Al 90at% Cu 10at%"
- 重量百分比: "Al-4wt%Cu"
- 中文描述: "铝硅合金 含7%硅"

返回格式（仅JSON，无其他文字）：
{
  "composition": {"元素符号": 比例数值, ...},
  "type": "mole_ratio" 或 "weight_pct" 或 "atom_pct",
  "confidence": 0.0到1.0
}

注意：
- 元素符号必须大写（如 AL, CU, SI, MG, ZN, FE, NI, TI, CR）
- 如果是摩尔比/原子比，直接给出比例（如 AL1CU0.5）
- 如果是重量百分比，给出百分比数值（如 AL93 SI7 表示 93wt%Al, 7wt%Si）
- 对于合金牌号，给出典型的标称成分（重量百分比）"""


def build_composition_prompt(user_input: str, available_comps: List[str]) -> str:
	"""构建成分解析用户提示词"""
	return f"""用户输入: "{user_input}"
数据库中可用的元素: {', '.join(available_comps)}

请解析此合金成分。仅返回JSON。"""


def parse_composition_response(raw_text: str,
                                available_comps: List[str]) -> Optional[Dict[str, float]]:
	"""解析成分解析响应，验证元素存在于可用列表中"""
	try:
		# 提取JSON对象
		json_match = re.search(r'\{.*\}', raw_text, re.DOTALL)
		if not json_match:
			return None

		data = json.loads(json_match.group())

		# 提取 composition 字段
		comp = data.get('composition')
		if not comp or not isinstance(comp, dict):
			return None

		# 构建大小写不敏感的映射
		valid_map = {c.upper(): c.upper() for c in available_comps}

		result = {}
		for elem, value in comp.items():
			elem_upper = elem.strip().upper()
			if elem_upper in valid_map and elem_upper not in ('VA', 'VACUUM'):
				try:
					val = float(value)
					if val > 0:
						result[valid_map[elem_upper]] = val
				except (ValueError, TypeError):
					continue

		return result if result else None
	except (json.JSONDecodeError, Exception):
		return None


# ============================================================================
# 3. 参数建议 (Parameter Suggestion)
# ============================================================================

PARAMETER_SUGGEST_SYSTEM = """你是一位计算热力学专家。给定合金组分和计算类型，建议最优的计算参数。

你需要根据合金体系的特点给出合理的参数建议：
- 温度范围应覆盖所有关键相变（熔化、固溶、析出等）
- 温度步长应在相变区域足够精细
- 网格点数应平衡精度与计算速度

返回格式（仅JSON，无其他文字）：
{
  "temp_min": 温度下限(K),
  "temp_max": 温度上限(K),
  "temp_step": 温度步长(K),
  "scan_points": 扫描点数(整数),
  "reasoning": "简要说明(中文)"
}"""


def build_parameter_prompt(components: List[str], calc_type: str,
                            current_params: Optional[Dict] = None) -> str:
	"""构建参数建议用户提示词"""
	calc_type_map = {
		"liquidus_solidus": "液相线/固相线计算",
		"solubility": "溶解度曲线计算",
		"pseudo_binary": "伪二元相图计算",
		"ternary": "三元等温截面计算",
		"surface": "液相面投影计算",
		"properties": "热力学性质计算",
	}
	calc_desc = calc_type_map.get(calc_type, calc_type)

	prompt = f"""合金组分: {', '.join(components)}
计算类型: {calc_desc}"""

	if current_params:
		prompt += f"\n当前参数: {json.dumps(current_params, ensure_ascii=False)}"

	prompt += "\n\n请建议最优的计算参数。仅返回JSON。"
	return prompt


def parse_parameter_response(raw_text: str) -> Optional[Dict[str, Any]]:
	"""解析参数建议响应，验证数值范围合理"""
	try:
		json_match = re.search(r'\{.*\}', raw_text, re.DOTALL)
		if not json_match:
			return None

		data = json.loads(json_match.group())

		result = {}

		# 温度下限 (200K - 5000K)
		if 'temp_min' in data:
			val = float(data['temp_min'])
			if 200 <= val <= 5000:
				result['temp_min'] = val

		# 温度上限
		if 'temp_max' in data:
			val = float(data['temp_max'])
			if 200 <= val <= 5000:
				result['temp_max'] = val

		# 验证 min < max
		if 'temp_min' in result and 'temp_max' in result:
			if result['temp_min'] >= result['temp_max']:
				return None

		# 温度步长 (0.1K - 100K)
		if 'temp_step' in data:
			val = float(data['temp_step'])
			if 0.1 <= val <= 100:
				result['temp_step'] = val

		# 扫描点数 (5 - 500)
		if 'scan_points' in data:
			val = int(data['scan_points'])
			if 5 <= val <= 500:
				result['scan_points'] = val

		# 说明
		if 'reasoning' in data:
			result['reasoning'] = str(data['reasoning'])

		return result if result else None
	except (json.JSONDecodeError, ValueError, Exception):
		return None


# ============================================================================
# 4. 结果解读 (Result Interpretation)
# ============================================================================

RESULT_INTERPRET_SYSTEM = """你是一位冶金学专家，正在分析热力学计算结果。
请提供简洁、专业的中文解读，重点关注：

1. 关键相变温度及其意义
2. 计算结果是否物理合理
3. 可能的计算异常或伪影
4. 对合金设计的实际启示

保持回复在200字以内。使用中文。
如果检测到异常，请明确指出并给出可能的原因和建议。"""


def build_result_interpret_prompt(result_dict: Dict, context: Dict) -> str:
	"""构建结果解读用户提示词"""
	# 简化结果数据，避免过大
	simplified = {
		'success': result_dict.get('success'),
		'mode': result_dict.get('mode'),
		'message': result_dict.get('message', ''),
	}

	# 提取关键数据
	if 'liquidus_data' in result_dict:
		ld = result_dict['liquidus_data']
		simplified['liquidus_data'] = {
			k: ld[k] for k in ['liquidus_K', 'solidus_K', 'range_K',
			                     'liquidus_C', 'solidus_C', 'range_C']
			if k in ld
		}

	if 'solubility_data' in result_dict:
		sd = result_dict['solubility_data']
		simplified['solubility_data'] = {
			k: sd[k] for k in ['temperature_K', 'solubility', 'matrix_phase']
			if k in sd
		}

	if 'properties_data' in result_dict:
		pd_data = result_dict['properties_data']
		simplified['properties_summary'] = {
			'n_points': len(pd_data.get('x_values', [])),
			'temperature': pd_data.get('temperature'),
		}

	prompt = f"""计算结果:
{json.dumps(simplified, ensure_ascii=False, indent=2)}

计算上下文:
- 组分: {context.get('components', '未知')}
- 计算类型: {context.get('calc_type', '未知')}
- 使用模型: {context.get('model', 'RKM')}

请分析此计算结果。"""
	return prompt


def parse_interpret_response(raw_text: str) -> Optional[str]:
	"""解析结果解读响应（纯文本，无需JSON解析）"""
	if not raw_text or len(raw_text.strip()) < 10:
		return None
	# 清理多余的引号和代码块标记
	text = raw_text.strip()
	text = re.sub(r'^```\w*\n?', '', text)
	text = re.sub(r'\n?```$', '', text)
	return text.strip() if text.strip() else None


# ============================================================================
# 5. 模型推荐 (Model Recommendation)
# ============================================================================

MODEL_RECOMMEND_SYSTEM = """你是一位计算热力学（CALPHAD）方法的专家。给定合金组分和计算类型，
推荐最佳的热力学外推模型。

可选模型：
- RKM: Redlich-Kister-Muggianu（默认模型，适用于大多数二元/三元体系）
- Muggianu: Muggianu对称模型（三元+体系，各组分相互作用对称）
- Toop: Toop非对称模型（三元+体系，一个组分化学性质显著不同）
- UEM1: 统一外推方法（多元体系，复杂相互作用，精度更高但计算较慢）

选择依据：
- 二元体系：通常用 RKM 即可
- 三元对称体系（如 Al-Cu-Zn, Fe-Co-Ni）：Muggianu
- 三元非对称体系（如 Al-Cr-Ni, Fe-C-Mn 中C为非金属）：Toop
- 四元及以上、或需要高精度：UEM1

返回格式（仅JSON，无其他文字）：
{
  "model": "RKM" 或 "Muggianu" 或 "Toop" 或 "UEM1",
  "reasoning": "简要说明(中文)",
  "confidence": 0.0到1.0
}"""


def build_model_recommend_prompt(components: List[str], calc_type: str) -> str:
	"""构建模型推荐用户提示词"""
	calc_type_map = {
		"liquidus_solidus": "液相线/固相线计算",
		"solubility": "溶解度计算",
		"pseudo_binary": "伪二元相图计算",
		"ternary": "三元等温截面计算",
		"surface": "液相面投影计算",
		"properties": "热力学性质计算",
	}
	calc_desc = calc_type_map.get(calc_type, calc_type)

	return f"""合金组分: {', '.join(components)}
计算类型: {calc_desc}
组分数量: {len([c for c in components if c.upper() not in ('VA', 'VACUUM')])}

请推荐最佳外推模型。仅返回JSON。"""


def parse_model_recommendation(raw_text: str) -> Optional[Dict[str, str]]:
	"""解析模型推荐响应，验证模型名合法"""
	try:
		json_match = re.search(r'\{.*\}', raw_text, re.DOTALL)
		if not json_match:
			return None

		data = json.loads(json_match.group())

		model = data.get('model', '').strip()
		valid_models = {'RKM', 'Muggianu', 'Toop', 'UEM1'}

		# 大小写不敏感匹配
		model_map = {m.upper(): m for m in valid_models}
		matched_model = model_map.get(model.upper())

		if not matched_model:
			return None

		return {
			'model': matched_model,
			'reasoning': str(data.get('reasoning', '')),
			'confidence': min(1.0, max(0.0, float(data.get('confidence', 0.5)))),
		}
	except (json.JSONDecodeError, ValueError, Exception):
		return None


# ============================================================================
# 6. 通用对话 (General Chat)
# ============================================================================

GENERAL_CHAT_SYSTEM = """你是一位合金热力学计算助手，集成在 UEM-PyCalphad 合金相图计算工具中。

你的能力：
- 解答计算热力学（CALPHAD方法）相关问题
- 解释合金相图特征
- 帮助用户理解计算结果
- 建议计算参数和模型选择
- 解释错误信息并给出解决建议

请用中文回答。保持简洁专业。"""


def build_chat_prompt(user_message: str, context: Optional[Dict] = None) -> str:
	"""构建通用对话用户提示词"""
	prompt = user_message

	if context:
		ctx_parts = []
		if context.get('components'):
			ctx_parts.append(f"当前体系组分: {', '.join(context['components'])}")
		if context.get('loaded_db'):
			ctx_parts.append(f"已加载数据库: {context['loaded_db']}")
		if context.get('available_phases'):
			ctx_parts.append(f"可用相数量: {len(context['available_phases'])}")
		if context.get('last_calc_type'):
			ctx_parts.append(f"上次计算类型: {context['last_calc_type']}")

		if ctx_parts:
			prompt = f"""[当前上下文]
{chr(10).join(ctx_parts)}

[用户提问]
{user_message}"""

	return prompt
