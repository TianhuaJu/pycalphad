#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
LLM 后端服务模块
=================

提供 LLM 提供商抽象和领域特定服务：
- 支持 OpenAI 兼容 API（OpenAI / Ollama / 其他兼容端点）
- 支持 Anthropic Messages API
- 支持自定义扩展端点
- 优雅降级：LLM 不可用时所有方法返回 None

用法：
    config = LLMConfig(provider='ollama', model='llama3')
    service = LLMService(config, logger=print)
    phases = service.predict_phases(['AL', 'CU'], all_phases, 'liquidus_solidus')
"""

import json
import threading
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple, Callable, Any

try:
	import httpx
	HTTPX_AVAILABLE = True
except ImportError:
	HTTPX_AVAILABLE = False
	# 尝试自动安装
	try:
		import subprocess, sys
		subprocess.check_call(
			[sys.executable, '-m', 'pip', 'install', 'httpx'],
			stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
			timeout=60
		)
		import httpx
		HTTPX_AVAILABLE = True
	except Exception:
		pass

from llm_prompts import (
	PHASE_PREDICTION_SYSTEM, build_phase_prediction_prompt, parse_phase_prediction,
	COMPOSITION_PARSE_SYSTEM, build_composition_prompt, parse_composition_response,
	PARAMETER_SUGGEST_SYSTEM, build_parameter_prompt, parse_parameter_response,
	RESULT_INTERPRET_SYSTEM, build_result_interpret_prompt, parse_interpret_response,
	MODEL_RECOMMEND_SYSTEM, build_model_recommend_prompt, parse_model_recommendation,
	GENERAL_CHAT_SYSTEM, build_chat_prompt,
)


# ============================================================================
# 配置数据类
# ============================================================================

@dataclass
class LLMConfig:
	"""LLM 配置，可序列化到 JSON 配置文件"""
	provider: str = 'ollama'       # 'openai' | 'anthropic' | 'ollama' | 'grok' | 'custom'
	api_key: str = ''
	base_url: str = 'http://localhost:11434/v1'
	model: str = 'llama3'
	temperature: float = 0.3
	max_tokens: int = 2048
	timeout: int = 30
	enabled: bool = False          # 默认关闭，需用户手动启用

	def to_dict(self) -> Dict:
		"""序列化为字典（用于 JSON 存储）"""
		return asdict(self)

	@classmethod
	def from_dict(cls, data: Dict) -> 'LLMConfig':
		"""从字典反序列化"""
		if not data:
			return cls()
		# 只取合法字段
		valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
		filtered = {k: v for k, v in data.items() if k in valid_fields}
		return cls(**filtered)


# 预置模型配置
PROVIDER_PRESETS = {
	'openai': {
		'base_url': 'https://api.openai.com/v1',
		'model': 'gpt-4o-mini',
		'models': ['gpt-4o-mini', 'gpt-4o', 'gpt-4.1-mini', 'gpt-4.1-nano'],
		'requires_key': True,
	},
	'anthropic': {
		'base_url': 'https://api.anthropic.com/v1',
		'model': 'claude-sonnet-4-20250514',
		'models': ['claude-sonnet-4-20250514', 'claude-haiku-4-20250514'],
		'requires_key': True,
	},
	'ollama': {
		'base_url': 'http://localhost:11434/v1',
		'model': 'llama3',
		'models': ['llama3', 'llama3.1', 'qwen2.5', 'gemma2', 'mistral', 'deepseek-r1'],
		'requires_key': False,
	},
	'grok': {
		'base_url': 'https://api.x.ai/v1',
		'model': 'grok-3-mini',
		'models': ['grok-3-mini', 'grok-3', 'grok-2'],
		'requires_key': True,
	},
	'deepseek': {
		'base_url': 'https://api.deepseek.com/v1',
		'model': 'deepseek-chat',
		'models': ['deepseek-chat', 'deepseek-reasoner'],
		'requires_key': True,
	},
	'siliconflow': {
		'base_url': 'https://api.siliconflow.cn/v1',
		'model': 'Qwen/Qwen2.5-7B-Instruct',
		'models': [],
		'requires_key': True,
	},
	'custom': {
		'base_url': '',
		'model': '',
		'models': [],
		'requires_key': True,
	},
}


def fetch_ollama_models(base_url: str = 'http://localhost:11434/v1',
                        timeout: int = 5) -> Tuple[List[str], str]:
	"""查询 Ollama 服务已安装的模型列表

	Parameters
	----------
	base_url : str
		Ollama 的 base URL（可含 /v1 后缀）
	timeout : int
		连接超时秒数

	Returns
	-------
	(models, message) : 模型名列表 + 状态信息
		失败时 models 为空列表
	"""
	if not HTTPX_AVAILABLE:
		return [], "httpx 未安装"

	# 去掉 /v1 后缀获取原生 Ollama 地址
	host = base_url.rstrip('/')
	if host.endswith('/v1'):
		host = host[:-3]

	try:
		with httpx.Client(timeout=timeout) as client:
			resp = client.get(f"{host}/api/tags")
			resp.raise_for_status()

		data = resp.json()
		models = []
		for m in data.get('models', []):
			name = m.get('name', '')
			# Ollama 返回 "qwen2.5:latest" 格式，取冒号前的短名
			short = name.split(':')[0] if name else ''
			if short and short not in models:
				models.append(short)
		return sorted(models), f"已发现 {len(models)} 个本地模型"

	except httpx.ConnectError:
		return [], f"无法连接 Ollama ({host})，请确认 ollama serve 已启动"
	except Exception as e:
		return [], f"查询失败: {e}"


def fetch_openai_models(base_url: str, api_key: str = '',
                        timeout: int = 10) -> Tuple[List[str], str]:
	"""查询 OpenAI 兼容端点的可用模型列表

	通过标准 GET /models 接口获取，适用于 OpenAI / DeepSeek /
	硅基流动 / 胜算云等所有 OpenAI 兼容 API。

	Parameters
	----------
	base_url : str
		API Base URL（如 https://api.deepseek.com/v1）
	api_key : str
		API Key（Bearer token）
	timeout : int
		连接超时秒数

	Returns
	-------
	(models, message) : 模型名列表 + 状态信息
	"""
	if not HTTPX_AVAILABLE:
		return [], "httpx 未安装"

	url = f"{base_url.rstrip('/')}/models"
	headers = {"Content-Type": "application/json"}
	if api_key:
		headers["Authorization"] = f"Bearer {api_key}"

	try:
		with httpx.Client(timeout=timeout) as client:
			resp = client.get(url, headers=headers)
			resp.raise_for_status()

		data = resp.json()
		models = []
		for m in data.get('data', []):
			mid = m.get('id', '')
			if mid and mid not in models:
				models.append(mid)
		return sorted(models), f"已发现 {len(models)} 个可用模型"

	except httpx.ConnectError:
		return [], f"无法连接 {base_url}，请检查地址"
	except httpx.HTTPStatusError as e:
		code = e.response.status_code
		if code == 401:
			return [], "API Key 无效，请检查密钥"
		return [], f"查询失败: HTTP {code}"
	except Exception as e:
		return [], f"查询失败: {e}"


# ============================================================================
# 抽象提供商基类
# ============================================================================

class LLMProvider(ABC):
	"""LLM 提供商抽象基类"""

	def __init__(self, config: LLMConfig):
		self.config = config

	@abstractmethod
	def chat(self, messages: List[Dict[str, str]], **kwargs) -> str:
		"""发送聊天请求，返回助手回复文本"""
		...

	@abstractmethod
	def test_connection(self) -> Tuple[bool, str]:
		"""测试连接，返回 (是否成功, 信息)"""
		...


# ============================================================================
# OpenAI 兼容提供商（覆盖 OpenAI / Ollama / 其他）
# ============================================================================

class OpenAICompatibleProvider(LLMProvider):
	"""
	OpenAI Chat Completions API 兼容提供商。
	适用于 OpenAI、Ollama、vLLM、LM Studio 等所有兼容端点。
	"""

	def chat(self, messages: List[Dict[str, str]], **kwargs) -> str:
		if not HTTPX_AVAILABLE:
			raise RuntimeError("httpx 未安装，请运行: pip install httpx")

		headers = {"Content-Type": "application/json"}
		if self.config.api_key:
			headers["Authorization"] = f"Bearer {self.config.api_key}"

		payload = {
			"model": self.config.model,
			"messages": messages,
			"temperature": kwargs.get('temperature', self.config.temperature),
			"max_tokens": kwargs.get('max_tokens', self.config.max_tokens),
		}

		url = f"{self.config.base_url.rstrip('/')}/chat/completions"

		with httpx.Client(timeout=self.config.timeout) as client:
			response = client.post(url, json=payload, headers=headers)
			response.raise_for_status()

		data = response.json()
		return data["choices"][0]["message"]["content"]

	def test_connection(self) -> Tuple[bool, str]:
		try:
			# 对 Ollama 先做健康检查（使用原生 API）
			if self.config.provider == 'ollama':
				return self._test_ollama()
			result = self.chat([{"role": "user", "content": "你好，请回复OK"}],
			                   max_tokens=20)
			return True, f"连接成功: {result[:50]}"
		except httpx.ConnectError:
			return False, f"连接失败: 无法连接到 {self.config.base_url}，请检查服务是否启动"
		except httpx.HTTPStatusError as e:
			code = e.response.status_code
			if code == 404:
				return False, f"连接失败: 模型 '{self.config.model}' 未找到或端点不存在"
			elif code == 401:
				return False, "连接失败: API Key 无效"
			return False, f"连接失败: HTTP {code}"
		except Exception as e:
			return False, f"连接失败: {e}"

	def _test_ollama(self) -> Tuple[bool, str]:
		"""Ollama 专用测试：先检查服务，再检查模型"""
		# 从 base_url 推导 Ollama 原生地址（去掉 /v1 后缀）
		base = self.config.base_url.rstrip('/')
		if base.endswith('/v1'):
			ollama_host = base[:-3]
		else:
			ollama_host = base

		try:
			with httpx.Client(timeout=10) as client:
				# 1. 检查 Ollama 服务是否运行
				try:
					resp = client.get(f"{ollama_host}/api/tags")
					resp.raise_for_status()
				except httpx.ConnectError:
					return False, (f"Ollama 服务未运行！\n"
					               f"请启动: ollama serve\n"
					               f"地址: {ollama_host}")
				except Exception:
					# 非标准 Ollama，回退到通用测试
					result = self.chat(
						[{"role": "user", "content": "你好，请回复OK"}],
						max_tokens=20)
					return True, f"连接成功: {result[:50]}"

				# 2. 检查模型是否已安装
				data = resp.json()
				installed = [m.get('name', '').split(':')[0]
				             for m in data.get('models', [])]
				model_name = self.config.model.split(':')[0]
				if installed and model_name not in installed:
					available = ', '.join(installed[:5])
					return False, (f"模型 '{self.config.model}' 未安装！\n"
					               f"请运行: ollama pull {self.config.model}\n"
					               f"已安装: {available}")

				# 3. 真正测试对话
				result = self.chat(
					[{"role": "user", "content": "你好，请回复OK"}],
					max_tokens=20)
				return True, f"连接成功: {result[:50]}"

		except httpx.ConnectError:
			return False, (f"Ollama 服务未运行！\n"
			               f"请启动: ollama serve\n"
			               f"地址: {ollama_host}")
		except httpx.HTTPStatusError as e:
			if e.response.status_code == 404:
				return False, (f"模型 '{self.config.model}' 未找到\n"
				               f"请运行: ollama pull {self.config.model}")
			return False, f"连接失败: HTTP {e.response.status_code}"
		except Exception as e:
			return False, f"连接失败: {e}"


# ============================================================================
# Anthropic 提供商
# ============================================================================

class AnthropicProvider(LLMProvider):
	"""Anthropic Messages API 提供商"""

	def chat(self, messages: List[Dict[str, str]], **kwargs) -> str:
		if not HTTPX_AVAILABLE:
			raise RuntimeError("httpx 未安装，请运行: pip install httpx")

		# 分离 system 消息
		system_content = ""
		user_messages = []
		for msg in messages:
			if msg["role"] == "system":
				system_content += msg["content"] + "\n"
			else:
				user_messages.append(msg)

		# 确保至少有一条用户消息
		if not user_messages:
			user_messages = [{"role": "user", "content": "请回复"}]

		headers = {
			"Content-Type": "application/json",
			"x-api-key": self.config.api_key,
			"anthropic-version": "2023-06-01",
		}

		payload = {
			"model": self.config.model,
			"max_tokens": kwargs.get('max_tokens', self.config.max_tokens),
			"messages": user_messages,
		}

		if system_content.strip():
			payload["system"] = system_content.strip()

		url = f"{self.config.base_url.rstrip('/')}/messages"

		with httpx.Client(timeout=self.config.timeout) as client:
			response = client.post(url, json=payload, headers=headers)
			response.raise_for_status()

		data = response.json()
		# Anthropic 返回 content 数组
		content_blocks = data.get("content", [])
		text_parts = [b["text"] for b in content_blocks if b.get("type") == "text"]
		return "\n".join(text_parts)

	def test_connection(self) -> Tuple[bool, str]:
		try:
			result = self.chat([{"role": "user", "content": "你好，请回复OK"}],
			                   max_tokens=20)
			return True, f"连接成功: {result[:50]}"
		except Exception as e:
			return False, f"连接失败: {e}"


# ============================================================================
# 自定义提供商（默认按 OpenAI 格式，可扩展）
# ============================================================================

class CustomProvider(OpenAICompatibleProvider):
	"""
	自定义提供商，继承 OpenAI 兼容格式。
	用户可以通过指定 base_url 连接任意兼容端点。
	"""
	pass


# ============================================================================
# 提供商工厂
# ============================================================================

def create_provider(config: LLMConfig) -> LLMProvider:
	"""根据配置创建合适的提供商实例"""
	provider_map = {
		'openai': OpenAICompatibleProvider,
		'ollama': OpenAICompatibleProvider,  # Ollama 兼容 OpenAI API
		'grok': OpenAICompatibleProvider,  # xAI Grok 兼容 OpenAI API
		'deepseek': OpenAICompatibleProvider,  # DeepSeek 兼容 OpenAI API
		'siliconflow': OpenAICompatibleProvider,  # 硅基流动兼容 OpenAI API
		'anthropic': AnthropicProvider,
		'custom': CustomProvider,
	}

	provider_cls = provider_map.get(config.provider, CustomProvider)
	return provider_cls(config)


# ============================================================================
# LLM 服务层 — 唯一对外接口
# ============================================================================

class LLMService:
	"""
	LLM 领域服务层。

	提供 5 个领域特定方法 + 1 个通用对话方法。
	所有方法在 LLM 不可用时返回 None（优雅降级）。
	线程安全。
	"""

	def __init__(self, config: Optional[LLMConfig] = None,
	             logger: Optional[Callable] = None):
		self.config = config or LLMConfig()
		self.logger = logger or (lambda msg: None)
		self._provider: Optional[LLMProvider] = None
		self._lock = threading.Lock()

		if self.config.enabled:
			self._init_provider()

	def _init_provider(self):
		"""初始化提供商（内部方法）"""
		try:
			self._provider = create_provider(self.config)
		except Exception as e:
			self.logger(f"LLM 提供商初始化失败: {e}")
			self._provider = None

	def update_config(self, new_config: LLMConfig):
		"""更新配置并重新初始化提供商"""
		with self._lock:
			self.config = new_config
			self._provider = None
			if new_config.enabled:
				self._init_provider()

	def is_available(self) -> bool:
		"""检查 LLM 是否可用"""
		return (self.config.enabled and
		        self._provider is not None and
		        HTTPX_AVAILABLE)

	def test_connection(self) -> Tuple[bool, str]:
		"""测试 LLM 连接"""
		if not HTTPX_AVAILABLE:
			return False, "httpx 未安装，请运行: pip install httpx"
		if not self.config.enabled:
			return False, "LLM 未启用"
		if self._provider is None:
			self._init_provider()
			if self._provider is None:
				return False, "提供商初始化失败"
		return self._provider.test_connection()

	def _call_llm(self, system_prompt: str, user_prompt: str,
	              **kwargs) -> Optional[str]:
		"""通用 LLM 调用包装（内部方法）"""
		if not self.is_available():
			return None

		try:
			messages = [
				{"role": "system", "content": system_prompt},
				{"role": "user", "content": user_prompt},
			]
			with self._lock:
				return self._provider.chat(messages, **kwargs)
		except Exception as e:
			self.logger(f"LLM 调用失败: {e}")
			return None

	# ------------------------------------------------------------------
	# 领域特定方法
	# ------------------------------------------------------------------

	def predict_phases(self, components: List[str],
	                   all_phases: List[str],
	                   calc_type: str = "general") -> Optional[List[str]]:
		"""
		预测计算所需的相列表。

		Returns: 推荐相列表，或 None（LLM 不可用/解析失败）
		"""
		prompt = build_phase_prediction_prompt(components, all_phases, calc_type)
		raw = self._call_llm(PHASE_PREDICTION_SYSTEM, prompt)
		if raw is None:
			return None

		result = parse_phase_prediction(raw, all_phases)
		if result:
			self.logger(f"AI 推荐 {len(result)} 个相: {', '.join(result)}")
		return result

	def parse_composition(self, user_input: str,
	                      available_comps: List[str]) -> Optional[Dict[str, float]]:
		"""
		用自然语言解析合金成分。

		Returns: 成分字典 {'AL': 1.0, 'CU': 0.5, ...}，或 None
		"""
		prompt = build_composition_prompt(user_input, available_comps)
		raw = self._call_llm(COMPOSITION_PARSE_SYSTEM, prompt)
		if raw is None:
			return None

		result = parse_composition_response(raw, available_comps)
		if result:
			comp_str = ', '.join(f"{k}:{v}" for k, v in result.items())
			self.logger(f"AI 解析合金成分: {user_input} → {comp_str}")
		return result

	def suggest_parameters(self, components: List[str],
	                       calc_type: str,
	                       current_params: Optional[Dict] = None) -> Optional[Dict]:
		"""
		建议计算参数。

		Returns: 参数字典 {temp_min, temp_max, temp_step, scan_points, reasoning}
		"""
		prompt = build_parameter_prompt(components, calc_type, current_params)
		raw = self._call_llm(PARAMETER_SUGGEST_SYSTEM, prompt)
		if raw is None:
			return None

		result = parse_parameter_response(raw)
		if result and 'reasoning' in result:
			self.logger(f"AI 参数建议: {result.get('reasoning', '')}")
		return result

	def interpret_result(self, result_dict: Dict,
	                     context: Dict) -> Optional[str]:
		"""
		解读计算结果。

		Returns: 中文分析文本，或 None
		"""
		prompt = build_result_interpret_prompt(result_dict, context)
		raw = self._call_llm(RESULT_INTERPRET_SYSTEM, prompt)
		if raw is None:
			return None

		return parse_interpret_response(raw)

	def recommend_model(self, components: List[str],
	                    calc_type: str) -> Optional[Dict[str, str]]:
		"""
		推荐外推模型。

		Returns: {'model': 'Toop', 'reasoning': '...', 'confidence': 0.8}
		"""
		prompt = build_model_recommend_prompt(components, calc_type)
		raw = self._call_llm(MODEL_RECOMMEND_SYSTEM, prompt)
		if raw is None:
			return None

		result = parse_model_recommendation(raw)
		if result:
			self.logger(f"AI 推荐模型: {result['model']} ({result.get('reasoning', '')})")
		return result

	def chat(self, user_message: str,
	         context: Optional[Dict] = None) -> Optional[str]:
		"""
		通用对话。

		Returns: 助手回复文本，或 None
		"""
		prompt = build_chat_prompt(user_message, context)
		return self._call_llm(GENERAL_CHAT_SYSTEM, prompt)
