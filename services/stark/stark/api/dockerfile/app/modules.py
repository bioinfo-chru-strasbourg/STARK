#!/usr/bin/env python

import json
import os

from config import MODULES_FILE, docker_stark
from security import MANDATORY_KEYS_PARAMS, MANDATORY_PARAMS
from functions import is_empty

def load_modules() -> dict:
    """Load the modules configuration from MODULES_FILE.

    Returns a dict keyed by module name (upper-cased). Each entry contains:
      - image (str): Docker image to use — server-side only, never overridable by clients
      - _description (str): Human-readable description
      - docker_extra_params (str): Extra Docker flags added by the module
      - use_stark_container_mount (bool): Whether to mount the STARK volumes
      - defaults (dict): Default values for queue, threads, memory, prioritize

    If the file does not exist it is created with the built-in STARK default.
    If the file is malformed the built-in default is returned.
    """

    default = {
        "STARK": {
            "_description": "Default STARK analysis module",
            "_enable": True,
            "_available": True,
            "image": docker_stark,
            "docker_extra_params": "",
            "use_stark_container_mount": True,
            "defaults": {
                "queue": "stark",
                "threads": None,
                "memory": None,
                "prioritize": False,
            },
        }
    }

    if not os.path.exists(MODULES_FILE):
        os.makedirs(os.path.dirname(MODULES_FILE), exist_ok=True)
        with open(MODULES_FILE, "w") as f:
            json.dump(default, f, indent=2)

    try:
        with open(MODULES_FILE, "r") as f:
            raw = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return default

    # Normalise: upper-case keys, filter out disabled modules and non-dict entries, but keep the rest of the config as-is for flexibility.
    modules: dict = {}
    for name, cfg in raw.items():
        if isinstance(cfg, dict) and ("_enable" not in cfg or cfg["_enable"]):

            # Key - Module name is upper-cased for case-insensitive matching, but the rest of the config is left as-is to allow flexibility in defining arbitrary extra fields if needed in the future without breaking the API. The mandatory keys for modules are filled with safe defaults if missing or None to avoid security issues. The "image" key is especially important to have a safe default as it defines the Docker image used for analyses and we don't want it to be accidentally misconfigured to an unsafe value.
            key = name.upper()

            # Safe security defaults values for parameters

            # If any of the mandatory keys are present but empty/None, remove them so they will be filled with the safe default values
            for k in MANDATORY_KEYS_PARAMS:
                if k in cfg and is_empty(cfg[k]):
                    del cfg[k]

            # If none of keys param are in module config, add stark image as default
            if not any(k in cfg for k in MANDATORY_KEYS_PARAMS):
                cfg["image"] = docker_stark

            # For each mandatory param, if it's missing, fill it with the default value
            for param, default_val in MANDATORY_PARAMS.items():
                if param not in cfg:
                    cfg[param] = default_val
            modules[key] = cfg

    return modules if modules else default


def resolve_module(json_input: dict) -> dict:
    """Return the module config for the 'module' key in json_input.

    Falls back to the first defined module (STARK by default) only when the key is absent.
    Raises ValueError for unknown module names so the caller can return HTTP 400.
    """

    modules = load_modules()
    first_module = next(iter(modules))

    requested = json_input.get("module")
    if requested is None:
        return modules[first_module]

    key = str(requested).upper()
    if key not in modules:
        raise ValueError(
            f"Unknown module '{requested}'. Available: {', '.join(modules)}"
        )
    return modules[key]


def apply_module_defaults(json_input: dict, module_cfg: dict) -> None:
    """Inject module defaults into json_input for keys not already provided by the client.

    Only fills in queue/threads/memory/prioritize when the client has not
    supplied a non-None/non-empty value for that key.
    Mutates json_input in place.
    """
    defaults = module_cfg.get("defaults", {})
    for key in ("queue", "threads", "memory", "prioritize"):
        default_val = defaults.get(key)
        if default_val is None:
            continue
        client_val = json_input.get(key)
        if client_val is None or (isinstance(client_val, str) and not client_val.strip()):
            json_input[key] = default_val
