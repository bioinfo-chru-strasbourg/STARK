#!/usr/bin/env python

import json
import os
import re
import subprocess
from typing import Optional

from config import QUEUES_FILE, ts, ts_savelist, ts_slots, ts_socket_env


def load_queues() -> dict:
    default = {
        "stark": {
            "savelist": ts_savelist,
            "slots": int(ts_slots),
            "description": "Default STARK analysis queue",
        },
        "light": {
            "savelist": "/ts-tmp-light",
            "slots": 2,
            "socket": "/tmp/ts-light.socket",
            "description": "Light STARK analysis queue",
        },
        "medium": {
            "savelist": "/ts-tmp-medium",
            "slots": 4,
            "socket": "/tmp/ts-medium.socket",
            "description": "Medium STARK analysis queue",
        },
        "gpu": {
            "savelist": "/ts-tmp-gpu",
            "slots": 4,
            "socket": "/tmp/ts-gpu.socket",
            "description": "GPU STARK analysis queue",
        },
    }
    if not os.path.exists(QUEUES_FILE):
        os.makedirs(os.path.dirname(QUEUES_FILE), exist_ok=True)
        with open(QUEUES_FILE, "w") as f:
            json.dump(default, f, indent=2)
    try:
        with open(QUEUES_FILE, "r") as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return default


def _resolve_queue_socket(name: str, cfg: dict, is_default: bool) -> str:
    """Return the TS_SOCKET path for this queue, or '' to not override TS_SOCKET."""
    if is_default:
        return ""
    explicit = cfg.get("socket", "")
    if explicit:
        return explicit
    safe_name = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    return f"/tmp/ts-{safe_name}.socket"


def get_queue_env(queue_name: Optional[str] = None) -> str:
    """Return the env prefix (TS_SOCKET + TS_SAVELIST + TS_SLOTS) for the given queue."""
    queues = load_queues()
    default_name = next(iter(queues))
    name = queue_name or default_name
    is_default = name == default_name
    if name in queues:
        cfg = queues[name]
    else:
        cfg = {"savelist": f"/ts-tmp-{name}", "slots": 1}
        is_default = False
    socket = _resolve_queue_socket(name, cfg, is_default)
    socket_part = f"TS_SOCKET={socket} " if socket else ""
    return f"{socket_part}TS_SAVELIST={cfg['savelist']} TS_SLOTS={cfg['slots']} "


def _prepare_queue_for_submission(queue_name: Optional[str] = None) -> tuple:
    """Return (ts_env_prefix, max_slots) and enforce the slot count on the ts daemon."""
    queues = load_queues()
    default_name = next(iter(queues))
    name = queue_name or default_name
    is_default = name == default_name
    if name in queues:
        cfg = queues[name]
    else:
        available = list(queues.keys())
        raise ValueError(
            f"Queue '{name}' is not defined. Available queues: {available}"
        )
    socket = _resolve_queue_socket(name, cfg, is_default)
    socket_part = f"TS_SOCKET={socket} " if socket else ""
    _ts_env = f"{socket_part}TS_SAVELIST={cfg['savelist']} TS_SLOTS={cfg['slots']} "
    slots = int(cfg.get("slots", 1))
    if ts and slots >= 1:
        subprocess.run(
            f"{_ts_env} {ts} -S {slots}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    return _ts_env, slots


def _queue_daemon_is_active(name: str, cfg: dict, is_default: bool) -> bool:
    """Return True if the ts daemon for this queue is running."""
    socket = _resolve_queue_socket(name, cfg, is_default)
    if not socket:
        if ts_socket_env:
            return os.path.exists(ts_socket_env)
        return True
    return os.path.exists(socket)
