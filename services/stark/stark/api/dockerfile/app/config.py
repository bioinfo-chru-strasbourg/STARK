#!/usr/bin/env python

import os

# --- JWT / Auth ---
SECRET_KEY = "a_very_secret_key_that_should_be_in_a_config_file"
ALGORITHM = "HS256"
USERS_FILE = "config/users.json"
QUEUES_FILE = "config/queues.json"
MODULES_FILE = "config/modules.json"
STARK_API_KEY = os.environ.get("STARK_API_KEY", "a_default_super_secret_api_key")

# --- Cluster / Peers ---
# Path to the peers config file (auto-created empty on first start if absent)
PEERS_FILE = "config/peers.json"
# Optional: explicit URL this node is reachable at (e.g. "http://192.168.1.10:8001")
# If not set, the node tries to self-identify by scanning peers via /whoami.
# If self-URL cannot be determined, the node always runs analyses locally.
STARK_API_SELF_URL = os.environ.get("STARK_API_SELF_URL", "")
STARK_API_SELF_PORT = int(os.environ.get("STARK_API_SELF_PORT", 8000))
local_port = STARK_API_SELF_PORT if STARK_API_SELF_URL else 8000

# --- Task Spooler ---
ts = os.environ.get("TS", "")
ts_savelist = os.environ.get("TS_SAVELIST", "/ts-tmp")
ts_slots = int(os.environ.get("TS_SLOTS", os.cpu_count()))
ts_socket_env = os.environ.get("TS_SOCKET", "")

# --- Shell ---
shell = os.environ.get("SHELL", "/bin/ash")

# --- Docker STARK image and parameters ---
docker_stark = os.environ.get("DOCKER_STARK_IMAGE", "stark")
docker_stark_container_mount = os.environ.get(
    "DOCKER_STARK_SERVICE_STARK_API_CONTAINER_MOUNT", ""
)
docker_stark_api_log_folder = os.environ.get(
    "DOCKER_STARK_SERVICE_STARK_API_LOG_FOLDER", "/STARK/services/stark/stark/api"
)
docker_stark_api_runs_folder = os.environ.get(
    "DOCKER_STARK_SERVICE_STARK_API_RUNS_FOLDER", "/STARK/input/runs"
)

# --- UI ---
refresh_interval_ms = int(os.environ.get("STARK_API_REFRESH_INTERVAL", "10")) * 1000

# --- Launch tab access control ---
# STARK_API_LAUNCH_ENABLED: "true" (default) shows the Launch tab to admins;
#   "false" hides it entirely for all users.
launch_enabled: bool = (
    os.environ.get("STARK_API_LAUNCH_ENABLED", "true").strip().lower() == "true"
)

# STARK_API_LAUNCH_MODES: comma-separated list of allowed sub-tabs.
#   Recognised values: run, docker, advanced  (all enabled by default)
_launch_modes_raw = os.environ.get("STARK_API_LAUNCH_MODES", "run,docker,advanced")
launch_modes: list = [
    m.strip()
    for m in _launch_modes_raw.split(",")
    if m.strip() in ("run", "docker", "advanced")
]
# Fallback: if the env var was set but produced no valid entries, keep all
if not launch_modes:
    launch_modes = ["run", "docker", "advanced"]
