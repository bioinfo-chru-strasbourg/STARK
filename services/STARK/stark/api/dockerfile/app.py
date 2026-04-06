#!/usr/bin/env python

from fastapi import (
    FastAPI,
    Request,
    Depends,
    HTTPException,
    status,
    Query,
    Form,
    Security,
)
from fastapi.responses import HTMLResponse, PlainTextResponse, JSONResponse
from fastapi.templating import Jinja2Templates
from fastapi.security import (
    OAuth2PasswordBearer,
    OAuth2PasswordRequestForm,
    APIKeyHeader,
)
from fastapi.staticfiles import StaticFiles
from jose import JWTError, jwt
from passlib.context import CryptContext
from pydantic import BaseModel
import os
import subprocess
import json
import random
import string
import time
import datetime
import re
import concurrent.futures
from typing import Optional, Union, List


# Command by JSON
# This will create an analysis with name "run_test" and run the command "RUN_TEST" in STARK Docker (default queue).
# Example: {"run": "RUN_TEST", "analysis_name": "run_test"}
# This will create an analysis with name "test_cmd" and run "echo hello world" directly in the shell (no Docker isolation).
# Example: {"command": "echo hello world", "analysis_name": "test_cmd"}
# This will run "echo hello from docker" inside a new container based on "alpine" (with predefined mounts).
# Example: {"command_docker": "sleep 10 && echo hello from docker", "image": "alpine", "analysis_name": "test_cmd_docker", "queue": "other"}
#
# All three modes accept an optional "queue" key to target a specific task-spooler queue
# defined in config/queues.json. When omitted, the first (default) queue is used.
# Example: {"run": "MY_RUN", "queue": "light"}
# Example: {"command": "python3 /scripts/heavy.py", "queue": "gpu", "analysis_name": "my_job"}
#
# config/queues.json is auto-created on first start from TS_SAVELIST/TS_SLOTS env vars.
# To configure multiple queues, edit it manually, e.g.:
# {
#   "stark": {"savelist": "/ts-tmp",      "slots": 1, "description": "Default STARK queue"},
#   "light": {"savelist": "/ts-tmp-light","slots": 4, "description": "Lightweight tasks"},
#   "gpu":   {"savelist": "/ts-tmp-gpu",  "slots": 2, "description": "GPU analyses"}
# }


# --- Configuration ---
SECRET_KEY = "a_very_secret_key_that_should_be_in_a_config_file"
ALGORITHM = "HS256"
USERS_FILE = "config/users.json"
QUEUES_FILE = "config/queues.json"
STARK_API_KEY = os.environ.get("STARK_API_KEY", "a_default_super_secret_api_key")

pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")

app = FastAPI()

app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates")


# --- Models ---
class Token(BaseModel):
    access_token: str
    token_type: str


class User(BaseModel):
    username: str
    groups: List[str] = []


# --- Authentication ---

api_key_header = APIKeyHeader(name="X-API-Key", auto_error=False)

oauth2_scheme = OAuth2PasswordBearer(tokenUrl="token", auto_error=False)


# --- STARK API LOGIC ---

# Task Spooler
# If TS is set, it will be used to schedule tasks with 'ts -L <analysisIDNAME> <command>'.
ts = os.environ.get("TS", "")
ts_savelist = os.environ.get("TS_SAVELIST", "/ts-tmp")
ts_slots = os.environ.get("TS_SLOTS", "1")
ts_socket_env = os.environ.get(
    "TS_SOCKET", ""
)  # container's default ts socket (used by the default queue)
# ts_env removed — use get_queue_env() which resolves per-task queue from config/queues.json

# Shell
# The shell used to run commands. Default to /bin/ash for better compatibility with minimal Docker images, but can be overridden by setting the SHELL environment variable (e.g. to /bin/bash if available and preferred).
shell = os.environ.get("SHELL", "/bin/ash")

# Docker STARK image and parameters
# The Docker image to use for analyses can be set via the DOCKER_STARK_IMAGE environment variable (default: "stark").
# The DOCKER_STARK_SERVICE_STARK_API_CONTAINER_MOUNT environment variable can be used to set predefined mount parameters for the Docker run command (e.g. "-v /data:/data"). This allows analyses to access specific host directories in a controlled way, without giving full access to the API container's filesystem. By default, no extra mounts are added.
# The DOCKER_STARK_SERVICE_STARK_API_LOG_FOLDER environment variable specifies where the API should store analysis JSON, info and output files. This folder should be mounted in the Docker image used for analyses to allow them to write their output and status. By default, it is set to "/STARK/services/stark/stark/api", which is a folder inside the API container that should be mounted as a volume in the Docker image used for analyses (e.g. with "-v /path/to/logs:/STARK/services/stark/stark/api").
# The DOCKER_STARK_SERVICE_STARK_API_RUNS_FOLDER environment variable specifies an optional folder where input runs can be stored and accessed by analyses. If an analysis is launched with a "run" parameter that matches a file or folder in this directory, it will be used as the run input and its hash will be computed for the analysis IDNAME. By default, it is set to "/STARK/input/runs", which can be mounted as a volume in the API container and used to provide input data for analyses.
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
refresh_interval_ms = int(os.environ.get("STARK_API_REFRESH_INTERVAL", "10")) * 1000


# --- Security: command and docker parameter validation ---

# Note: these regex patterns are best-effort and may not catch all possible destructive commands or parameters, but they provide a basic layer of protection against common dangerous patterns. Always review and test the patterns to ensure they fit your specific use case and threat model.
# The command validation focuses on common destructive patterns like 'rm -rf /', 'dd if=/dev/zero of=/dev/sda', 'mkfs', 'shred', fork bombs, and redirections to block devices. The docker parameter validation focuses on flags that grant elevated privileges or access to the host system, such as '--privileged', '--cap-add=ALL', '--pid=host', and volume mounts of sensitive paths. The image name validation ensures that only simple, well-formed image names are accepted, without shell metacharacters that could lead to command injection.

# Blocklist for destructive shell command patterns (case-insensitive, best-effort)
_DANGEROUS_COMMAND_PATTERNS = [
    re.compile(
        r"\brm\s+.*(-[a-zA-Z]*[rR][a-zA-Z]*\s.*-[a-zA-Z]*[fF]|-[a-zA-Z]*[fF][a-zA-Z]*\s.*-[a-zA-Z]*[rR]|--recursive|--force).*[/]",
        re.IGNORECASE,
    ),
    re.compile(
        r"\brm\s+-[a-zA-Z]*[rRfF]{2,}", re.IGNORECASE
    ),  # rm -rf, rm -fr, rm -Rf ...
    re.compile(r"\bdd\b.*\bof=/dev/[a-z]", re.IGNORECASE),  # dd to block device
    re.compile(r"\bmkfs\b", re.IGNORECASE),  # filesystem format
    re.compile(r"\bwipefs\b", re.IGNORECASE),  # wipe filesystem signatures
    re.compile(r"\bshred\b", re.IGNORECASE),  # secure file deletion
    re.compile(r":\s*\(\s*\)\s*\{.*;\s*\}.*;"),  # fork bomb
    re.compile(r">\s*/dev/[sh]d[a-z]"),  # redirect to block device
    re.compile(
        r"\bchmod\s+.*[0-7]*7[0-7][0-7]\s+/", re.IGNORECASE
    ),  # chmod 777+ on root
]

# Blocklist for dangerous docker extra parameters
_DANGEROUS_DOCKER_PARAMS_PATTERNS = [
    re.compile(r"--privileged", re.IGNORECASE),
    re.compile(
        r"--cap-add\s+(ALL|SYS_ADMIN|SYS_PTRACE|SYS_MODULE|SYS_RAWIO|SYS_BOOT|NET_ADMIN)",
        re.IGNORECASE,
    ),
    re.compile(r"--pid[=\s]+host", re.IGNORECASE),
    re.compile(r"--userns[=\s]+host", re.IGNORECASE),
    re.compile(r"--security-opt\s+seccomp=unconfined", re.IGNORECASE),
    re.compile(r"--device\s+/dev/[sh]d[a-z]", re.IGNORECASE),  # raw block device
    re.compile(r"-v\s*/\s*:", re.IGNORECASE),  # mount host root
    re.compile(r"--volume\s*/\s*:", re.IGNORECASE),
    re.compile(
        r"-v\s*/(?:etc|root|proc|sys|dev)\b", re.IGNORECASE
    ),  # sensitive host paths
    re.compile(r"--volume\s*/(?:etc|root|proc|sys|dev)\b", re.IGNORECASE),
]

# Valid docker image name: registry/name:tag — no shell metacharacters
_VALID_IMAGE_RE = re.compile(r"^[a-zA-Z0-9][a-zA-Z0-9_.:\-/@]*$")


# --- Functions ---


# This function checks the API key provided in the "X-API-Key" header against the expected value. If the key is valid, it returns the key; otherwise, it returns None, allowing the request to proceed to user authentication if no valid API key is provided. This way, both service calls with a valid API key and user-authenticated calls can access the API, while invalid API keys are rejected.
# Note: In a production environment, you would typically want to implement a more robust API key management system, possibly with multiple keys, expiration, and revocation capabilities. This example uses a single static key for simplicity.
def get_api_key(api_key: str = Security(api_key_header)):
    if api_key == STARK_API_KEY:
        return api_key
    else:
        return None

# This function loads users from a JSON file, supporting both a new format (where each user is an object with "password" and "groups" fields) and a legacy format (where each user is just a plain password string). If the users file does not exist or is invalid, it creates it with a default user ("stark" with password "password" and admin group) and a guest user. This allows for flexible user management while maintaining backward compatibility with older configurations.
# The expected format of the users JSON file is:
# {
#     "stark": {"password": "secret", "groups": ["admin", "users"]},
#     "john":  {"password": "pass",   "groups": ["users"]}
# }
def load_users():
    """Loads users from the JSON file.
    
    Expected format:
    {
        "stark": {"password": "secret", "groups": ["admin", "users"]},
        "john":  {"password": "pass",   "groups": ["users"]}
    }
    Legacy format (plain password string) is also supported.
    """
    # default = {"stark": {"password": "password", "groups": ["admin"]}}
    default = {
        "stark": {"password": "password", "groups": ["admin"]},
        "guest": {"password": "guest", "groups": []},
    }
    if not os.path.exists(USERS_FILE):
        os.makedirs(os.path.dirname(USERS_FILE), exist_ok=True)
        with open(USERS_FILE, "w") as f:
            json.dump(default, f, indent=2)
    try:
        with open(USERS_FILE, "r") as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return default


def load_queues() -> dict:
    """Load queue definitions from config/queues.json.

    Each entry: {"savelist": "/ts-tmp", "slots": 1, "description": "..."}
    The first entry is the default queue used when no 'queue' key is provided.
    Auto-created from TS_SAVELIST/TS_SLOTS env vars if the file is missing,
    which preserves backward-compatible behaviour.

    The optional 'socket' key sets the TS_SOCKET path for the daemon of that queue,
    which is what makes queues truly independent. Without it, non-default queues
    get an auto-derived socket at /tmp/ts-<name>.socket.
    The default (first) queue never overrides TS_SOCKET — it uses the container's value.

    Example config/queues.json:
    {
        "stark": {"savelist": "/ts-tmp",     "slots": 1, "description": "Default STARK queue"},
        "light": {"savelist": "/ts-tmp-light","slots": 4, "socket": "/tmp/ts-light.socket", "description": "Light STARK queue"},
        "medium": {"savelist": "/ts-tmp-medium","slots": 2, "socket": "/tmp/ts-medium.socket", "description": "Medium STARK queue"},
        "gpu":   {"savelist": "/ts-tmp-gpu",  "slots": 1, "socket": "/tmp/ts-gpu.socket", "description": "GPU STARK queue"}
    }
    """
    default = {
        "stark": {
            "savelist": ts_savelist,
            "slots": int(ts_slots),
            "description": "Default STARK analysis queue",
        },
        "light": {
            "savelist": "/ts-tmp-light",
            "slots": 4,
            "socket": "/tmp/ts-light.socket",
            "description": "Light STARK analysis queue",
        },
        "medium": {
            "savelist": "/ts-tmp-medium",
            "slots": 2,
            "socket": "/tmp/ts-medium.socket",
            "description": "Medium STARK analysis queue",
        },
        "gpu": {
            "savelist": "/ts-tmp-gpu",
            "slots": 1,
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
    """Return the TS_SOCKET path for this queue, or '' to not override TS_SOCKET.

    - Default queue: never override TS_SOCKET; let ts use whatever TS_SOCKET is
      already set in the container environment (backward compatible).
    - Non-default queues: use 'socket' from config if present, otherwise
      auto-derive /tmp/ts-<name>.socket to ensure daemon isolation.
    """
    if is_default:
        return ""  # keep container's TS_SOCKET unchanged
    explicit = cfg.get("socket", "")
    if explicit:
        return explicit
    safe_name = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    return f"/tmp/ts-{safe_name}.socket"


def get_queue_env(queue_name: Optional[str] = None) -> str:
    """Return the env prefix (TS_SOCKET + TS_SAVELIST + TS_SLOTS) for the given queue.

    Each non-default queue gets a unique TS_SOCKET so its daemon is fully isolated.
    The default queue never overrides TS_SOCKET — it uses the container's env value.
    Use _prepare_queue_for_submission() when submitting tasks (enforces slot count).
    """
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


def _prepare_queue_for_submission(queue_name: Optional[str] = None) -> str:
    """Return the env prefix for the queue AND enforce the slot count on the ts daemon.

    task-spooler identifies daemons by TS_SOCKET. Each non-default queue gets its
    own TS_SOCKET path so its daemon is fully isolated from other queues.
    Because TS_SLOTS is only read at daemon startup, we call 'ts -S <slots>' on
    every submission to enforce the configured concurrency at all times.
    """
    queues = load_queues()
    default_name = next(iter(queues))
    name = queue_name or default_name
    is_default = name == default_name
    if name in queues:
        cfg = queues[name]
    else:
        available = list(queues.keys())
        raise ValueError(
            f"Queue '{name}' is not defined. " f"Available queues: {available}"
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
    return _ts_env


def _queue_daemon_is_active(name: str, cfg: dict, is_default: bool) -> bool:
    """Return True if the ts daemon for this queue is running.

    The daemon is identifed by its TS_SOCKET file. For the default queue (which
    does not override TS_SOCKET), we check ts_socket_env if known, or assume
    active (the container's ts daemon should always be running).
    For non-default queues, we check the auto-derived or configured socket file.
    """
    socket = _resolve_queue_socket(name, cfg, is_default)
    if not socket:
        # Default queue: check TS_SOCKET from env if available, otherwise assume active
        if ts_socket_env:
            return os.path.exists(ts_socket_env)
        return True
    return os.path.exists(socket)


# This function authenticates a user by checking the provided username and password against the loaded users. It supports both the new format (where each user is an object with "password" and "groups") and the legacy format (where each user is just a plain password string). If authentication is successful, it returns a User object with the username and groups; otherwise, it returns False. This allows for flexible authentication while maintaining compatibility with older user configurations.
def authenticate_user(username: str, password: str):
    """Authenticates a user. Supports both legacy (plain string) and new (dict) formats."""
    users_db = load_users()
    if username not in users_db:
        return False
    entry = users_db[username]
    # New format: {"password": "...", "groups": [...]}
    if isinstance(entry, dict):
        if entry.get("password") == password:
            return User(username=username, groups=entry.get("groups", []))
    # Legacy format: plain password string
    elif entry == password:
        return User(username=username, groups=[])
    return False

# This function creates a JWT access token containing the provided data (e.g., username and groups) and signs it with the secret key. The token can then be used for authenticated requests to the API. In a production environment, you would typically want to include an expiration time in the token and possibly other claims, but this example focuses on the basic functionality of encoding user information in a JWT.
def create_access_token(data: dict):
    to_encode = data.copy()
    encoded_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)
    return encoded_jwt

# This function retrieves the current user based on the provided JWT token in the "Authorization" header. It decodes the token, extracts the username, and checks it against the loaded users. If the token is valid and the user exists, it returns a User object with the username and groups; otherwise, it raises an HTTP 401 Unauthorized exception. This function is used as a dependency in routes that require user authentication.
async def get_current_user(token: str = Depends(oauth2_scheme)):
    if not token:
        return None
    credentials_exception = HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Could not validate credentials",
        headers={"WWW-Authenticate": "Bearer"},
    )
    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
        username: str = payload.get("sub")
        if username is None:
            raise credentials_exception
    except JWTError:
        raise credentials_exception

    users_db = load_users()
    if username not in users_db:
        raise credentials_exception
    entry = users_db[username]
    groups = entry.get("groups", []) if isinstance(entry, dict) else []
    return User(username=username, groups=groups)

# This function checks for a valid API key in the "X-API-Key" header and returns "service" if the key is valid, allowing service calls to access the API without user authentication. If no valid API key is provided, it falls back to user authentication using the JWT token. If the user is authenticated successfully, it returns the User object; otherwise, it raises an HTTP 401 Unauthorized exception. This allows both service calls with a valid API key and user-authenticated calls to access the API, while rejecting invalid API keys.
async def get_current_user_or_service(
    api_key: str = Security(api_key_header), user: User = Depends(get_current_user)
):
    if api_key and api_key == STARK_API_KEY:
        return "service"
    if user:
        return user
    raise HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED, detail="Not authenticated"
    )


# This function generates a random string of letters and digits of a specified length (default is 6). It is used to create unique identifiers for analyses and runs, ensuring that each analysis has a distinct IDNAME even if the same analysis name is used multiple times. The generated string consists of uppercase and lowercase letters and digits, providing a large pool of possible combinations to minimize the risk of collisions.
def randomStringDigits(stringLength=6):
    """Generate a random string of letters and digits"""
    lettersAndDigits = string.ascii_letters + string.digits
    return "".join(random.choice(lettersAndDigits) for i in range(stringLength))


# This function attempts to extract the analysisIDNAME from the output of the 'ts -i' command, which provides information about a scheduled task. It looks for patterns that indicate the analysisIDNAME, either from the Docker container name (if the task is a Docker-based analysis) or from the output file path (if the task is a plain-command analysis). If it finds a match, it returns the extracted analysisIDNAME; otherwise, it returns None. This function is useful for correlating tasks in the task spooler with their corresponding analyses based on their IDNAME.
def _extract_analysis_name(ts_info_stdout: str) -> Optional[str]:
    """Extract analysisIDNAME from 'ts -i' output.
    Works for docker-based tasks (--name flag) and plain-command tasks (output file path).
    """
    # Docker-based tasks: docker run --name <analysisIDNAME>
    m = re.search(r"--name\s+(\S+)", ts_info_stdout)
    if m:
        return m.group(1)
    # Plain-command tasks: the output redirect contains the analysisIDNAME
    m = re.search(
        r"(STARK\.[A-Za-z0-9]+\.ID-[A-Za-z0-9]+-NAME-\S+?)\.output", ts_info_stdout
    )
    if m:
        return m.group(1)
    return None


# This function validates a shell command by checking it against a list of known destructive patterns. If the command matches any of the patterns, it raises a ValueError with a message indicating that the command is rejected and specifying which pattern was matched. This provides a basic layer of protection against potentially harmful commands being executed through the API, while still allowing a wide range of safe commands to be used.
# Note: The patterns used for validation are best-effort and may not catch all possible destructive commands, so it's important to review and test the patterns to ensure they fit your specific use case and threat model. Additionally, consider implementing further security measures such as running commands in a restricted environment or using containerization to mitigate risks.
def _validate_command(command: str) -> None:
    """Raise ValueError if the command matches a known destructive pattern."""
    for pattern in _DANGEROUS_COMMAND_PATTERNS:
        if pattern.search(command):
            raise ValueError(
                f"Command rejected: matches a dangerous pattern ({pattern.pattern!r})"
            )


# This function validates the extra parameters provided for the Docker run command by checking them against a list of known dangerous patterns that could grant elevated privileges or access to the host system. If any of the dangerous patterns are detected in the provided parameters, it raises a ValueError with a message indicating that the parameters are rejected and specifying which pattern was matched. This helps to ensure that only safe and controlled Docker parameters are allowed when launching analyses through the API.
# Note: The patterns used for validation are best-effort and may not catch all possible dangerous parameters, so it's important to review and test the patterns to ensure they fit your specific use case and threat model. Additionally, consider implementing further security measures such as running Docker containers with a restricted set of capabilities or using container security profiles to mitigate risks.
def _validate_docker_extra_params(params: str) -> None:
    """Raise ValueError if docker_extra_params contains privilege-escalation flags."""
    for pattern in _DANGEROUS_DOCKER_PARAMS_PATTERNS:
        if pattern.search(params):
            raise ValueError(
                f"docker_extra_params rejected: dangerous flag detected ({pattern.pattern!r})"
            )


# This function validates a Docker image name by checking it against a regex pattern that allows only alphanumeric characters, underscores, dots, colons, dashes, slashes, and at signs. If the image name contains any characters that do not match this pattern, it raises a ValueError with a message indicating that the image name is invalid. This helps to prevent command injection vulnerabilities by ensuring that only well-formed image names are accepted when launching analyses through the API.
# Note: The regex pattern used for validation is a basic check and may not cover all valid Docker image name formats, but it provides a reasonable level of protection against common injection patterns. Always review and test the pattern to ensure it fits your specific use case and threat model. Additionally, consider implementing further security measures such as running Docker containers with a restricted set of capabilities or using container security profiles to mitigate risks.
def _validate_image(image: str) -> None:
    """Raise ValueError if the Docker image name contains shell metacharacters."""
    if not _VALID_IMAGE_RE.match(image):
        raise ValueError(
            f"Invalid image name: {image!r}. Only alphanumeric, '_', '.', ':', '-', '/' characters are allowed."
        )


# This function sanitizes an analysis name by replacing any characters that are not alphanumeric, underscores, dots, or dashes with underscores. It also truncates the sanitized name to 64 characters to comply with Docker container name limits. If the resulting sanitized name is empty, it returns "UNKNOWN". This ensures that analysis names used in Docker container names are valid and do not contain characters that could cause issues when creating containers.
# Note: The choice of allowed characters and the truncation length are based on Docker's naming rules, but you may want to adjust the sanitization logic to fit your specific requirements or constraints. Always review and test the sanitization function to ensure it behaves as expected in your use case.
def _sanitize_analysis_name(name: str) -> str:
    """Sanitize an analysis name: keep only [A-Za-z0-9._-], replace anything else with '_'.
    Truncated to 64 characters to satisfy Docker container name limits."""
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    return sanitized[:64] if sanitized else "UNKNOWN"


# This function builds the analysisIDNAME and analysesRUNNAME based on the provided JSON input. It generates a unique analysisID using a random string, and determines the runID and runMD5 based on the "run" parameter in the input. If the "run" parameter points to an existing file or directory, it computes a hash of its contents to create a unique identifier. If an "analysis_name" is provided in the input, it uses that as the runID after sanitization. Finally, it constructs the analysisIDNAME in the format "STARK.{analysisID}.ID-{runMD5}-NAME-{runID}" and returns both the analysisIDNAME and analysesRUNNAME. This function is crucial for ensuring that each analysis has a unique and descriptive identifier based on its input parameters, which helps with tracking and managing analyses in the system.
# Note: The hashing of the run input is designed to create a unique identifier for the analysis based on the content of the input, which can help with caching and avoiding duplicate analyses. However, be mindful of the performance implications of hashing large files or directories, and consider implementing additional checks or limits if necessary. Additionally, the sanitization of the analysis name ensures that it can be safely used in Docker container names, but you may want to further customize the naming scheme to fit your specific requirements or constraints. Always review and test the function to ensure it behaves as expected in your use case.
def _build_analysis_idname(json_input: dict) -> tuple:
    """Build analysisIDNAME and analysesRUNNAME from json_input.
    Returns (analysisIDNAME, analysesRUNNAME)."""
    analysesID = randomStringDigits(12)
    runID = "UNKNOWN"
    runMD5 = randomStringDigits(41)

    if "run" in json_input:
        name = json_input["run"].split(":")[0]
        runID = os.path.basename(name)

        runFolder = ""
        if os.path.isdir(name):
            runFolder = name
        elif os.path.isdir(os.path.join(docker_stark_api_runs_folder, name)):
            runFolder = os.path.join(docker_stark_api_runs_folder, name)

        if runFolder:
            myCmd = f"find {runFolder} -maxdepth 1 -type f -print0 | xargs -0 sha1sum | cut -b-40 | sha1sum | awk '{{print $1}}'"
        else:
            myCmd = f"echo {name} | sha1sum | awk '{{print $1}}'"

        runMD5 = (
            subprocess.run(myCmd, shell=True, stdout=subprocess.PIPE)
            .stdout.decode("utf-8")
            .strip()
        )

        if "analysis_name" in json_input:
            runID = _sanitize_analysis_name(str(json_input["analysis_name"]))

    elif "analysis_name" in json_input:
        runID = _sanitize_analysis_name(str(json_input["analysis_name"]))
        runMD5 = randomStringDigits(41)

    analysesRUNNAME = runID
    analysesNAME = f"ID-{runMD5}-NAME-{runID}"
    analysisIDNAME = f"STARK.{analysesID}.{analysesNAME}"
    return analysisIDNAME, analysesRUNNAME


# This function builds and queues a STARK Docker analysis based on the provided JSON input. It first converts the input to a JSON string, then generates a unique analysisIDNAME and analysesRUNNAME using the _build_analysis_idname function. It prepares the Docker run parameters, including mounts and container name, and creates file paths for storing the analysis input, output, and status information. The function writes the JSON input to a file that will be used by the Docker container. It then constructs a shell command to run the Docker container with the specified parameters, redirecting output to a file and writing status information upon completion. Finally, it executes the command using subprocess and returns the analysisIDNAME if successful, or raises a RuntimeError if no task ID is returned by task-spooler. This function is central to launching analyses through the API in a controlled and trackable manner.
# Note: The function assumes that the Docker image specified in the configuration is set up to read the analysis input from the provided JSON file and write its output and status to the specified locations. Ensure that the Docker image and the API are configured correctly to allow for this interaction. Additionally, consider implementing error handling and logging to capture any issues that may arise during the execution of the Docker container or the task scheduling process. Always review and test the function to ensure it behaves as expected in your use case.
def queue_analysis(json_input: dict) -> str:
    """Build and queue a STARK Docker analysis from a JSON input dict.
    Returns the analysis IDNAME string on success, or raises RuntimeError."""
    json_dump = json.dumps(json_input)

    analysisIDNAME, analysesRUNNAME = _build_analysis_idname(json_input)

    docker_name = f" --name {analysisIDNAME} "
    docker_parameters = f" --rm {docker_stark_container_mount} {docker_name} "

    analysisFOLDER = docker_stark_api_log_folder
    analysisFILE = os.path.join(analysisFOLDER, f"{analysisIDNAME}.json")
    analysisINFOFILE = os.path.join(analysisFOLDER, f"{analysisIDNAME}.info")
    analysisOUTPUTFILE = os.path.join(analysisFOLDER, f"{analysisIDNAME}.output")

    with open(analysisFILE, "w") as f:
        f.write(json_dump)

    _ts_env = _prepare_queue_for_submission(json_input.get("queue"))
    ts_cmd = f"{_ts_env}{ts} -L {analysisIDNAME}" if ts else ""
    myCmd = (
        f'{ts_cmd} sh -c "'
        f"docker run {docker_parameters} {docker_stark} "
        f"--analysis_name={analysesRUNNAME} --analysis={analysisFILE} "
        f"> {analysisOUTPUTFILE} 2>&1 "
        f"&& (echo 'done' > {analysisINFOFILE} && exit 0) "
        f"|| (echo 'failed' > {analysisINFOFILE} && exit 1)\""
    )

    task_id = (
        subprocess.run(myCmd, shell=True, stdout=subprocess.PIPE)
        .stdout.decode("utf-8")
        .strip()
    )
    if not task_id:
        raise RuntimeError("task-spooler returned no task ID")
    return analysisIDNAME


# This function queues a raw shell command for execution, based on the "command" key in the provided JSON input. It first validates the command against known destructive patterns to ensure it is safe to execute. It then generates a unique analysisIDNAME using the _build_analysis_idname function, and prepares file paths for storing the command input, output, and status information. The function writes the JSON input to a file that can be used by the command if needed. It constructs a shell command to execute the provided command, redirecting output to a file and writing status information upon completion. Finally, it executes the command using subprocess and returns the analysisIDNAME if successful, or raises a RuntimeError if no task ID is returned by task-spooler. This function allows for flexible execution of arbitrary shell commands while maintaining a level of safety through validation and controlled execution.
# Note: Executing raw shell commands can be dangerous, so it's crucial to ensure that the validation function is robust and that the API is properly secured to prevent unauthorized access. Additionally, consider implementing further security measures such as running commands in a restricted environment or using containerization to mitigate risks. Always review and test the function to ensure it behaves as expected in your use case.
def queue_command(json_input: dict) -> str:
    """Queue a raw shell command (json key 'command').
    Returns the analysisIDNAME string on success, or raises RuntimeError."""
    command = json_input["command"]
    _validate_command(command)
    analysisIDNAME, _ = _build_analysis_idname(json_input)

    analysisFOLDER = docker_stark_api_log_folder
    analysisFILE = os.path.join(analysisFOLDER, f"{analysisIDNAME}.json")
    analysisINFOFILE = os.path.join(analysisFOLDER, f"{analysisIDNAME}.info")
    analysisOUTPUTFILE = os.path.join(analysisFOLDER, f"{analysisIDNAME}.output")

    with open(analysisFILE, "w") as f:
        f.write(json.dumps(json_input))

    _ts_env = _prepare_queue_for_submission(json_input.get("queue"))
    ts_cmd = f"{_ts_env}{ts} -L {analysisIDNAME}" if ts else ""
    myCmd = (
        f'{ts_cmd} sh -c "'
        f"{command} "
        f"> {analysisOUTPUTFILE} 2>&1 "
        f"&& (echo 'done' > {analysisINFOFILE} && exit 0) "
        f"|| (echo 'failed' > {analysisINFOFILE} && exit 1)\""
    )

    task_id = (
        subprocess.run(myCmd, shell=True, stdout=subprocess.PIPE)
        .stdout.decode("utf-8")
        .strip()
    )
    if not task_id:
        raise RuntimeError("task-spooler returned no task ID")
    return analysisIDNAME


# This function queues a Docker command for execution, based on the "command_docker" key in the provided JSON input. It validates the Docker image name and any extra parameters to ensure they are safe to use. It generates a unique analysisIDNAME using the _build_analysis_idname function, and prepares file paths for storing the command input, output, and status information. The function writes the JSON input to a file that can be used by the Docker container if needed. It constructs a shell command to run the specified Docker command with the provided parameters, redirecting output to a file and writing status information upon completion. Finally, it executes the command using subprocess and returns the analysisIDNAME if successful, or raises a RuntimeError if no task ID is returned by task-spooler. This function allows for flexible execution of commands within Docker containers while maintaining a level of safety through validation and controlled execution.
# Note: Running commands in Docker containers can provide an additional layer of isolation, but it's still important to ensure that the validation functions are robust and that the API is properly secured to prevent unauthorized access. Additionally, consider implementing further security measures such as running Docker containers with a restricted set of capabilities or using container security profiles to mitigate risks. Always review and test the function to ensure it behaves as expected in your use case.
def queue_command_docker(json_input: dict) -> str:
    """Queue a docker command with predefined mount parameters (json key 'command_docker').
    JSON keys:
      - command_docker      (required): command to run inside the container
      - image               (required): docker image to use
      - docker_extra_params (optional): free-text extra docker run parameters (e.g. --entrypoint ...)
    Returns the analysisIDNAME string on success, or raises RuntimeError."""
    command_docker = json_input["command_docker"]
    image = json_input.get("image")
    if not image:
        raise ValueError("'image' is required when using 'command_docker'")
    _validate_image(image)
    docker_extra_params = json_input.get("docker_extra_params", "")
    if docker_extra_params:
        _validate_docker_extra_params(docker_extra_params)

    analysisIDNAME, _ = _build_analysis_idname(json_input)

    docker_name = f"--name {analysisIDNAME}"
    docker_parameters = (
        f"--rm {docker_stark_container_mount} {docker_name} {docker_extra_params}"
    )

    analysisFOLDER = docker_stark_api_log_folder
    analysisFILE = os.path.join(analysisFOLDER, f"{analysisIDNAME}.json")
    analysisINFOFILE = os.path.join(analysisFOLDER, f"{analysisIDNAME}.info")
    analysisOUTPUTFILE = os.path.join(analysisFOLDER, f"{analysisIDNAME}.output")

    with open(analysisFILE, "w") as f:
        f.write(json.dumps(json_input))

    _ts_env = _prepare_queue_for_submission(json_input.get("queue"))
    ts_cmd = f"{_ts_env}{ts} -L {analysisIDNAME}" if ts else ""
    myCmd = (
        f'{ts_cmd} sh -c "'
        f"docker run {docker_parameters} {image} {command_docker} "
        f"> {analysisOUTPUTFILE} 2>&1 "
        f"&& (echo 'done' > {analysisINFOFILE} && exit 0) "
        f"|| (echo 'failed' > {analysisINFOFILE} && exit 1)\""
    )

    task_id = (
        subprocess.run(myCmd, shell=True, stdout=subprocess.PIPE)
        .stdout.decode("utf-8")
        .strip()
    )
    if not task_id:
        raise RuntimeError("task-spooler returned no task ID")
    return analysisIDNAME


# --- Routes ---


# The /token route allows users to obtain a JWT access token by providing their username and password. It uses the OAuth2PasswordRequestForm to parse the credentials from the request, authenticates the user, and if successful, creates and returns an access token that can be used for authenticated requests to the API. If authentication fails, it raises an HTTP 401 Unauthorized exception with a message indicating that the username or password is incorrect. This route is essential for enabling user authentication and securing the API endpoints that require a valid token.
# Note: In a production environment, you would typically want to include additional claims in the JWT token (such as an expiration time) and implement a more robust authentication and token management system. Always review and test the authentication flow to ensure it behaves as expected in your use case.
@app.post("/token", response_model=Token)
async def login_for_access_token(form_data: OAuth2PasswordRequestForm = Depends()):
    user = authenticate_user(form_data.username, form_data.password)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect username or password",
            headers={"WWW-Authenticate": "Bearer"},
        )
    access_token = create_access_token(
        data={"sub": user.username, "groups": user.groups}
    )
    return {"access_token": access_token, "token_type": "bearer"}


# The /me route returns the current user's information, including their username and groups. It uses the get_current_user dependency to retrieve the authenticated user based on the provided JWT token. If the user is not authenticated, it raises an HTTP 401 Unauthorized exception with a message indicating that authentication is required. This route allows users to verify their authentication status and view their associated information.
# Note: This route requires a valid JWT token in the "Authorization" header to access the user's information. Always review and test the authentication flow to ensure it behaves as expected in your use case.
@app.get("/me")
async def get_me(user: User = Depends(get_current_user)):
    """Returns the current user's info (username and groups)."""
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED, detail="Not authenticated"
        )
    return {"username": user.username, "groups": user.groups}


# The root route ("/") serves the main HTML page of the application. It uses the Jinja2Templates to render the "index.html" template, passing in the request object, a version string based on the modification time of the "static/script.js" file (to help with cache busting), and a refresh interval in milliseconds. This route provides the user interface for interacting with the API and viewing analysis results.
# Note: Ensure that the "index.html" template and the associated static files are properly set up in the "templates" and "static" directories, respectively. The versioning mechanism for the script.js file helps to ensure that users receive the latest version of the JavaScript code when it changes, which can be important for functionality and security updates. Always review and test the route to ensure it behaves as expected in your use case.
@app.get("/", response_class=HTMLResponse)
async def read_root(request: Request):
    import os

    version = str(int(os.path.getmtime("static/script.js")))
    return templates.TemplateResponse(
        "index.html",
        {
            "request": request,
            "version": version,
            "refresh_interval_ms": refresh_interval_ms,
        },
    )


# The /analysis route allows authorized users to launch a new analysis by providing a JSON payload with the necessary parameters. It checks if the user is authorized (either an admin user or a service call with a valid API key) and then processes the input to queue the analysis using one of the queue functions (queue_command, queue_command_docker, or queue_analysis) based on the keys present in the input. If the analysis is successfully queued, it returns the analysisIDNAME; otherwise, it returns an error message with a 400 status code. This route is central to enabling users to launch analyses through the API while enforcing access control and input validation.
# Note: Ensure that the input JSON is properly structured and that the necessary keys are included for the type of analysis being launched. The access control checks help to ensure that only authorized users can launch analyses, which is important for maintaining the security and integrity of the system. Always review and test the route to ensure it behaves as expected in your use case.
@app.post("/analysis")
async def stark_launch(
    request: Request,
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    # Restrict launch to admin users (service calls are always allowed)
    if authorized != "service":
        if not isinstance(authorized, User) or "admin" not in authorized.groups:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin group required to launch analyses",
            )
    json_input = await request.json()
    try:
        if "command" in json_input:
            analysisIDNAME = queue_command(json_input)
        elif "command_docker" in json_input:
            analysisIDNAME = queue_command_docker(json_input)
        else:
            analysisIDNAME = queue_analysis(json_input)
        return PlainTextResponse(content=analysisIDNAME, status_code=200)
    except Exception as e:
        return PlainTextResponse(content=f"KO: {e}", status_code=400)


# The /list route retrieves the list of all tasks currently in the task spooler queue. It executes the 'ts -l' command to get the list of tasks, parses the output to extract relevant information about each task (such as ID, state, E-Level, times, and run name), and returns this information as a JSON response. If there are no tasks or if an error occurs while executing the command, it returns an empty list. This route allows users to monitor the status of their queued analyses and commands.
# Note: The parsing of the 'ts -l' output is designed to be flexible to accommodate variations in the output format, but it may need adjustments if the output structure changes significantly. Always review and test the route to ensure it behaves as expected in your use case, especially with different versions of task spooler that may have different output formats. Additionally, consider implementing error handling and logging to capture any issues that may arise during the execution of the command or the parsing of its output.
@app.get("/list")
async def list_task():
    """
    List all tasks from all configured queues.

    Usage:

    # GOOD from external docker container
    curl -s -X GET -H 'Content-Type: application/json' http://localhost:4200/list | python3 -m json.tool
    # GOOD from internal docker container
    curl -s -X GET -H 'Content-Type: application/json' http://stark-module-stark-submodule-stark-service-api:8000/list | python3 -m json.tool

    """
    line_regex = re.compile(
        r"^(?P<id>\d+)\s+"
        r"(?P<state>\w+)\s+"
        r"(?P<output>\S+)\s+"
        r"(?P<elevel>\S*)\s*"  # Optional E-Level
        r"(?P<times>[\d./\s-]*)\s+"  # Optional Times
        r"(?P<command>.*)$"
    )
    all_tasks = []
    queues = load_queues()
    default_name = next(iter(queues))
    for queue_name, cfg in queues.items():
        is_default = queue_name == default_name
        if not _queue_daemon_is_active(queue_name, cfg, is_default):
            continue
        socket = _resolve_queue_socket(queue_name, cfg, is_default)
        socket_part = f"TS_SOCKET={socket} " if socket else ""
        _ts_env = f"{socket_part}TS_SAVELIST={cfg['savelist']} TS_SLOTS={cfg['slots']} "
        command = f"{_ts_env} {ts} -l"
        try:
            result = subprocess.run(
                command,
                shell=True,
                executable=shell,
                capture_output=True,
                text=True,
                check=False,
            )
            if result.returncode != 0 and (
                "tasks" in result.stderr.lower() or not result.stdout
            ):
                continue
            lines = result.stdout.strip().split("\n")
            if len(lines) > 1:
                for line in lines[1:]:
                    match = line_regex.match(line)
                    if match:
                        data = match.groupdict()
                        run_name_match = re.search(r"-NAME-([^\s\]]+)", data["command"])
                        run_name = run_name_match.group(1) if run_name_match else "N/A"
                        all_tasks.append(
                            {
                                "id": int(data["id"].strip()),
                                "state": data["state"].strip(),
                                "elevel": (
                                    int(data["elevel"].strip())
                                    if data["elevel"].isdigit()
                                    else None
                                ),
                                "times": (
                                    float(data["times"].strip().split("/")[0])
                                    if data["times"] != ""
                                    else None
                                ),
                                "run_name": run_name,
                                "queue": queue_name,
                            }
                        )
        except FileNotFoundError:
            raise HTTPException(status_code=500, detail=f"Command not found: {ts}")
    return JSONResponse(content=all_tasks)


# The /queue route allows authorized users to either retrieve the current task spooler queue or perform various actions on individual tasks (such as getting info, logs, analysis details, killing, prioritizing, removing, or swapping tasks). The specific action is determined by the "action" query parameter, and the target task is specified by the "id" query parameter when applicable. This route provides a flexible interface for managing tasks in the task spooler while enforcing access control to ensure that only authorized users can perform potentially destructive actions.
# Note: The implementation of the various actions (info, log, analysis, kill, prioritize, remove, swap) is not included in the provided code snippet, so you would need to implement the logic for each action based on the task spooler commands and the structure of your system. Additionally, ensure that the access control checks are properly enforced to prevent unauthorized users from performing destructive actions on tasks. Always review and test the route to ensure it behaves as expected in your use case, especially with regard to the different actions and their effects on the task spooler queue.
@app.get("/queue")
async def queue(
    action: str = Query(
        "list",
        enum=[
            "list",
            "info",
            "log",
            "analysis",
            "kill",
            "prioritize",
            "remove",
            "swap",
        ],
    ),
    id: Optional[str] = Query(None),
    queue: Optional[str] = Query(
        None, description="Queue name (default: first configured queue)"
    ),
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """
    Get the task spooler queue or perform an action on a task.

    The optional 'queue' parameter selects which queue to operate on.
    When omitted, the default (first) queue from config/queues.json is used,
    which preserves backward-compatible behaviour.
    For 'action=list', all queues are aggregated and each task carries a 'queue' field.

    Usage:

    # List all queues
    curl -s -X GET -H 'Content-Type: application/json' -H "X-API-Key: a_default_super_secret_api_key" http://localhost:4200/queue?action=list | python3 -m json.tool
    # Kill a task in a specific queue
    curl -s "http://localhost:4200/queue?action=kill&id=3&queue=light" -H "X-API-Key: a_default_super_secret_api_key"
    """
    _ts_env = get_queue_env(queue)

    if action == "list":
        line_regex = re.compile(
            r"^(?P<id>\d+)\s+"
            r"(?P<state>\w+)\s+"
            r"(?P<output>\S+)\s+"
            r"(?P<elevel>\S*)\s*"  # Optional E-Level
            r"(?P<times>[\d./\s-]*)\s+"  # Optional Times
            r"(?P<command>.*)$"
        )
        all_tasks = []
        queues_all = load_queues()
        default_name_q = next(iter(queues_all))
        queues_to_scan = (
            {queue: queues_all.get(queue, {"savelist": f"/ts-tmp-{queue}", "slots": 1})}
            if queue
            else queues_all
        )
        for queue_name, cfg in queues_to_scan.items():
            is_default_q = queue_name == default_name_q
            if not _queue_daemon_is_active(queue_name, cfg, is_default_q):
                continue
            socket_q = _resolve_queue_socket(queue_name, cfg, is_default_q)
            socket_part_q = f"TS_SOCKET={socket_q} " if socket_q else ""
            _q_env = (
                f"{socket_part_q}TS_SAVELIST={cfg['savelist']} TS_SLOTS={cfg['slots']} "
            )
            command = f"{_q_env} {ts} -l"
            try:
                result = subprocess.run(
                    command,
                    shell=True,
                    executable=shell,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                if result.returncode != 0 and (
                    "tasks" in result.stderr.lower() or not result.stdout
                ):
                    continue
                lines = result.stdout.strip().split("\n")
                if len(lines) > 1:
                    for line in lines[1:]:
                        match = line_regex.match(line)
                        if match:
                            data = match.groupdict()
                            run_name_match = re.search(
                                r"-NAME-([^\s\]]+)", data["command"]
                            )
                            run_name = (
                                run_name_match.group(1) if run_name_match else "N/A"
                            )
                            all_tasks.append(
                                {
                                    "id": data["id"].strip(),
                                    "state": data["state"].strip(),
                                    "output": data["output"].strip(),
                                    "elevel": (
                                        "SUCCESS"
                                        if data["elevel"].strip() == "0"
                                        else (
                                            "FAILED: " + data["elevel"].strip()
                                            if data["elevel"].isdigit()
                                            else ""
                                        )
                                    ),
                                    "times": "",
                                    "run_name": run_name,
                                    "queue": queue_name,
                                    "_ts_env": _q_env,
                                }
                            )
            except FileNotFoundError:
                raise HTTPException(status_code=500, detail=f"Command not found: {ts}")

        # Enrich running and finished tasks with elapsed/finished time (parallel)
        def fetch_time(task):
            state = task["state"].lower()
            if state not in ("running", "finished"):
                return
            try:
                info_result = subprocess.run(
                    f"{task['_ts_env']} {ts} -i {task['id']}",
                    shell=True,
                    executable=shell,
                    capture_output=True,
                    text=True,
                    check=False,
                    timeout=5,
                )
                stdout = info_result.stdout
                if state == "finished":
                    time_match = re.search(r"Time run:\s*([\d.]+)s", stdout)
                    task["times"] = time_match.group(1) if time_match else ""
                else:
                    start_match = re.search(r"Start time:\s*(.+)", stdout)
                    if start_match:
                        start_str = start_match.group(1).strip()
                        try:
                            start_dt = datetime.datetime.strptime(
                                start_str, "%a %b %d %H:%M:%S %Y"
                            )
                            elapsed = time.time() - start_dt.timestamp()
                            task["times"] = f"{elapsed:.1f}"
                        except ValueError:
                            task["times"] = ""
                    else:
                        task["times"] = ""
            except Exception:
                task["times"] = ""

        tasks_to_enrich = [
            t for t in all_tasks if t["state"].lower() in ("running", "finished")
        ]
        if tasks_to_enrich:
            with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
                executor.map(fetch_time, tasks_to_enrich)

        # Strip internal field before returning
        for t in all_tasks:
            t.pop("_ts_env", None)

        return JSONResponse(content=all_tasks)
    else:
        # Restrict destructive actions to admin group
        restricted_actions = {"kill", "prioritize", "remove", "swap"}
        if action in restricted_actions and authorized != "service":
            if not isinstance(authorized, User) or "admin" not in authorized.groups:
                raise HTTPException(
                    status_code=status.HTTP_403_FORBIDDEN,
                    detail="Admin group required for this action",
                )

        # Handling for other actions like info, log, kill
        if not id:
            raise HTTPException(
                status_code=400, detail="Task ID is required for this action"
            )

        action_map = {
            "info": "-i",
            "kill": "-k",
            "prioritize": "-u",
            "remove": "-r",
            "swap": "-U",
        }

        # For kill, also stop the Docker container if it is running
        if action == "kill":
            info_cmd = f"{_ts_env} {ts} -i {id}"
            info_result = subprocess.run(
                info_cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=False,
                timeout=10,
            )
            container_match = re.search(r"--name\s+(\S+)", info_result.stdout)
            if container_match:
                container_name = container_match.group(1)
                subprocess.run(
                    f"docker stop {container_name}",
                    shell=True,
                    capture_output=True,
                    text=True,
                    check=False,
                    timeout=30,
                )

        # For analysis: return the original JSON used to launch the task
        elif action == "analysis":
            info_cmd = f"{_ts_env} {ts} -i {id}"
            info_result = subprocess.run(
                info_cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=False,
                timeout=10,
            )
            analysis_name = _extract_analysis_name(info_result.stdout)
            if not analysis_name:
                return PlainTextResponse(
                    "Could not determine task name from ts info.", status_code=404
                )
            json_file = os.path.join(
                docker_stark_api_log_folder, f"{analysis_name}.json"
            )
            if not os.path.isfile(json_file):
                return PlainTextResponse(
                    f"JSON file not found: {json_file}", status_code=404
                )
            try:
                with open(json_file, "r", errors="replace") as f:
                    content = f.read()
                try:
                    content = json.dumps(
                        json.loads(content), indent=2, ensure_ascii=False
                    )
                except json.JSONDecodeError:
                    pass
                return PlainTextResponse(content)
            except OSError as e:
                return PlainTextResponse(
                    f"Error reading JSON file: {e}", status_code=500
                )

        # For log: return the output file
        if action == "log":
            info_cmd = f"{_ts_env} {ts} -i {id}"
            info_result = subprocess.run(
                info_cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=False,
                timeout=10,
            )
            analysis_name = _extract_analysis_name(info_result.stdout)
            if not analysis_name:
                return PlainTextResponse(
                    "Could not determine task name from ts info.", status_code=404
                )
            output_file = os.path.join(
                docker_stark_api_log_folder, f"{analysis_name}.output"
            )
            if not os.path.isfile(output_file):
                return PlainTextResponse(
                    f"Output file not found: {output_file}", status_code=404
                )
            try:
                with open(output_file, "r", errors="replace") as f:
                    return PlainTextResponse(f.read())
            except OSError as e:
                return PlainTextResponse(
                    f"Error reading output file: {e}", status_code=500
                )

        ts_action = action_map.get(action)
        if not ts_action:
            raise HTTPException(status_code=400, detail=f"Invalid action: {action}")

        command = f"{_ts_env} {ts} {ts_action} {id}"
        try:
            result = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                check=False,
                timeout=10,
            )
            response_text = result.stdout if result.stdout else result.stderr
            return PlainTextResponse(response_text)
        except subprocess.TimeoutExpired:
            return PlainTextResponse(f"Command timed out.", status_code=504)
        except FileNotFoundError:
            return PlainTextResponse(f"Command not found: {ts}", status_code=500)


# The /relaunch/{ts_id} route allows authorized users to re-queue a finished task using its original JSON file. It retrieves the analysis name from the task spooler information for the specified task ID, locates the corresponding JSON file in the log folder, and then uses the queue_analysis function to launch a new analysis with the same parameters. This route is useful for quickly re-running analyses without having to manually prepare the input again, while still enforcing access control to ensure that only authorized users can perform this action.
# Note: Ensure that the task ID provided in the URL corresponds to a finished task and that the original JSON file is still available in the log folder. The access control checks help to ensure that only authorized users can relaunch analyses, which is important for maintaining the security and integrity of the system. Always review and test the route to ensure it behaves as expected in your use case, especially with regard to the retrieval of task information and the handling of the JSON file for relaunching the analysis. Additionally, consider implementing error handling and logging to capture any issues that may arise during the retrieval of task information or the relaunching of the analysis.
@app.post("/relaunch/{ts_id}")
async def relaunch_task(
    ts_id: str,
    queue: Optional[str] = Query(
        None,
        description="Queue name the task belongs to (default: first configured queue)",
    ),
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """Re-queue a finished task using its original JSON file.
    Pass the same 'queue' parameter used when the task was launched."""
    if authorized != "service":
        if not isinstance(authorized, User) or "admin" not in authorized.groups:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin group required to relaunch analyses",
            )

    _ts_env = get_queue_env(queue)

    # Retrieve the analysis name from ts -i
    info_result = subprocess.run(
        f"{_ts_env} {ts} -i {ts_id}",
        shell=True,
        capture_output=True,
        text=True,
        check=False,
        timeout=10,
    )
    analysis_name = _extract_analysis_name(info_result.stdout)
    if not analysis_name:
        raise HTTPException(
            status_code=404, detail="Could not determine task name from ts info."
        )

    json_file = os.path.join(docker_stark_api_log_folder, f"{analysis_name}.json")
    if not os.path.isfile(json_file):
        raise HTTPException(status_code=404, detail=f"JSON file not found: {json_file}")

    try:
        with open(json_file, "r") as f:
            json_input = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        raise HTTPException(status_code=500, detail=f"Error reading JSON file: {e}")

    try:
        if "command" in json_input:
            analysisIDNAME = queue_command(json_input)
        elif "command_docker" in json_input:
            analysisIDNAME = queue_command_docker(json_input)
        else:
            analysisIDNAME = queue_analysis(json_input)
        return PlainTextResponse(content=analysisIDNAME, status_code=200)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Relaunch failed: {e}")


# --- Main ---

if __name__ == '__main__':
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
