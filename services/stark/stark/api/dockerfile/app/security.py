#!/usr/bin/env python

import re
import shlex
from typing import Optional, Dict

# --- Blocklist: destructive shell command patterns ---
_DANGEROUS_COMMAND_PATTERNS = [
    re.compile(
        r"\brm\s+.*(-[a-zA-Z]*[rR][a-zA-Z]*\s.*-[a-zA-Z]*[fF]|-[a-zA-Z]*[fF][a-zA-Z]*\s.*-[a-zA-Z]*[rR]|--recursive|--force).*[/]",
        re.IGNORECASE,
    ),
    re.compile(r"\brm\s+-[a-zA-Z]*[rRfF]{2,}", re.IGNORECASE),
    re.compile(r"\bdd\b.*\bof=/dev/[a-z]", re.IGNORECASE),
    re.compile(r"\bmkfs\b", re.IGNORECASE),
    re.compile(r"\bwipefs\b", re.IGNORECASE),
    re.compile(r"\bshred\b", re.IGNORECASE),
    re.compile(r":\s*\(\s*\)\s*\{.*;\s*\}.*;"),
    re.compile(r">\s*/dev/[sh]d[a-z]"),
    re.compile(r"\bchmod\s+.*[0-7]*7[0-7][0-7]\s+/", re.IGNORECASE),
]

# --- Blocklist: dangerous docker extra parameter patterns ---
_DANGEROUS_DOCKER_PARAMS_PATTERNS = [
    re.compile(r"--privileged", re.IGNORECASE),
    re.compile(
        r"--cap-add\s+(ALL|SYS_ADMIN|SYS_PTRACE|SYS_MODULE|SYS_RAWIO|SYS_BOOT|NET_ADMIN)",
        re.IGNORECASE,
    ),
    re.compile(r"--pid[=\s]+host", re.IGNORECASE),
    re.compile(r"--userns[=\s]+host", re.IGNORECASE),
    re.compile(r"--security-opt\s+seccomp=unconfined", re.IGNORECASE),
    re.compile(r"--device\s+/dev/[sh]d[a-z]", re.IGNORECASE),
    re.compile(r"-v\s*/\s*:", re.IGNORECASE),
    re.compile(r"--volume\s*/\s*:", re.IGNORECASE),
    re.compile(r"-v\s*/(?:etc|root|proc|sys|dev)\b", re.IGNORECASE),
    re.compile(r"--volume\s*/(?:etc|root|proc|sys|dev)\b", re.IGNORECASE),
]

# Valid docker image name: registry/name:tag — no shell metacharacters
_VALID_IMAGE_RE = re.compile(r"^[a-zA-Z0-9][a-zA-Z0-9_.:\-/@]*$")

# Valid Docker container name or ID: alphanumeric start, then alphanumeric / _ . -
_VALID_CONTAINER_NAME_RE = re.compile(r"^[a-zA-Z0-9][a-zA-Z0-9_.\-]*$")

# Forbidden docker flags that can lead to privilege escalation or security issues
_FORBIDDEN_DOCKER_FLAGS = {
    "--name": True,  # Requires a value, but we block any use of --name to prevent container name conflicts
    "--rm": False,  # We enforce --rm for cleanup, so we block any use of --rm to prevent conflicts
    "-d": False,
    "--detach": False,
    "--restart": True,  # Requires a value, but we block any use of --restart to prevent unintended container persistence
    "--init": False,
    "--hostname": True,  # Requires a value, but we block any use of --hostname to prevent spoofing
    "-ti": False,
    "-t": False,
    "--tty": False,
    "-i": False,
    "--interactive": False,
}

# Keys that are mandatory in the module config for security reasons and will be filled with safe defaults if missing/None in the config file. The presence of the mandatory keys is crucial for ensuring that each module has a defined Docker image, service name, and container name, which are essential for the secure operation of the analysis execution. The mandatory params will be filled with safe defaults to ensure secure operation even if the module is misconfigured or missing these fields.
MANDATORY_KEYS_PARAMS = ("module", "image", "service", "container", "endpoint")
MANDATORY_PARAMS = {
    "docker_extra_params": "",
    "use_stark_container_mount": False,
    "command_prefix": "",
    "command_postfix": "",
    "api_key_variable": None,
    "bearer_token_variable": None,
}
# Sensitive keys include both the mandatory keys that must be present in the module config for security reasons, and the mandatory params that will be filled with safe defaults if missing to ensure secure operation even with misconfigured modules.
SENSITIVE_KEYS = MANDATORY_KEYS_PARAMS + tuple(MANDATORY_PARAMS.keys())

# Keys allowed from client input for analyses, to prevent abuse via extra fields in the JSON input that are not expected by the API and could be used to inject malicious content if not properly handled. This is a whitelist of known safe fields that the API accepts from the client for launching analyses. Any field outside of this set will be ignored if sent by the client.
ALLOWED_CLIENT_KEYS = {
    "command",
    "analysis_name",
    "run",
    "queue",
    "threads",
    "memory",
    "prioritize",
}

def _safe_split(s: str) -> list:
    """Split user input into tokens safely (no eval, no shell)."""
    return shlex.split(s) if s else []


def _validate_command(command: str) -> None:
    """Raise ValueError if the command matches a known destructive pattern."""
    for pattern in _DANGEROUS_COMMAND_PATTERNS:
        if pattern.search(command):
            raise ValueError(
                f"Command rejected: matches a dangerous pattern ({pattern.pattern!r})"
            )


def _validate_docker_extra_params(params: str) -> None:
    """Raise ValueError if docker_extra_params contains privilege-escalation flags."""
    for pattern in _DANGEROUS_DOCKER_PARAMS_PATTERNS:
        if pattern.search(params):
            raise ValueError(
                f"docker_extra_params rejected: dangerous flag detected ({pattern.pattern!r})"
            )


def _validate_image(image: str) -> None:
    """Raise ValueError if the Docker image name contains shell metacharacters."""
    if not _VALID_IMAGE_RE.match(image):
        raise ValueError(
            f"Invalid image name: {image!r}. Only alphanumeric, '_', '.', ':', '-', '/' characters are allowed."
        )


def _validate_container_name(name: str) -> None:
    """Raise ValueError if the container name or ID contains shell metacharacters."""
    if not _VALID_CONTAINER_NAME_RE.match(name):
        raise ValueError(
            f"Invalid container name: {name!r}. "
            "Only alphanumeric, '_', '.', '-' characters are allowed, starting with alphanumeric."
        )


def _sanitize_analysis_name(name: str) -> str:
    """Sanitize an analysis name: keep only [A-Za-z0-9._-], truncate to 80 chars."""
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    return sanitized[:80] if sanitized else "UNKNOWN"


def _sanitize_docker_extra_params(
    extra_params: str, extra_forbidden: Optional[Dict[str, bool]] = None
) -> str:
    """
    Removes forbidden Docker flags from the parameter string.
    Correctly handles flags with values (e.g., `--name value` or `--name=value`).
    Allows extending the list of forbidden flags dynamically.
    """

    if not extra_params:
        return extra_params

    # Base forbidden list (existing)
    forbidden = dict(_FORBIDDEN_DOCKER_FLAGS)

    # Merge optional forbidden flags
    if extra_forbidden:
        forbidden.update(extra_forbidden)

    tokens = shlex.split(extra_params)
    cleaned = []
    i = 0

    while i < len(tokens):
        token = tokens[i]

        # Check for "--flag=value" form
        matched_assignment = False
        for flag, takes_value in forbidden.items():
            if token.startswith(flag + "="):
                matched_assignment = True
                break

        if matched_assignment:
            # Skip "--flag=value"
            i += 1
            continue

        # Check for flag as separate token
        if token in forbidden:
            takes_value = forbidden[token]
            i += 1

            # Skip next token if this flag takes a value
            if takes_value and i < len(tokens):
                i += 1
            continue

        # Allowed token -> keep it
        cleaned.append(token)
        i += 1

    return " ".join(cleaned)
