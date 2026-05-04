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


def _validate_docker_command(cmd: str):
    # Interdire redirections & opérateurs shell
    if any(op in cmd for op in [";", "|", "&&", "||", "`", "$(", ">", "<"]):
        raise ValueError("command_docker contains unsafe shell operators")


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


def _sanitize_analysis_name(name: str) -> str:
    """Sanitize an analysis name: keep only [A-Za-z0-9._-], truncate to 64 chars."""
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    return sanitized[:64] if sanitized else "UNKNOWN"


def _sanitize_docker_extra_params(extra_params: str) -> str:
    tokens = shlex.split(extra_params)
    filtered = []

    for t in tokens:
        if any(t == f or t.startswith(f + "=") for f in _FORBIDDEN_DOCKER_FLAGS):
            continue
        filtered.append(t)

    return " ".join(filtered)


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

        # Allowed token → keep it
        cleaned.append(token)
        i += 1

    return " ".join(cleaned)
