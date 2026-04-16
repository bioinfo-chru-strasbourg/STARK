#!/usr/bin/env python

import re

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


def _sanitize_analysis_name(name: str) -> str:
    """Sanitize an analysis name: keep only [A-Za-z0-9._-], truncate to 64 chars."""
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    return sanitized[:64] if sanitized else "UNKNOWN"
