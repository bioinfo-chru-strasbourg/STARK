#!/usr/bin/env python
"""Archive-level actions: read log/JSON, relaunch, delete.

All operations work from the ``analysis_id_name`` (e.g.
``STARK.cIXZMaOhfOad.ID-PEdPKFxK1E-NAME-DNA_TEST_small``) which is the
stem of the files written by ``tasks.py`` into ``docker_stark_api_log_folder``.

Security notes
--------------
* ``analysis_id_name`` is validated against a strict allowlist regex before any
  filesystem operation.  Shell metacharacters and path traversal are rejected.
* Delete removes **only** the four known extensions (``.json``,
  ``.json.stark_analysis``, ``.info``, ``.output``) inside
  ``docker_stark_api_log_folder``.  No ``os.remove`` is ever called on a
  caller-supplied path; the path is always constructed server-side.
"""

import json
import os
import re
from typing import Union

from fastapi import APIRouter, Depends, HTTPException, status  # pyright: ignore[reportMissingImports]
from fastapi.responses import JSONResponse, PlainTextResponse  # pyright: ignore[reportMissingImports]

from authentication import get_current_user_or_service
from config import docker_stark_api_log_folder
from models import User
from queues import load_queues
from functions import _safe_mtime

router = APIRouter()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Strict allowlist: STARK.<alphanum>.<ID-…-NAME-…> — no slashes, no dots
# Strict allowlist: STARK.<alphanum>.<trailing-segment> — no slashes or
# shell metacharacters; the trailing segment may contain alphanumerics and
# the characters ``.``, ``_``, ``@``, ``=``, and ``-``.
_ANALYSIS_ID_RE = re.compile(r"^STARK\.[A-Za-z0-9]+\.[A-Za-z0-9_.@=-]+$")


def _validate_analysis_id(analysis_id_name: str) -> None:
    """Raise HTTPException 400 if the id is not in the expected format."""
    if not _ANALYSIS_ID_RE.match(analysis_id_name):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid analysis_id_name format: '{analysis_id_name}'",
        )


def _log_folder_path(analysis_id_name: str, ext: str) -> str:
    """Build an absolute path inside docker_stark_api_log_folder.

    The resulting path is verified to stay within the log folder (defence
    against any residual path-traversal attempt).
    """
    base = os.path.join(docker_stark_api_log_folder, f"{analysis_id_name}{ext}")
    real_base = os.path.realpath(base)
    real_folder = os.path.realpath(docker_stark_api_log_folder)
    if not real_base.startswith(real_folder + os.sep) and real_base != real_folder:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Path traversal detected.",
        )
    return base


@router.get("/archives")
async def list_archives():
    """List all historical tasks based on .json parameter files in the log folder.

    Each entry is built from the `.json` parameter file written at task submission
    time and the companion `.info` file written at task completion.

    Returns up to 1000 most recent entries (by file modification time, newest first).
    """
    import glob as _glob

    pattern = os.path.join(docker_stark_api_log_folder, "STARK.*.json")
    files = sorted(_glob.glob(pattern), key=os.path.getmtime, reverse=True)

    max_entries = 100000

    result = []
    for json_path in files[:max_entries]:
        base = os.path.splitext(os.path.basename(json_path))[0]
        info_path = json_path.replace(".json", ".info")

        # Parse run_name from filename (-NAME-<run_name>)
        m = re.search(r"-NAME-(.+)$", base)
        run_name = m.group(1) if m else base

        # Read original parameters
        try:
            with open(json_path, "r", errors="replace") as f:
                params = json.load(f)
        except Exception:
            params = {}

        # Queue
        queues_all = load_queues()
        default_name_q = next(iter(queues_all))
        queue = params.get("queue", default_name_q)

        # Threads
        threads = params.get("threads")

        # Determine status from .info file
        status = "unknown"
        if os.path.exists(info_path):
            try:
                with open(info_path, "r", errors="replace") as f:
                    info_content = f.read().strip()
                if info_content == "finished":
                    status = "finished"
                elif info_content == "failed":
                    status = "failed"
                else:
                    status = info_content
            except Exception:
                status = "unknown"

        # Times
        # Start time
        mtime = _safe_mtime(json_path)
        # End time
        end_date = _safe_mtime(info_path)
        # Execution time (from file modification times)
        exec_time = (end_date - mtime) if end_date is not None else None

        # Append to results
        result.append(
            {
                "analysis_id_name": base,
                "run_name": run_name,
                "queue": queue,
                "threads": threads,
                "status": status,
                "mtime": mtime,
                "end_date": end_date,
                "exec_time": exec_time,
            }
        )

    return JSONResponse(content=result)


# ---------------------------------------------------------------------------
# Read-only actions (no auth required — mirrors GET /queue behaviour)
# ---------------------------------------------------------------------------

@router.get("/archive/{analysis_id_name}/log")
async def archive_log(analysis_id_name: str):
    """Return the ``.output`` log file for an archived analysis."""
    _validate_analysis_id(analysis_id_name)
    path = _log_folder_path(analysis_id_name, ".output")
    if not os.path.isfile(path):
        return PlainTextResponse(
            f"Output file not found for '{analysis_id_name}'", status_code=404
        )
    try:
        with open(path, "r", errors="replace") as f:
            return PlainTextResponse(f.read())
    except OSError as e:
        return PlainTextResponse(f"Error reading output file: {e}", status_code=500)


@router.get("/archive/{analysis_id_name}/json")
async def archive_json(analysis_id_name: str):
    """Return the ``.json`` parameters file for an archived analysis (pretty-printed)."""
    _validate_analysis_id(analysis_id_name)
    path = _log_folder_path(analysis_id_name, ".json")
    if not os.path.isfile(path):
        return PlainTextResponse(
            f"JSON file not found for '{analysis_id_name}'", status_code=404
        )
    try:
        with open(path, "r", errors="replace") as f:
            content = f.read()
        try:
            content = json.dumps(json.loads(content), indent=2, ensure_ascii=False)
        except json.JSONDecodeError:
            pass
        return PlainTextResponse(content)
    except OSError as e:
        return PlainTextResponse(f"Error reading JSON file: {e}", status_code=500)


# ---------------------------------------------------------------------------
# Admin actions
# ---------------------------------------------------------------------------

def _require_admin(authorized: Union[User, str]) -> None:
    if authorized != "service":
        if not isinstance(authorized, User) or "admin" not in authorized.groups:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin group required",
            )


@router.post("/archive/{analysis_id_name}/relaunch")
async def archive_relaunch(
    analysis_id_name: str,
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """Re-submit an archived analysis from its saved ``.json`` parameters.

    The function detects whether the original task was a *module* analysis
    (presence of a ``command`` key) or a standard STARK/docker analysis, and
    routes to the appropriate queue function — exactly like ``POST /analysis``
    and ``POST /analysis/module`` do.

    No task-spooler ID is involved (the task no longer exists in the spooler).
    """
    _require_admin(authorized)
    _validate_analysis_id(analysis_id_name)

    path = _log_folder_path(analysis_id_name, ".json")
    if not os.path.isfile(path):
        return PlainTextResponse(
            f"JSON file not found for '{analysis_id_name}'", status_code=404
        )

    try:
        with open(path, "r") as f:
            json_input = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        return PlainTextResponse(f"Error reading JSON file: {e}", status_code=500)

    # Route to the correct task function
    from routers.analysis import _run_locally

    try:
        analysis_id_name_new = await _run_locally(json_input)
        return PlainTextResponse(content=analysis_id_name_new, status_code=200)
    except Exception as e:
        return PlainTextResponse(content=f"Relaunch failed: {e}", status_code=400)


@router.delete("/archive/{analysis_id_name}")
async def archive_delete(
    analysis_id_name: str,
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """Delete all files associated with an archived analysis.

    Removes only the four known extensions inside ``docker_stark_api_log_folder``:
    * ``.json``
    * ``.json.stark_analysis``
    * ``.info``
    * ``.output``

    Returns a summary of what was deleted and what was not found.
    """
    _require_admin(authorized)
    _validate_analysis_id(analysis_id_name)

    extensions = [".json", ".json.stark_analysis", ".info", ".output"]
    deleted = []
    not_found = []

    for ext in extensions:
        path = _log_folder_path(analysis_id_name, ext)
        if os.path.isfile(path):
            try:
                os.remove(path)
                deleted.append(ext)
            except OSError as e:
                return PlainTextResponse(
                    f"Error deleting '{ext}': {e}", status_code=500
                )
        else:
            not_found.append(ext)

    if not deleted:
        return PlainTextResponse(
            f"No files found for '{analysis_id_name}'", status_code=404
        )

    return JSONResponse(
        content={"deleted": deleted, "not_found": not_found},
        status_code=200,
    )
