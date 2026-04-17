#!/usr/bin/env python

import concurrent.futures
import datetime
import json
import os
import re
import subprocess
import time
from typing import Optional, Union

from fastapi import APIRouter, Depends, HTTPException, Query, status  # pyright: ignore[reportMissingImports]
from fastapi.responses import JSONResponse, PlainTextResponse  # pyright: ignore[reportMissingImports]

from authentication import get_current_user_or_service
from config import docker_stark_api_log_folder, shell, ts
from models import User
from queues import (
    _queue_daemon_is_active,
    _resolve_queue_socket,
    get_queue_env,
    load_queues,
)
from tasks import _extract_analysis_name, _read_task_slots

router = APIRouter()

_LINE_REGEX = re.compile(
    r"^(?P<id>\d+)\s+"
    r"(?P<state>\w+)\s+"
    r"(?P<output>\S+)\s+"
    r"(?P<elevel>\S*)\s*"
    r"(?P<times>[\d./\s-]*)\s+"
    r"(?P<command>.*)$"
)


def _iter_queue_tasks(queue_name: str, cfg: dict, is_default: bool) -> list:
    """Run 'ts -l' on a single queue and return parsed task dicts."""
    socket = _resolve_queue_socket(queue_name, cfg, is_default)
    socket_part = f"TS_SOCKET={socket} " if socket else ""
    _q_env = f"{socket_part}TS_SAVELIST={cfg['savelist']} TS_SLOTS={cfg['slots']} "

    try:
        result = subprocess.run(
            f"{_q_env} {ts} -l",
            shell=True,
            executable=shell,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        raise HTTPException(status_code=500, detail=f"Command not found: {ts}")

    if result.returncode != 0 and (
        "tasks" in result.stderr.lower() or not result.stdout
    ):
        return []

    tasks = []
    lines = result.stdout.strip().split("\n")
    for line in lines[1:]:
        match = _LINE_REGEX.match(line)
        if not match:
            continue
        data = match.groupdict()
        run_name_match = re.search(r"-NAME-([^\s\]]+)", data["command"])
        run_name = re.sub(
            r"\.(output|json|info)$",
            "",
            run_name_match.group(1) if run_name_match else "N/A",
        )
        task_slots = _read_task_slots(data["command"], int(cfg.get("slots", 1)))
        tasks.append(
            {
                "id": data["id"].strip(),
                "state": data["state"].strip(),
                "output": data["output"].strip(),
                "elevel": (
                    "SUCCESS"
                    if data["elevel"].strip() == "0"
                    else (
                        "FAILED: " + data["elevel"].strip()
                        if data[
                            "elevel"
                        ].strip()  # any non-zero, non-empty → failed (handles "signal:15" etc.)
                        else ""
                    )
                ),
                "times": "",
                "run_name": run_name,
                "queue": queue_name,
                "task_slots": task_slots,
                "queue_slots": int(cfg.get("slots", 1)),
                "_ts_env": _q_env,
            }
        )
    return tasks


def _enrich_task_times(task: dict) -> None:
    """Fetch elapsed/finished time for a running or finished task (in-place)."""
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


@router.get("/list")
async def list_task():
    all_tasks = []
    queues = load_queues()
    default_name = next(iter(queues))
    for queue_name, cfg in queues.items():
        is_default = queue_name == default_name
        if not _queue_daemon_is_active(queue_name, cfg, is_default):
            continue
        tasks = _iter_queue_tasks(queue_name, cfg, is_default)
        all_tasks.extend(tasks)

    # Enrich times for running/finished tasks (same as /queue?action=list)
    tasks_to_enrich = [
        t for t in all_tasks if t["state"].lower() in ("running", "finished")
    ]
    if tasks_to_enrich:
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
            executor.map(_enrich_task_times, tasks_to_enrich)

    # return (elevel && elevel.startsWith('FAILED')) ? 'state-finished-failed' : 'state-finished-success';

    result = []
    for t in all_tasks:
        result.append(
            {
                "id": int(t["id"]),
                # "state": t["state"],
                "state": (
                    "failed"
                    if t["state"].lower() == "finished"
                    and t["elevel"]
                    and t["elevel"].startswith("FAILED")
                    else t["state"].lower()
                ),
                # elevel is already "SUCCESS" / "FAILED: X" / "" from _iter_queue_tasks
                "elevel": t["elevel"] or None,
                "times": (float(t["times"].split("/")[0]) if t["times"] else None),
                "run_name": t["run_name"],
                "queue": t["queue"],
                "task_slots": t["task_slots"],
                "queue_slots": t["queue_slots"],
            }
        )
    return JSONResponse(content=result)


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

    max_entries = 1000

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

        result.append(
            {
                "analysis_id_name": base,
                "run_name": run_name,
                "queue": queue,
                "threads": threads,
                "status": status,
                "mtime": os.path.getmtime(json_path),
            }
        )

    return JSONResponse(content=result)


@router.get("/queue")
async def queue_action(
    action: str = Query(
        "list",
        enum=["list", "info", "log", "analysis", "kill", "prioritize", "remove", "swap"],
    ),
    id: Optional[str] = Query(None),
    queue: Optional[str] = Query(
        None, description="Queue name (default: first configured queue)"
    ),
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    _ts_env = get_queue_env(queue)

    if action == "list":
        queues_all = load_queues()
        default_name_q = next(iter(queues_all))
        queues_to_scan = (
            {queue: queues_all.get(queue, {"savelist": f"/ts-tmp-{queue}", "slots": 1})}
            if queue
            else queues_all
        )
        all_tasks = []
        for queue_name, cfg in queues_to_scan.items():
            is_default_q = queue_name == default_name_q
            if not _queue_daemon_is_active(queue_name, cfg, is_default_q):
                continue
            all_tasks.extend(_iter_queue_tasks(queue_name, cfg, is_default_q))

        tasks_to_enrich = [
            t for t in all_tasks if t["state"].lower() in ("running", "finished")
        ]
        if tasks_to_enrich:
            with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
                executor.map(_enrich_task_times, tasks_to_enrich)

        for t in all_tasks:
            t.pop("_ts_env", None)

        return JSONResponse(content=all_tasks)

    # --- Destructive / single-task actions ---
    restricted_actions = {"kill", "prioritize", "remove", "swap"}
    if action in restricted_actions and authorized != "service":
        if not isinstance(authorized, User) or "admin" not in authorized.groups:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin group required for this action",
            )

    if not id:
        raise HTTPException(
            status_code=400, detail="Task ID is required for this action"
        )

    if action == "kill":
        info_result = subprocess.run(
            f"{_ts_env} {ts} -i {id}",
            shell=True,
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
        container_match = re.search(r"--name\s+(\S+)", info_result.stdout)
        if container_match:
            subprocess.run(
                f"docker stop {container_match.group(1)}",
                shell=True,
                capture_output=True,
                text=True,
                check=False,
                timeout=30,
            )

    elif action == "analysis":
        info_result = subprocess.run(
            f"{_ts_env} {ts} -i {id}",
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
        json_file = os.path.join(docker_stark_api_log_folder, f"{analysis_name}.json")
        if not os.path.isfile(json_file):
            return PlainTextResponse(
                f"JSON file not found: {json_file}", status_code=404
            )
        try:
            with open(json_file, "r", errors="replace") as f:
                content = f.read()
            try:
                content = json.dumps(json.loads(content), indent=2, ensure_ascii=False)
            except json.JSONDecodeError:
                pass
            return PlainTextResponse(content)
        except OSError as e:
            return PlainTextResponse(f"Error reading JSON file: {e}", status_code=500)

    elif action == "log":
        info_result = subprocess.run(
            f"{_ts_env} {ts} -i {id}",
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
            return PlainTextResponse(f"Error reading output file: {e}", status_code=500)

    action_map = {"info": "-i", "kill": "-k", "prioritize": "-u", "remove": "-r", "swap": "-U"}
    ts_action = action_map.get(action)
    if not ts_action:
        raise HTTPException(status_code=400, detail=f"Invalid action: {action}")

    try:
        result = subprocess.run(
            f"{_ts_env} {ts} {ts_action} {id}",
            shell=True,
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
        if result.returncode != 0:
            return PlainTextResponse(
                result.stderr or f"Command failed (exit {result.returncode})",
                status_code=500,
            )
        return PlainTextResponse(
            result.stdout or f"OK ({action}#{id})", status_code=200
        )
    except subprocess.TimeoutExpired:
        return PlainTextResponse("Command timed out.", status_code=504)
    except FileNotFoundError:
        return PlainTextResponse(f"Command not found: {ts}", status_code=500)
