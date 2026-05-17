#!/usr/bin/env python

import json
import os
import subprocess
from typing import Optional, Union

from fastapi import APIRouter, Depends, HTTPException, Query, Request, status  # pyright: ignore[reportMissingImports]
from fastapi.responses import PlainTextResponse  # pyright: ignore[reportMissingImports]

import httpx  # pyright: ignore[reportMissingImports]

from authentication import get_current_user_or_service
from config import STARK_API_KEY, docker_stark_api_log_folder, ts
from models import User
from nodes import (
    compute_best_node,
    forward_request,
    get_local_metrics,
    get_nodes_metrics,
    resolve_self_url,
    load_nodes,
)
from queues import get_queue_env, load_queues
from modules import apply_module_defaults, resolve_module
from tasks import (
    _extract_analysis_name,
    queue_analysis,
    queue_command_docker,
    queue_command_docker_compose,
    queue_module_analysis,
)

router = APIRouter()


async def _run_locally(json_input: dict) -> str:
    """Dispatch json_input to the appropriate queue function and return the IDNAME."""

    # if "command" in json_input:
    #     return queue_command(json_input)
    # el
    if "command_docker" in json_input:
        return queue_command_docker(json_input)
    elif "command_docker_compose" in json_input:
        return queue_command_docker_compose(json_input)
    else:
        # Resolve module (raises ValueError for unknown names -> HTTP 400 in caller)
        module_cfg = resolve_module(json_input)
        # Apply module defaults for keys not already provided by the client
        apply_module_defaults(json_input, module_cfg)
        return queue_analysis(
            json_input,
            image=module_cfg["image"],
            module_docker_extra_params=module_cfg["docker_extra_params"],
            use_stark_container_mount=module_cfg["use_stark_container_mount"],
        )


async def _run_locally_module(json_input: dict) -> str:
    """Dispatch a module analysis locally (CLI command forwarded to container)."""
    module_cfg = resolve_module(json_input)
    apply_module_defaults(json_input, module_cfg)
    return queue_module_analysis(
        json_input,
        image=module_cfg["image"],
        module_docker_extra_params=module_cfg["docker_extra_params"],
        use_stark_container_mount=module_cfg["use_stark_container_mount"],
    )


@router.post("/analysis")
async def stark_launch(
    request: Request,
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """Launch a STARK analysis, command, or docker-command.

    Routing logic (transparent to the caller):
      1. If the request carries X-STARK-Forwarded (already routed once) ->
         run locally immediately to prevent loops.
      2. If no nodes are configured -> run locally.
      3. Collect metrics from all nodes + self asynchronously.
      4. Select the node with the most available capacity for the requested queue.
      5. If the best node is this node (or self-URL is unknown) -> run locally.
      6. Otherwise -> forward the request transparently to the best node.
    """
    if authorized != "service":
        if not isinstance(authorized, User) or "admin" not in authorized.groups:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin group required to launch analyses",
            )

    # --- Anti-loop guard: already forwarded once -> run locally immediately ---
    already_forwarded = request.headers.get("X-STARK-Forwarded", "0") == "1"

    try:
        json_input = await request.json()
    except Exception:
        json_input = {}

    # Queue
    queues = load_queues()
    if not queues:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="No queues configured",
        )
    default_queue = next(iter(queues))
    requested_queue = json_input.get("queue")
    if isinstance(requested_queue, str) and requested_queue.strip() != "":
        queue_name = requested_queue
    else:
        queue_name = default_queue

    json_input["queue"] = queue_name  # ensure queue is set for forwarded request

    # --- Routing (only on first hop when nodes are configured) ---
    if not already_forwarded and load_nodes():

        # Gather metrics: nodes + self
        nodes_metrics = await get_nodes_metrics()
        self_url = await resolve_self_url()
        if self_url:
            nodes_metrics[self_url] = get_local_metrics()
        target = compute_best_node(queue_name, nodes_metrics, json_input)

        if target:
            try:
                return await forward_request(target, request)
            except Exception:
                pass  # forward failed -> fall through to local execution

    # --- Run locally ---
    try:
        analysis_id_name = await _run_locally(json_input)
        return PlainTextResponse(content=analysis_id_name, status_code=200)
    except Exception as e:
        return PlainTextResponse(content=f"Launch failed: {e}", status_code=400)


@router.post("/analysis/module")
async def module_launch(
    request: Request,
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """Launch an analysis using a configured module.

    The ``command`` is always forwarded as direct CLI args to the container:
        docker run IMAGE <command>

    JSON input in the UI is a convenience — it is converted to a CLI string
    client-side before submission.

    Payload fields:
      - module        (str, optional): module name; defaults to the first configured module.
      - command       (str, required): CLI command string, e.g. ``--run=MY_RUN --sample_filter=S1,S2``.
      - analysis_name (str, optional): human-readable label for the analysis.
      - queue         (str, optional): queue name (defaults to first configured queue).
      - threads       (int, optional): number of CPU threads / task-spooler slots.
      - memory        (str, optional): Docker memory limit, e.g. ``4G``.
      - prioritize    (bool, optional): move to front of the queue.
    """
    if authorized != "service":
        if not isinstance(authorized, User) or "admin" not in authorized.groups:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin group required to launch analyses",
            )

    try:
        json_input = await request.json()
    except Exception:
        json_input = {}

    if not json_input.get("command"):
        return PlainTextResponse(
            content="'command' field is required", status_code=400
        )

    # Resolve module
    try:
        module_cfg = resolve_module(json_input)
    except ValueError as e:
        return PlainTextResponse(content=str(e), status_code=400)

    apply_module_defaults(json_input, module_cfg)

    # Queue
    queues = load_queues()
    if not queues:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="No queues configured",
        )
    default_queue = next(iter(queues))
    requested_queue = json_input.get("queue")
    if isinstance(requested_queue, str) and requested_queue.strip() != "":
        queue_name = requested_queue
    else:
        queue_name = default_queue
    json_input["queue"] = queue_name

    # --- Routing (same node logic as /analysis) ---
    already_forwarded = request.headers.get("X-STARK-Forwarded", "0") == "1"
    if not already_forwarded and load_nodes():
        nodes_metrics = await get_nodes_metrics()
        self_url = await resolve_self_url()
        if self_url:
            nodes_metrics[self_url] = get_local_metrics()
        target = compute_best_node(queue_name, nodes_metrics, json_input)
        if target:
            try:
                return await forward_request(target, request)
            except Exception:
                pass

    # --- Run locally ---
    try:
        analysis_id_name = queue_module_analysis(
            json_input,
            image=module_cfg["image"],
            module_docker_extra_params=module_cfg["docker_extra_params"],
            use_stark_container_mount=module_cfg["use_stark_container_mount"],
        )
        return PlainTextResponse(content=analysis_id_name, status_code=200)
    except Exception as e:
        return PlainTextResponse(content=f"Launch failed: {e}", status_code=400)


@router.post("/relaunch/{ts_id}")
async def relaunch_task(
    ts_id: str,
    queue: Optional[str] = Query(
        None,
        description="Queue name the task belongs to (default: first configured queue)",
    ),
    prioritize: Optional[bool] = Query(
        None,
        description="Whether to prioritize the task (default: False)",
    ),
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    if authorized != "service":
        if not isinstance(authorized, User) or "admin" not in authorized.groups:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin group required to relaunch analyses",
            )

    _ts_env = get_queue_env(queue)

    try:
        info_result = subprocess.run(
            f"{_ts_env} {ts} -i {ts_id}",
            shell=True,
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
    except subprocess.TimeoutExpired:
        raise HTTPException(
            status_code=504, detail="Failed to retrieve task info: command timed out"
        )
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to retrieve task info: internal error ({e})",
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

    if prioritize is not None:
        json_input["prioritize"] = prioritize

    try:
        relaunch_result = subprocess.run(
            f"{_ts_env} {ts} -r {ts_id}",
            shell=True,
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
        if relaunch_result.returncode != 0 or relaunch_result.stderr:
            return PlainTextResponse(
                content=f"Relaunch failed: {relaunch_result.stderr or 'unknown error'}",
                status_code=500,
            )
        # Verify the task was actually removed (ts -r can exit 0 silently on invalid IDs)
        verify_result = subprocess.run(
            f"{_ts_env} {ts} -i {ts_id}",
            shell=True,
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
        if verify_result.returncode == 0 and verify_result.stdout.strip():
            return PlainTextResponse(
                content=f"Relaunch failed: task {ts_id} was not removed from queue",
                status_code=500,
            )
    except subprocess.TimeoutExpired:
        return PlainTextResponse(
            content="Relaunch failed: command timed out",
            status_code=504,
        )
    except Exception as e:
        return PlainTextResponse(
            content=f"Relaunch failed: internal error ({e})",
            status_code=500,
        )

    # Queue
    queues = load_queues()
    if not queues:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="No queues configured",
        )
    default_queue = next(iter(queues))
    requested_queue = json_input.get("queue")
    if isinstance(requested_queue, str) and requested_queue.strip() != "":
        queue_name = requested_queue
    else:
        queue_name = default_queue

    json_input["queue"] = queue_name  # ensure queue is set for forwarded request

    # Write the updated JSON input back to the file for accurate forwarding if needed
    try:
        with open(json_file, "w") as f:
            json.dump(json_input, f, indent=2)
    except OSError as e:
        return PlainTextResponse(
            content=f"Relaunch failed: error writing JSON file ({e})",
            status_code=500,
        )

    # Route to the best node (same logic as a new analysis submission).
    if load_nodes():

        # Gather metrics: nodes + self
        nodes_metrics = await get_nodes_metrics()
        self_url = await resolve_self_url()
        if self_url:
            nodes_metrics[self_url] = get_local_metrics()
        target = compute_best_node(queue_name, nodes_metrics, json_input)

        if target:
            try:
                async with httpx.AsyncClient(timeout=15.0) as client:
                    resp = await client.post(
                        f"{target}/analysis",
                        json=json_input,
                        headers={"X-API-Key": STARK_API_KEY, "X-STARK-Forwarded": "1"},
                    )
                return PlainTextResponse(
                    content=resp.text, status_code=resp.status_code
                )
            except Exception:
                pass  # forward failed -> run locally

    # --- Run locally ---
    try:
        analysis_id_name = await _run_locally(json_input)
        return PlainTextResponse(content=analysis_id_name, status_code=200)
    except Exception as e:
        return PlainTextResponse(content=f"Relaunch failed: {e}", status_code=500)
