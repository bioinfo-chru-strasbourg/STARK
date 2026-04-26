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
from peers import (
    compute_best_peer,
    forward_request,
    get_local_metrics,
    get_peers_metrics,
    get_self_url,
    load_peers,
)
from queues import get_queue_env, load_queues
from tasks import (
    _extract_analysis_name,
    queue_analysis,
    queue_command,
    queue_command_docker,
    queue_command_docker_compose,
)

router = APIRouter()


async def _run_locally(json_input: dict) -> str:
    """Dispatch json_input to the appropriate queue function and return the IDNAME."""

    # Retrieve specific parameters for resources, and create a separate json/dict
    # to pass to the queue functions (to avoid passing irrelevant parameters).

    # specific_params_list = ["threads", "memory", "prioritize"]

    # specific_params = {
    #     param: json_input[param]
    #     for param in specific_params_list
    #     if param in json_input
    # }

    if "command" in json_input:
        return queue_command(json_input)
    elif "command_docker" in json_input:
        return queue_command_docker(json_input)
    elif "command_docker_compose" in json_input:
        return queue_command_docker_compose(json_input)
    else:
        return queue_analysis(json_input)


@router.post("/analysis")
async def stark_launch(
    request: Request,
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """Launch a STARK analysis, command, or docker-command.

    Routing logic (transparent to the caller):
      1. If the request carries X-STARK-Forwarded (already routed once) ->
         run locally immediately to prevent loops.
      2. If no peers are configured -> run locally.
      3. Collect metrics from all peers + self asynchronously.
      4. Select the peer with the most available capacity for the requested queue.
      5. If the best peer is this node (or self-URL is unknown) -> run locally.
      6. Otherwise -> forward the request transparently to the best peer.
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

    # --- Routing (only on first hop when peers are configured) ---
    if not already_forwarded and load_peers():
        queues = load_queues()
        default_queue = next(iter(queues))
        queue_name = json_input.get("queue") or default_queue

        # Gather metrics: peers + self
        peers_metrics = await get_peers_metrics()
        self_url = get_self_url()
        if self_url:
            peers_metrics[self_url] = get_local_metrics()
        target = compute_best_peer(queue_name, peers_metrics, json_input)

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

    if prioritize is not None:
        json_input["prioritize"] = prioritize
        print(
            f"Set prioritize={prioritize} for task {ts_id} based on query parameter {_ts_env}."
        )
        info_result = subprocess.run(
            f"{_ts_env} {ts} -r {ts_id}",
            shell=True,
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )

    # Route to the best peer (same logic as a new analysis submission).
    if load_peers():
        queues = load_queues()
        default_queue = next(iter(queues))
        queue_name = json_input.get("queue") or default_queue

        # Gather metrics: peers + self
        peers_metrics = await get_peers_metrics()
        self_url = get_self_url()
        if self_url:
            peers_metrics[self_url] = get_local_metrics()
        target = compute_best_peer(queue_name, peers_metrics, json_input)

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
