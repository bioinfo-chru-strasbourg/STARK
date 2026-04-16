#!/usr/bin/env python

import json
import os
import subprocess
from typing import Optional, Union

from fastapi import APIRouter, Depends, HTTPException, Query, Request, status  # pyright: ignore[reportMissingImports]
from fastapi.responses import PlainTextResponse  # pyright: ignore[reportMissingImports]

from authentication import get_current_user_or_service
from config import docker_stark_api_log_folder, ts
from models import User
from queues import get_queue_env
from tasks import (
    _extract_analysis_name,
    queue_analysis,
    queue_command,
    queue_command_docker,
)

router = APIRouter()


@router.post("/analysis")
async def stark_launch(
    request: Request,
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
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


@router.post("/relaunch/{ts_id}")
async def relaunch_task(
    ts_id: str,
    queue: Optional[str] = Query(
        None,
        description="Queue name the task belongs to (default: first configured queue)",
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
