#!/usr/bin/env python

import json
import os
import random
import re
import string
import subprocess
from typing import Optional

from config import (
    docker_stark,
    docker_stark_api_log_folder,
    docker_stark_api_runs_folder,
    docker_stark_container_mount,
    ts,
)
from queues import _prepare_queue_for_submission
from security import (
    _sanitize_analysis_name,
    _validate_command,
    _validate_docker_extra_params,
    _validate_image,
)


def random_string_digits(string_length: int = 6) -> str:
    letters_and_digits = string.ascii_letters + string.digits
    return "".join(random.choice(letters_and_digits) for _ in range(string_length))


def _extract_analysis_name(ts_info_stdout: str) -> Optional[str]:
    """Extract analysisIDNAME from 'ts -i' output."""
    m = re.search(r"--name\s+(\S+)", ts_info_stdout)
    if m:
        return m.group(1)
    m = re.search(
        r"(STARK\.[A-Za-z0-9]+\.ID-[A-Za-z0-9]+-NAME-\S+?)\.output", ts_info_stdout
    )
    if m:
        return m.group(1)
    return None


def _resolve_task_slots(json_input: dict, max_slots: int) -> int:
    """Return the number of ts slots (-N) a task should consume."""
    raw = json_input.get("threads")
    if raw is None:
        return max_slots
    try:
        n = int(raw)
        if 0 <= n <= max_slots:
            return n
        return max_slots
    except (ValueError, TypeError):
        return max_slots


def _read_task_slots(command_str: str, max_slots: int) -> Optional[int]:
    """Resolve the effective slot count for a task from its stored JSON."""
    m = re.search(r">\s*(\S+)\.output\s+2>&1", command_str)
    if not m:
        return None
    json_path = m.group(1) + ".json"
    try:
        with open(json_path) as f:
            data = json.load(f)
        return _resolve_task_slots(data, max_slots)
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return None


def _build_analysis_idname(json_input: dict) -> tuple:
    """Build analysisIDNAME and analysesRUNNAME from json_input."""
    analyses_id = random_string_digits(12)
    run_id = "UNKNOWN"
    run_md5 = random_string_digits(41)

    if "run" in json_input:
        name = json_input["run"].split(":")[0]
        run_id = os.path.basename(name)

        run_folder = ""
        if os.path.isdir(name):
            run_folder = name
        elif os.path.isdir(os.path.join(docker_stark_api_runs_folder, name)):
            run_folder = os.path.join(docker_stark_api_runs_folder, name)

        if run_folder:
            my_cmd = f"find {run_folder} -maxdepth 1 -type f -print0 | xargs -0 sha1sum | cut -b-40 | sha1sum | awk '{{print $1}}'"
        else:
            my_cmd = f"echo {name} | sha1sum | awk '{{print $1}}'"

        run_md5 = (
            subprocess.run(my_cmd, shell=True, stdout=subprocess.PIPE)
            .stdout.decode("utf-8")
            .strip()
        )

        if "analysis_name" in json_input:
            run_id = _sanitize_analysis_name(str(json_input["analysis_name"]))

    elif "analysis_name" in json_input:
        run_id = _sanitize_analysis_name(str(json_input["analysis_name"]))
        run_md5 = random_string_digits(41)

    analyses_run_name = run_id
    analyses_name = f"ID-{run_md5}-NAME-{run_id}"
    analysis_id_name = f"STARK.{analyses_id}.{analyses_name}"
    return analysis_id_name, analyses_run_name


def _queue_task(my_cmd: str, prioritize: bool = False, ts_cmd: str = "") -> str:
    task_id = (
        subprocess.run(my_cmd, shell=True, stdout=subprocess.PIPE)
        .stdout.decode("utf-8")
        .strip()
    )
    if not task_id:
        raise RuntimeError("task-spooler returned no task ID")

    if prioritize:
        # Move the task to the front of the queue
        _ = (
            subprocess.run(f"{ts_cmd} -u {task_id}", shell=True, stdout=subprocess.PIPE)
            .stdout.decode("utf-8")
            .strip()
        )
    return task_id


def queue_analysis(json_input: dict) -> str:
    """Build and queue a STARK Docker analysis. Returns analysisIDNAME or raises RuntimeError."""
    analysis_id_name, analyses_run_name = _build_analysis_idname(json_input)

    # Docker parameters

    # Use analysis_id_name as Docker container name to be able to identify the container corresponding to a task-spooler task, and stop it if needed when the task is killed.
    docker_name = f" --name {analysis_id_name} "

    # Extra docker parameters
    docker_extra_params = json_input.get("docker_extra_params", "")
    if docker_extra_params:
        _validate_docker_extra_params(docker_extra_params)

    # Combine Docker parameters, ensuring --name is included and --rm is set for cleanup.
    docker_parameters = (
        f"--rm {docker_name} {docker_extra_params} {docker_stark_container_mount}"
    )
    # if json_input.get("use_stark_container_mount", True):
    #     docker_parameters += f" {docker_stark_container_mount}"

    analysis_folder = docker_stark_api_log_folder
    analysis_file = os.path.join(analysis_folder, f"{analysis_id_name}.json")
    analysis_file_stark = os.path.join(
        analysis_folder, f"{analysis_id_name}.json.stark_analysis"
    )
    analysis_info_file = os.path.join(analysis_folder, f"{analysis_id_name}.info")
    analysis_output_file = os.path.join(analysis_folder, f"{analysis_id_name}.output")

    _ts_env, _max_slots = _prepare_queue_for_submission(json_input.get("queue"))
    _task_slots = _resolve_task_slots(json_input, _max_slots)

    threads = (_task_slots == 0) and _max_slots or _task_slots
    json_input["threads"] = threads

    # CPU/Threads
    if "--cpus" not in docker_parameters:
        docker_parameters += f" --cpus={threads} "

    # Memory
    if (
        "--memory" not in docker_parameters
        and "memory" in json_input
        and json_input.get("memory", None)
    ):
        docker_parameters += f" --memory={json_input.get('memory', '')} "

    # Prioritize
    prioritize = json_input.get("prioritize", False)

    # Write the JSON
    with open(analysis_file, "w") as f:
        f.write(json.dumps(json_input))

    # Remove host-side task scheduling parameters from the JSON passed to the container.
    json_input_for_container = json_input.copy()
    forbidden_params = [
        "queue",
        "docker_extra_params",
        "use_stark_container_mount",
        "memory",
        "prioritize",
    ]
    for param in forbidden_params:
        json_input_for_container.pop(param, None)

    # Write the final JSON input for the container, which may be used for metrics and debugging.
    with open(analysis_file_stark, "w") as f:
        f.write(json.dumps(json_input_for_container))

    ts_cmd = f"{_ts_env}{ts} -N {_task_slots} -L {analysis_id_name}" if ts else ""
    my_cmd = (
        f'{ts_cmd} bash -c "'
        f"trap 'docker stop {analysis_id_name} 2>/dev/null; echo failed > {analysis_info_file}' TERM INT; "
        f"docker run {docker_parameters} {docker_stark} "
        f"--analysis_name={analyses_run_name} --analysis={analysis_file_stark} "
        f"> {analysis_output_file} 2>&1 "
        f"&& (echo 'finished' > {analysis_info_file} && exit 0) "
        f"|| (echo 'failed' > {analysis_info_file} && exit 1)\""
    )

    # Queue the task and get the task ID from task-spooler
    _ = _queue_task(my_cmd, prioritize=prioritize, ts_cmd=f"{_ts_env}{ts}")

    return analysis_id_name


def queue_command(json_input: dict) -> str:
    """Queue a raw shell command (json key 'command'). Returns analysisIDNAME or raises RuntimeError."""
    command = json_input["command"]
    _validate_command(command)
    analysis_id_name, _ = _build_analysis_idname(json_input)

    analysis_folder = docker_stark_api_log_folder
    analysis_file = os.path.join(analysis_folder, f"{analysis_id_name}.json")
    analysis_info_file = os.path.join(analysis_folder, f"{analysis_id_name}.info")
    analysis_output_file = os.path.join(analysis_folder, f"{analysis_id_name}.output")

    _ts_env, _max_slots = _prepare_queue_for_submission(json_input.get("queue"))
    _task_slots = _resolve_task_slots(json_input, _max_slots)

    with open(analysis_file, "w") as f:
        f.write(json.dumps(json_input))

    ts_cmd = f"{_ts_env}{ts} -N {_task_slots} -L {analysis_id_name}" if ts else ""
    my_cmd = (
        f'{ts_cmd} bash -c "'
        f"trap 'echo failed > {analysis_info_file}' TERM INT; "
        f"{command} "
        f"> {analysis_output_file} 2>&1 "
        f"&& (echo 'finished' > {analysis_info_file} && exit 0) "
        f"|| (echo 'failed' > {analysis_info_file} && exit 1)\""
    )

    # Queue the task and get the task ID from task-spooler
    _ = _queue_task(
        my_cmd, prioritize=json_input.get("prioritize", False), ts_cmd=f"{_ts_env}{ts}"
    )

    return analysis_id_name


def queue_command_docker(json_input: dict) -> str:
    """Queue a docker command with predefined mount parameters (json key 'command_docker').
    Returns analysis_id_name or raises RuntimeError."""
    command_docker = json_input["command_docker"]
    image = json_input.get("image")
    if not image:
        raise ValueError("'image' is required when using 'command_docker'")
    _validate_image(image)
    docker_extra_params = json_input.get("docker_extra_params", "")
    if docker_extra_params:
        _validate_docker_extra_params(docker_extra_params)

    analysis_id_name, _ = _build_analysis_idname(json_input)

    docker_name = f"--name {analysis_id_name}"
    docker_parameters = f"--rm {docker_name} {docker_extra_params}"
    if json_input.get("use_stark_container_mount", True):
        docker_parameters += f" {docker_stark_container_mount}"

    analysis_folder = docker_stark_api_log_folder
    analysis_file = os.path.join(analysis_folder, f"{analysis_id_name}.json")
    analysis_info_file = os.path.join(analysis_folder, f"{analysis_id_name}.info")
    analysis_output_file = os.path.join(analysis_folder, f"{analysis_id_name}.output")

    _ts_env, _max_slots = _prepare_queue_for_submission(json_input.get("queue"))
    _task_slots = _resolve_task_slots(json_input, _max_slots)

    threads = (_task_slots == 0) and _max_slots or _task_slots
    json_input["threads"] = threads

    with open(analysis_file, "w") as f:
        f.write(json.dumps(json_input))

    # CPU/Threads
    if "--cpus" not in docker_parameters:
        docker_parameters += f" --cpus={threads} "

    # Memory
    if (
        "--memory" not in docker_parameters
        and "memory" in json_input
        and json_input.get("memory", None)
    ):
        docker_parameters += f" --memory={json_input.get('memory', '')} "

    ts_cmd = f"{_ts_env}{ts} -N {_task_slots} -L {analysis_id_name}" if ts else ""
    my_cmd = (
        f'{ts_cmd} bash -c "'
        f"trap 'docker stop {analysis_id_name} 2>/dev/null; echo failed > {analysis_info_file}' TERM INT; "
        f"docker run {docker_parameters} {image} {command_docker} "
        f"> {analysis_output_file} 2>&1 "
        f"&& (echo 'finished' > {analysis_info_file} && exit 0) "
        f"|| (echo 'failed' > {analysis_info_file} && exit 1)\""
    )

    # Queue the task and get the task ID from task-spooler
    _ = _queue_task(
        my_cmd, prioritize=json_input.get("prioritize", False), ts_cmd=f"{_ts_env}{ts}"
    )

    return analysis_id_name


def queue_command_docker_compose(json_input: dict) -> str:
    """Queue a docker-compose command with predefined parameters (json key 'command_docker_compose').
    Returns analysis_id_name or raises RuntimeError."""

    # Command
    command_docker_compose = json_input["command_docker_compose"]

    # Service
    service = json_input.get("service", json_input.get("image"))
    if not service:
        raise ValueError("'service' is required when using 'command_docker_compose'")
    # _validate_image(service)

    # docker configuration file
    docker_compose_file = json_input.get("docker_compose_file")
    if not docker_compose_file:
        raise ValueError(
            "'docker_compose_file' is required when using 'command_docker_compose'"
        )
    if not os.path.isfile(docker_compose_file):
        raise ValueError(f"'docker_compose_file' does not exist: {docker_compose_file}")

    # Extra docker parameters
    docker_extra_params = json_input.get("docker_extra_params", "")
    if docker_extra_params:
        _validate_docker_extra_params(docker_extra_params)

    analysis_id_name, _ = _build_analysis_idname(json_input)

    docker_name = f"--name {analysis_id_name}"
    docker_parameters = f"--rm {docker_name} {docker_extra_params}"
    if json_input.get("use_stark_container_mount", True):
        docker_parameters += f" {docker_stark_container_mount}"

    analysis_folder = docker_stark_api_log_folder
    analysis_file = os.path.join(analysis_folder, f"{analysis_id_name}.json")
    analysis_info_file = os.path.join(analysis_folder, f"{analysis_id_name}.info")
    analysis_output_file = os.path.join(analysis_folder, f"{analysis_id_name}.output")

    # Touch output
    with open(analysis_output_file, "w") as f:
        f.write("")

    _ts_env, _max_slots = _prepare_queue_for_submission(json_input.get("queue"))
    _task_slots = _resolve_task_slots(json_input, _max_slots)

    threads = (_task_slots == 0) and _max_slots or _task_slots
    json_input["threads"] = threads

    with open(analysis_file, "w") as f:
        f.write(json.dumps(json_input))

    ts_cmd = f"{_ts_env}{ts} -N {_task_slots} -L {analysis_id_name}" if ts else ""

    my_cmd = (
        f'{ts_cmd} bash -c "'
        f"trap 'docker stop {analysis_id_name} 2>/dev/null; echo failed > {analysis_info_file}' TERM INT; "
        f"docker-compose -f {docker_compose_file} run {docker_parameters} {service} {command_docker_compose} "
        f"> {analysis_output_file} 2>&1 "
        f"&& (echo 'finished' > {analysis_info_file} && exit 0) "
        f"|| (echo 'failed' > {analysis_info_file} && exit 1)\""
    )

    # Queue the task and get the task ID from task-spooler
    _ = _queue_task(
        my_cmd, prioritize=json_input.get("prioritize", False), ts_cmd=f"{_ts_env}{ts}"
    )

    return analysis_id_name
