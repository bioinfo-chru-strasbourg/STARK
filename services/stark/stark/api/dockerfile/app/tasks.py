#!/usr/bin/env python

import json
import os
import random
import re
import shlex
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
    _sanitize_docker_extra_params,
    _validate_image,
    _safe_split,
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
        if 1 <= n <= max_slots:
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


# def _build_analysis_idname(json_input: dict) -> tuple:
#     """Build analysisIDNAME and analysesRUNNAME from json_input."""
#     analyses_id = random_string_digits(12)
#     run_id = "UNKNOWN"
#     run_md5 = random_string_digits(41)

#     if "run" in json_input:
#         name = json_input["run"].split(":")[0]
#         run_id = os.path.basename(name)

#         run_folder = ""
#         if os.path.isdir(name):
#             run_folder = name
#         elif os.path.isdir(os.path.join(docker_stark_api_runs_folder, name)):
#             run_folder = os.path.join(docker_stark_api_runs_folder, name)

#         if run_folder:
#             my_cmd = f"find {run_folder} -maxdepth 1 -type f -print0 | xargs -0 sha1sum | cut -b-40 | sha1sum | awk '{{print $1}}'"
#         else:
#             my_cmd = f"echo {name} | sha1sum | awk '{{print $1}}'"

#         run_md5 = (
#             subprocess.run(my_cmd, shell=True, stdout=subprocess.PIPE)
#             .stdout.decode("utf-8")
#             .strip()
#         )

#         if "analysis_name" in json_input:
#             run_id = _sanitize_analysis_name(str(json_input["analysis_name"]))

#     elif "analysis_name" in json_input:
#         run_id = _sanitize_analysis_name(str(json_input["analysis_name"]))
#         run_md5 = random_string_digits(41)

#     analyses_run_name = run_id
#     analyses_name = f"ID-{run_md5}-NAME-{run_id}"
#     analysis_id_name = f"STARK.{analyses_id}.{analyses_name}"
#     return analysis_id_name, analyses_run_name


def _build_analysis_idname(json_input: dict) -> tuple:
    """Build analysisIDNAME and analysesRUNNAME from json_input."""

    # Analysis ID
    analyses_id = random_string_digits(12)
    # run_id = "UNKNOWN_" + str(random_string_digits(10))
    # run_md5 = random_string_digits(41)
    import datetime

    # RUN ID
    if "analysis_name" in json_input:
        run_id = _sanitize_analysis_name(str(json_input["analysis_name"]))
    elif "run" in json_input:
        run_id = _sanitize_analysis_name(str(json_input["run"].split(":")[0]))
    elif "command" in json_input:
        if isinstance(json_input["command"], str):
            # Extract run name from the command string if possible, using the parameter "--run=<run_id> " (e.g. "--run=MY_RUN" or "--run MY_RUN")
            analysis_name_match = re.search(r"--analysis_name[=\s]+(\S+)", json_input["command"])
            run_id_match = re.search(r"--run[=\s]+(\S+)", json_input["command"])
            if analysis_name_match:
                run_id = _sanitize_analysis_name(analysis_name_match.group(1))
            elif run_id_match:
                run_id = _sanitize_analysis_name(run_id_match.group(1))
            else:
                run_id = _sanitize_analysis_name(str(json_input["command"].split(":")[0]))
        elif isinstance(json_input["command"], dict) and len(json_input["command"]) > 0:
            # Extract run name from the command dict if possible, using the key "run" (e.g. {"run": "MY_RUN"})
            analysis_name = json_input["command"].get("analysis_name")
            run = json_input["command"].get("run")
            if analysis_name:
                run_id = _sanitize_analysis_name(str(analysis_name))
            elif run:
                run_id = _sanitize_analysis_name(str(run))
            else:
                run_id = _sanitize_analysis_name(str(json_input["command"].get("run", "UNKNOWN")))
        else:
            run_id = (
                datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                + "_"
                + "UNKNOWN_"
                + str(random_string_digits(10))
            )
    else:
        # run_id = "UNKNOWN_" + str(random_string_digits(10))
        run_id = (
            datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            + "_"
            + "UNKNOWN_"
            + str(random_string_digits(10))
        )

    # MD5
    run_md5 = random_string_digits(10)

    # if "run" in json_input:
    #     name = json_input["run"].split(":")[0]
    #     run_id = os.path.basename(name)

    #     run_folder = ""
    #     if os.path.isdir(name):
    #         run_folder = name
    #     elif os.path.isdir(os.path.join(docker_stark_api_runs_folder, name)):
    #         run_folder = os.path.join(docker_stark_api_runs_folder, name)

    #     if run_folder:
    #         my_cmd = f"find {run_folder} -maxdepth 1 -type f -print0 | xargs -0 sha1sum | cut -b-40 | sha1sum | awk '{{print $1}}'"
    #     else:
    #         my_cmd = f"echo {name} | sha1sum | awk '{{print $1}}'"

    #     run_md5 = (
    #         subprocess.run(my_cmd, shell=True, stdout=subprocess.PIPE)
    #         .stdout.decode("utf-8")
    #         .strip()
    #     )

    #     if "analysis_name" in json_input:
    #         run_id = _sanitize_analysis_name(str(json_input["analysis_name"]))

    # elif "analysis_name" in json_input:
    #     run_id = _sanitize_analysis_name(str(json_input["analysis_name"]))
    #     run_md5 = random_string_digits(41)

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
        try:
            reprioritize_result = subprocess.run(
                f"{ts_cmd} -u {task_id}",
                shell=True,
                capture_output=True,
                text=True,
                timeout=10,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                f"timed out after 10 seconds while reprioritizing task {task_id} "
                f"with '{ts_cmd} -u {task_id}'"
            ) from exc
        if reprioritize_result.returncode != 0:
            error_output = (
                reprioritize_result.stderr or reprioritize_result.stdout or ""
            ).strip()
            if error_output:
                # ts returns "cannot be urged" when the task is already running or finished:
                # this is expected and not an error.
                if "cannot be urged" in error_output:
                    pass
                else:
                    raise RuntimeError(
                        f"failed to reprioritize task {task_id} with '{ts_cmd} -u {task_id}': {error_output} | stdout: {reprioritize_result.stdout.strip()} | stderr: {reprioritize_result.stderr.strip()}"
                    )
            # No output and non-zero exit: task is likely already running or finished,
            # reprioritization is not applicable but not an error.
    return task_id


def queue_analysis(
    json_input: dict,
    image: str = "",
    module_docker_extra_params: str = "",
    use_stark_container_mount: bool = True,
) -> str:
    """Build and queue a STARK Docker analysis. Returns analysisIDNAME or raises RuntimeError.

    Args:
        json_input: The analysis parameters from the client.
        image: Docker image to use. Resolved from modules.json by the caller;
               falls back to the ``docker_stark`` config variable when empty.
        module_docker_extra_params: Extra Docker flags defined by the module
               (e.g. ``--shm-size=4g``). Sanitized before use.
        use_stark_container_mount: Whether to append the STARK volume mount
               string to the Docker parameters.
    """
    if not image:
        image = docker_stark

    analysis_id_name, analyses_run_name = _build_analysis_idname(json_input)

    # Docker parameters

    # Use analysis_id_name as Docker container name to be able to identify the container corresponding to a task-spooler task, and stop it if needed when the task is killed.
    docker_name = f" --name {analysis_id_name} "

    # Extra docker parameters from the module definition (server-side)
    if module_docker_extra_params:
        _validate_docker_extra_params(module_docker_extra_params)
        module_docker_extra_params = _sanitize_docker_extra_params(module_docker_extra_params)

    # Extra docker parameters from the client request (ignored for module analyses — kept for
    # backward-compat with direct callers but stripped before reaching the container)
    docker_extra_params = ""

    # Combine Docker parameters, ensuring --name is included and --rm is set for cleanup.
    mount = docker_stark_container_mount if use_stark_container_mount else ""
    docker_parameters = (
        f"--rm {docker_name} {module_docker_extra_params} {mount}"
    )

    # Analysis file paths
    analysis_folder = docker_stark_api_log_folder
    analysis_file = os.path.join(analysis_folder, f"{analysis_id_name}.json")
    analysis_file_stark = os.path.join(
        analysis_folder, f"{analysis_id_name}.json.stark_analysis"
    )
    analysis_info_file = os.path.join(analysis_folder, f"{analysis_id_name}.info")
    analysis_output_file = os.path.join(analysis_folder, f"{analysis_id_name}.output")

    # Task spooler configuration
    _ts_env, _max_slots = _prepare_queue_for_submission(json_input.get("queue"))
    _task_slots = _resolve_task_slots(json_input, _max_slots)

    # Threads
    threads = (_task_slots == 0) and _max_slots or _task_slots
    json_input["threads"] = threads

    # CPU/Threads
    docker_extra_params = _sanitize_docker_extra_params(
        docker_extra_params, extra_forbidden={"--cpus": True}
    )
    docker_parameters += f" --cpus={threads} "

    # Memory
    if "memory" in json_input and json_input.get("memory", None):
        docker_extra_params = _sanitize_docker_extra_params(
            docker_extra_params, extra_forbidden={"--memory": True, "-m": True}
        )
        memory = str(json_input.get("memory", "")).strip()
        if not re.fullmatch(r"\d+(?:[bBkKmMgG])?", memory):
            raise ValueError(f"Invalid memory value: {memory}")
        docker_parameters += f" --memory={memory} "

    # Prioritize
    prioritize = json_input.get("prioritize", False)

    # Write the JSON
    with open(analysis_file, "w") as f:
        f.write(json.dumps(json_input))

    # Remove host-side task scheduling parameters from the JSON passed to the container.
    json_input_for_container = json_input.copy()
    forbidden_params = [
        "module",
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

    # Build the docker run command with the module image.
    docker_cmd = ["docker", "run"]
    docker_cmd += _safe_split(docker_parameters)
    docker_cmd.append(image)
    docker_cmd_str = " ".join(shlex.quote(tok) for tok in docker_cmd)

    # Task spooler command
    ts_cmd = f"{_ts_env}{ts} -N {_task_slots} -L {analysis_id_name}" if ts else ""

    # Task spooler command: we wrap the docker run in a bash command that traps termination signals to stop the container and mark the task as failed, and writes "finished" or "failed" to an info file based on the exit status of the docker command. This allows us to track the status of the task and ensure cleanup if it's killed.
    my_cmd = (
        f'{ts_cmd} bash -c "'
        f"trap 'docker stop {analysis_id_name} 2>/dev/null; echo failed > {analysis_info_file}' TERM INT; "
        f"{docker_cmd_str} "
        f"--analysis_name={analyses_run_name} --analysis={analysis_file_stark} "
        f"> {analysis_output_file} 2>&1 "
        f"&& (echo 'finished' > {analysis_info_file} && exit 0) "
        f"|| (echo 'failed' > {analysis_info_file} && exit 1)\""
    )

    # Queue the task and get the task ID from task-spooler
    _ = _queue_task(my_cmd, prioritize=prioritize, ts_cmd=f"{_ts_env}{ts}")

    return analysis_id_name


def _parse_command_to_dict(command: str) -> dict:
    """Parse a CLI command string into a dict of key/value pairs.

    Examples:
        "--run=MY_RUN --sample_filter=S1,S2" -> {"run": "MY_RUN", "sample_filter": "S1,S2"}
        "--flag"                              -> {"flag": True}
        "--key value"                         -> {"key": "value"}
    """
    result = {}
    try:
        tokens = shlex.split(command)
    except ValueError:
        tokens = command.split()
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok.startswith("--"):
            if "=" in tok:
                k, v = tok[2:].split("=", 1)
                result[k] = v
            elif i + 1 < len(tokens) and not tokens[i + 1].startswith("-"):
                result[tok[2:]] = tokens[i + 1]
                i += 1
            else:
                result[tok[2:]] = True
        i += 1
    return result


def queue_module_analysis(
    json_input: dict,
    image: str,
    module_docker_extra_params: str = "",
    use_stark_container_mount: bool = True,
) -> str:
    """Queue a module analysis by passing the ``command`` field as direct CLI args.

    The ``command`` value is always forwarded as-is to the container:
        docker run IMAGE <command>

    The caller is responsible for converting any JSON-based input to a CLI
    string before calling this function (see ``jsonCommandToCli`` on the
    client side).

    Args:
        json_input: Must contain a ``command`` key with the CLI string.
        image: Docker image to use (resolved from modules.json by the caller).
        module_docker_extra_params: Extra Docker flags defined by the module.
        use_stark_container_mount: Whether to append the STARK volume mount.
    """
    if not image:
        image = docker_stark

    analysis_id_name, _ = _build_analysis_idname(json_input)

    docker_name = f" --name {analysis_id_name} "

    if module_docker_extra_params:
        _validate_docker_extra_params(module_docker_extra_params)
        module_docker_extra_params = _sanitize_docker_extra_params(module_docker_extra_params)

    mount = docker_stark_container_mount if use_stark_container_mount else ""
    docker_parameters = f"--rm {docker_name} {module_docker_extra_params} {mount}"

    # Analysis file paths
    analysis_folder = docker_stark_api_log_folder
    analysis_file = os.path.join(analysis_folder, f"{analysis_id_name}.json")
    analysis_info_file = os.path.join(analysis_folder, f"{analysis_id_name}.info")
    analysis_output_file = os.path.join(analysis_folder, f"{analysis_id_name}.output")

    # Task spooler configuration
    _ts_env, _max_slots = _prepare_queue_for_submission(json_input.get("queue"))
    _task_slots = _resolve_task_slots(json_input, _max_slots)

    threads = (_task_slots == 0) and _max_slots or _task_slots
    json_input["threads"] = threads

    docker_parameters += f" --cpus={threads} "

    if "memory" in json_input and json_input.get("memory", None):
        memory = str(json_input.get("memory", "")).strip()
        if not re.fullmatch(r"\d+(?:[bBkKmMgG])?", memory):
            raise ValueError(f"Invalid memory value: {memory}")
        docker_parameters += f" --memory={memory} "

    prioritize = json_input.get("prioritize", False)

    # Write the full JSON (with all fields including command) for traceability
    with open(analysis_file, "w") as f:
        f.write(json.dumps(json_input))

    # Build the docker run command
    docker_cmd = ["docker", "run"]
    docker_cmd += _safe_split(docker_parameters)
    docker_cmd.append(image)
    docker_cmd_str = " ".join(shlex.quote(tok) for tok in docker_cmd)

    ts_cmd = f"{_ts_env}{ts} -N {_task_slots} -L {analysis_id_name}" if ts else ""

    raw_command = str(json_input.get("command", "")).strip()
    _validate_command(raw_command)
    command_tokens = _safe_split(raw_command)
    command_str = " ".join(shlex.quote(t) for t in command_tokens)

    my_cmd = (
        f'{ts_cmd} bash -c "'
        f"trap 'docker stop {analysis_id_name} 2>/dev/null; echo failed > {analysis_info_file}' TERM INT; "
        f"{docker_cmd_str} {command_str} "
        f"> {analysis_output_file} 2>&1 "
        f"&& (echo 'finished' > {analysis_info_file} && exit 0) "
        f"|| (echo 'failed' > {analysis_info_file} && exit 1)\""
    )

    _ = _queue_task(my_cmd, prioritize=prioritize, ts_cmd=f"{_ts_env}{ts}")

    return analysis_id_name


def queue_command_docker(json_input: dict) -> str:
    """Queue a docker command with predefined mount parameters (json key 'command_docker').
    Returns analysis_id_name or raises RuntimeError."""

    # Command
    command_docker = json_input["command_docker"]

    # Image
    image = json_input.get("image")
    if not image:
        raise ValueError("'image' is required when using 'command_docker'")
    _validate_image(image)

    # Extra docker parameters
    docker_extra_params = json_input.get("docker_extra_params", "")
    if docker_extra_params:
        _validate_docker_extra_params(docker_extra_params)
        docker_extra_params = _sanitize_docker_extra_params(docker_extra_params)

    # Build analysis_id_name
    analysis_id_name, _ = _build_analysis_idname(json_input)

    # Docker parameters: use analysis_id_name as Docker container name to be able to identify the container corresponding to a task-spooler task, and stop it if needed when the task is killed.
    docker_name = f"--name {analysis_id_name}"
    docker_parameters = f"--rm {docker_name} {docker_extra_params}"
    if json_input.get("use_stark_container_mount", False):
        docker_parameters += f" {docker_stark_container_mount}"

    # Analysis file paths
    analysis_folder = docker_stark_api_log_folder
    analysis_file = os.path.join(analysis_folder, f"{analysis_id_name}.json")
    analysis_info_file = os.path.join(analysis_folder, f"{analysis_id_name}.info")
    analysis_output_file = os.path.join(analysis_folder, f"{analysis_id_name}.output")

    # Task spooler configuration
    _ts_env, _max_slots = _prepare_queue_for_submission(json_input.get("queue"))
    _task_slots = _resolve_task_slots(json_input, _max_slots)

    # Threads
    threads = (_task_slots == 0) and _max_slots or _task_slots
    json_input["threads"] = threads

    # Write the JSON input for the container, which may be used for metrics and debugging.
    with open(analysis_file, "w") as f:
        f.write(json.dumps(json_input))

    # CPU/Threads
    docker_extra_params = _sanitize_docker_extra_params(
        docker_extra_params, extra_forbidden={"--cpus": True}
    )
    docker_parameters += f" --cpus={threads} "

    # Memory
    if "memory" in json_input and json_input.get("memory", None):
        docker_extra_params = _sanitize_docker_extra_params(
            docker_extra_params, extra_forbidden={"--memory": True, "-m": True}
        )
        memory = str(json_input.get("memory", "")).strip()
        if not re.fullmatch(r"\d+(?:[bBkKmMgG])?", memory):
            raise ValueError(f"Invalid memory value: {memory}")
        docker_parameters += f" --memory={memory} "

    # Security checks for docker-compose command: we allow more freedom than for raw docker commands, but we still want to block some obviously dangerous patterns.
    docker_cmd = [
        "docker",
        "run",
    ]
    docker_cmd += _safe_split(docker_parameters)
    docker_cmd.append(image)
    docker_cmd += _safe_split(command_docker)
    docker_cmd_str = " ".join(shlex.quote(tok) for tok in docker_cmd)

    # Task spooler command
    ts_cmd = f"{_ts_env}{ts} -N {_task_slots} -L {analysis_id_name}" if ts else ""

    # Task spooler command: we wrap the docker run in a bash command that traps termination signals to stop the container and mark the task as failed, and writes "finished" or "failed" to an info file based on the exit status of the docker command. This allows us to track the status of the task and ensure cleanup if it's killed.
    my_cmd = (
        f'{ts_cmd} bash -c "'
        f"trap 'docker stop {analysis_id_name} 2>/dev/null; echo failed > {analysis_info_file}' TERM INT; "
        f"{docker_cmd_str} "
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
        docker_extra_params = _sanitize_docker_extra_params(docker_extra_params)

    # Analysis name
    analysis_id_name, _ = _build_analysis_idname(json_input)

    # Docker parameters: use analysis_id_name as Docker container name to be able to identify the container corresponding to a task-spooler task, and stop it if needed when the task is killed.
    docker_name = f"--name {analysis_id_name}"
    docker_parameters = f"--rm {docker_name} {docker_extra_params}"
    if json_input.get("use_stark_container_mount", False):
        docker_parameters += f" {docker_stark_container_mount}"

    # Analysis file paths
    analysis_folder = docker_stark_api_log_folder
    analysis_file = os.path.join(analysis_folder, f"{analysis_id_name}.json")
    analysis_info_file = os.path.join(analysis_folder, f"{analysis_id_name}.info")
    analysis_output_file = os.path.join(analysis_folder, f"{analysis_id_name}.output")

    # Touch output
    with open(analysis_output_file, "w") as f:
        f.write("")

    # Task spooler configuration
    _ts_env, _max_slots = _prepare_queue_for_submission(json_input.get("queue"))
    _task_slots = _resolve_task_slots(json_input, _max_slots)

    # Threads
    threads = (_task_slots == 0) and _max_slots or _task_slots
    json_input["threads"] = threads

    # Write the JSON input for the container, which may be used for metrics and debugging.
    with open(analysis_file, "w") as f:
        f.write(json.dumps(json_input))

    # Task spooler command

    # Security checks for docker-compose command: we allow more freedom than for raw docker commands, but we still want to block some obviously dangerous patterns.
    docker_cmd = [
        "docker-compose",
        "-f",
        docker_compose_file,
        "run",
    ]
    docker_cmd += _safe_split(docker_parameters)
    docker_cmd.append(service)
    docker_cmd += _safe_split(command_docker_compose)
    docker_cmd_str = " ".join(shlex.quote(tok) for tok in docker_cmd)

    # Task spooler command
    ts_cmd = f"{_ts_env}{ts} -N {_task_slots} -L {analysis_id_name}" if ts else ""

    # Task spooler command: we wrap the docker-compose run in a bash command that traps termination signals to stop the container and mark the task as failed, and writes "finished" or "failed" to an info file based on the exit status of the docker command. This allows us to track the status of the task and ensure cleanup if it's killed.
    my_cmd = (
        f'{ts_cmd} bash -c "'
        f"trap 'docker stop {analysis_id_name} 2>/dev/null; echo failed > {analysis_info_file}' TERM INT; "
        f"{docker_cmd_str} "
        f"> {analysis_output_file} 2>&1 "
        f"&& (echo 'finished' > {analysis_info_file} && exit 0) "
        f"|| (echo 'failed' > {analysis_info_file} && exit 1)\""
    )

    # Queue the task and get the task ID from task-spooler
    _ = _queue_task(
        my_cmd, prioritize=json_input.get("prioritize", False), ts_cmd=f"{_ts_env}{ts}"
    )

    return analysis_id_name
