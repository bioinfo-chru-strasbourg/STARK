#!/usr/bin/env python

import json
import os
import random
import re
import shlex
import string
import subprocess
from dataclasses import dataclass
from typing import Optional

from config import (
    docker_stark,
    docker_stark_api_log_folder,
    docker_stark_container_mount,
    ts,
    ts_timeout,
)
from queues import _prepare_queue_for_submission
from security import (
    _sanitize_analysis_name,
    _validate_command,
    _validate_container_name,
    _validate_docker_extra_params,
    _sanitize_docker_extra_params,
    _validate_image,
    _safe_split,
    ALLOWED_CLIENT_KEYS,
    SENSITIVE_KEYS
)
from modules import resolve_module, apply_module_defaults


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
        # run_id = _sanitize_analysis_name(str(json_input["run"].split(":")[0]))
        run_id = _sanitize_analysis_name(
            str(json_input["run"].split(":")[0].split("/")[-1])
        )
    elif "command" in json_input:
        if isinstance(json_input["command"], str):
            # Extract run name from the command string if possible, using the parameter "--run=<run_id> " (e.g. "--run=MY_RUN" or "--run MY_RUN")
            analysis_name_match = re.search(
                r"--analysis_name[=\s]+(\S+)", json_input["command"]
            )
            run_id_match = re.search(r"--run[=\s]+(\S+)", json_input["command"])
            if analysis_name_match:
                run_id = _sanitize_analysis_name(analysis_name_match.group(1))
            elif run_id_match:
                # run_id = _sanitize_analysis_name(run_id_match.group(1))
                run_id = _sanitize_analysis_name(
                    str(run_id_match.group(1).split(":")[0].split("/")[-1])
                )
            else:
                run_id = _sanitize_analysis_name(
                    str(json_input["command"].split(":")[0])
                )
        elif isinstance(json_input["command"], dict) and len(json_input["command"]) > 0:
            # Extract run name from the command dict if possible, using the key "run" (e.g. {"run": "MY_RUN"})
            analysis_name = json_input["command"].get("analysis_name")
            run = json_input["command"].get("run")
            if analysis_name:
                run_id = _sanitize_analysis_name(str(analysis_name))
            elif run:
                run_id = _sanitize_analysis_name(str(run.split(":")[0].split("/")[-1]))
            else:
                run_id = _sanitize_analysis_name(
                    str(
                        json_input["command"]
                        .get("run", "UNKNOWN")
                        .split(":")[0]
                        .split("/")[-1]
                    )
                )
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
                timeout=ts_timeout,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                f"timed out after {ts_timeout} seconds while reprioritizing task {task_id} "
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


# ---------------------------------------------------------------------------
# Private helpers shared by all queue_* functions
# ---------------------------------------------------------------------------


@dataclass
class _TaskContext:
    """Internal value object produced by _setup_task_context."""

    analysis_id_name: str
    analyses_run_name: str
    analysis_folder: str
    analysis_file: str
    analysis_info_file: str
    analysis_output_file: str
    ts_env: str
    max_slots: int
    task_slots: int
    threads: int


def _setup_task_context(json_input: dict) -> "_TaskContext":
    """Build analysis ID, file paths, and task-spooler context from json_input.

    Side-effect: sets ``json_input["threads"]`` to the resolved thread count.
    """
    analysis_id_name, analyses_run_name = _build_analysis_idname(json_input)
    analysis_folder = docker_stark_api_log_folder
    analysis_file = os.path.join(analysis_folder, f"{analysis_id_name}.json")
    analysis_info_file = os.path.join(analysis_folder, f"{analysis_id_name}.info")
    analysis_output_file = os.path.join(analysis_folder, f"{analysis_id_name}.output")
    _ts_env, _max_slots = _prepare_queue_for_submission(json_input.get("queue"))
    _task_slots = _resolve_task_slots(json_input, _max_slots)
    threads = (_task_slots == 0) and _max_slots or _task_slots
    json_input["threads"] = threads
    return _TaskContext(
        analysis_id_name=analysis_id_name,
        analyses_run_name=analyses_run_name,
        analysis_folder=analysis_folder,
        analysis_file=analysis_file,
        analysis_info_file=analysis_info_file,
        analysis_output_file=analysis_output_file,
        ts_env=_ts_env,
        max_slots=_max_slots,
        task_slots=_task_slots,
        threads=threads,
    )


def _build_docker_resource_flags(threads: int, json_input: dict) -> str:
    """Return ' --cpus=N [--memory=M] ' resource flag string for docker parameters."""
    flags = f" --cpus={threads} "
    if json_input.get("memory"):
        memory = str(json_input["memory"]).strip()
        if not re.fullmatch(r"\d+(?:[bBkKmMgG])?", memory):
            raise ValueError(f"Invalid memory value: {memory}")
        flags += f" --memory={memory} "
    return flags


def _build_ts_bash_wrapper(
    ts_cmd: str,
    inner_cmd_str: str,
    analysis_id_name: str,
    analysis_info_file: str,
    analysis_output_file: str,
) -> str:
    """Wrap inner_cmd_str in a bash trap for signal handling and status reporting.

    The running container is stopped via ``docker stop`` on TERM/INT.
    Use ``_build_ts_bash_wrapper_exec`` when the command runs inside an existing
    container (docker exec) where stopping the container is not appropriate.
    """
    return (
        f'{ts_cmd} bash -c "'
        f"trap 'docker stop {analysis_id_name} 2>/dev/null; echo failed > {analysis_info_file}' TERM INT; "
        f"{inner_cmd_str} "
        f"> {analysis_output_file} 2>&1 "
        f"&& (echo 'finished' > {analysis_info_file} && exit 0) "
        f"|| (echo 'failed' > {analysis_info_file} && exit 1)\""
    )


def _build_ts_bash_wrapper_exec(
    ts_cmd: str,
    inner_cmd_str: str,
    analysis_info_file: str,
    analysis_output_file: str,
) -> str:
    """Wrap a ``docker exec`` command in a bash trap for signal handling.

    Unlike ``_build_ts_bash_wrapper``, the signal handler kills only the exec
    process (by PID) rather than stopping the target container, which is a
    long-running service that must not be stopped.
    """
    return (
        f'{ts_cmd} bash -c "'
        f"{inner_cmd_str} > {analysis_output_file} 2>&1 & PID=\\$! ; "
        f"trap 'kill \\$PID 2>/dev/null; echo killed > {analysis_info_file}; exit 1' TERM INT; "
        f"wait \\$PID; RC=\\$?; "
        f"if [ \\$RC -eq 0 ]; then echo finished > {analysis_info_file}; else echo failed > {analysis_info_file}; fi; "
        f'exit \\$RC"'
    )


def queue_module_analysis(json_input: dict) -> str:
    """Queue a module analysis by passing the ``command`` field as direct CLI args.
    """

    # Resolve module
    module_name = json_input.get("module", "UNKNOWN")
    module_cfg = resolve_module(json_input)
    apply_module_defaults(json_input, module_cfg)

    # Filter with allowed keys for security, the rest of the config will be taken from module config to ensure security of sensitive parameters like "image" and "docker_extra_params" that could be misconfigured in the JSON input if we let them come from there, even if the module config is safe. The "image" parameter is especially important to have a safe default as it defines the Docker image used for analyses and we don't want it to be accidentally misconfigured to an unsafe value. The "docker_extra_params" and "use_stark_container_mount" parameters are also important to enforce from the module config to prevent potential abuse via unsafe extra Docker parameters or mounts if they are exposed in the JSON input.
    json_input = {k: v for k, v in json_input.items() if k in ALLOWED_CLIENT_KEYS}

    # Systematically fix sensitive parameters from module config to avoid security issues, even if the module config is misconfigured with unsafe values. This ensures that the API will enforce safe values for these critical parameters regardless of the module configuration, which is important for security since some of these parameters (e.g. "image") can have a big impact on the security of the system if set to an unsafe value. The "image" parameter is especially important to have a safe default as it defines the Docker image used for analyses and we don't want it to be accidentally misconfigured to an unsafe value. The "docker_extra_params" and "use_stark_container_mount" parameters are also important to enforce from the module config to prevent potential abuse via unsafe extra Docker parameters or mounts if they are exposed in the JSON input.
    for k in SENSITIVE_KEYS:
        if k in module_cfg:
            json_input[k] = module_cfg[k]

    # Add defaults if not in json_input, but only for the standard JSON input fields, not for the extra module config fields to allow flexibility in module configuration without affecting the expected structure of the JSON input. The module config may contain arbitrary extra fields that are not part of the standard JSON input but can be used for other purposes (e.g. defining extra Docker parameters or mounts) without affecting the expected structure of the JSON input that the API processes for launching analyses. The "defaults" key in the module config is specifically designed to allow setting default values for standard JSON input fields if they are not provided by the client, while still allowing the module config to contain extra fields for other purposes without affecting the JSON input structure.
    if "defaults" in module_cfg:
        for k, v in module_cfg["defaults"].items():
            if k not in json_input:
                json_input[k] = v

    # Add module name in JSON input for reference, even if it's not used directly in the analysis execution since the module config is already applied to set the relevant parameters for execution. This is just for reference and debugging purposes to keep track of which module was resolved for this analysis.
    json_input["module"] = module_name

    # The module config may contain the following keys (in addition to the standard JSON input fields)
    if "container" in json_input:
        return queue_command_docker_exec(json_input)    # docker exec
    elif "service" in json_input:
        return queue_command_docker_compose(json_input) # docker-compose run
    elif "image" in json_input:
        return queue_command_docker(json_input)         # docker run (client image)
    else:
        return queue_analysis(json_input)               # STARK JSON analysis


def queue_analysis(json_input: dict) -> str:
    """Build and queue a STARK Docker analysis. Returns analysisIDNAME or raises RuntimeError."""

    # Resolve module
    module_cfg = resolve_module(json_input)
    apply_module_defaults(json_input, module_cfg)
    image = module_cfg.get("image") or docker_stark
    module_docker_extra_params = module_cfg.get("docker_extra_params", "")
    use_stark_container_mount = module_cfg.get("use_stark_container_mount", False)

    # Task context (also sets json_input["threads"])
    ctx = _setup_task_context(json_input)

    # Docker parameters
    docker_name = f" --name {ctx.analysis_id_name} "
    if module_docker_extra_params:
        _validate_docker_extra_params(module_docker_extra_params)
        module_docker_extra_params = _sanitize_docker_extra_params(
            module_docker_extra_params
        )
    mount = docker_stark_container_mount if use_stark_container_mount else ""
    docker_parameters = (
        f"--rm {docker_name} {module_docker_extra_params} {mount}"
        + _build_docker_resource_flags(ctx.threads, json_input)
    )

    prioritize = json_input.get("prioritize", False)

    # Write the JSON
    with open(ctx.analysis_file, "w") as f:
        f.write(json.dumps(json_input))

    # Write the container-specific JSON (stripped of host-side scheduling params)
    analysis_file_stark = os.path.join(
        ctx.analysis_folder, f"{ctx.analysis_id_name}.json.stark_analysis"
    )
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
    with open(analysis_file_stark, "w") as f:
        f.write(json.dumps(json_input_for_container))

    # Build docker command
    docker_cmd = ["docker", "run"] + _safe_split(docker_parameters) + [image]
    docker_cmd_str = " ".join(shlex.quote(tok) for tok in docker_cmd)
    inner_cmd = (
        f"{docker_cmd_str} "
        f"--analysis_name={ctx.analyses_run_name} --analysis={analysis_file_stark}"
    )

    ts_cmd = (
        f"{ctx.ts_env}{ts} -N {ctx.task_slots} -L {ctx.analysis_id_name}" if ts else ""
    )
    my_cmd = _build_ts_bash_wrapper(
        ts_cmd,
        inner_cmd,
        ctx.analysis_id_name,
        ctx.analysis_info_file,
        ctx.analysis_output_file,
    )

    _ = _queue_task(my_cmd, prioritize=prioritize, ts_cmd=f"{ctx.ts_env}{ts}")
    return ctx.analysis_id_name


def queue_command_docker_exec(json_input: dict) -> str:
    """Queue a command inside an already-running container via ``docker exec``.

    Returns analysis_id_name or raises RuntimeError.

    The command is read from the ``command`` key of json_input.
    Unlike ``queue_command_docker``, no new container is created. The target
    container must already be running. On TERM/INT the exec process is killed
    by PID; the container itself is left untouched.

    ``docker exec`` does not support ``--cpus`` or ``--memory``; those flags
    are ignored. Thread count only affects the task-spooler slot count (-N).

    JSON keys:
        command             (str, required)  : command to run inside the container.
        container           (str, required)  : name or ID of the running container.
        docker_extra_params (str, optional)  : exec flags, e.g. ``-e VAR=val``,
                                               ``-w /workdir``, ``-u user``.
        analysis_name       (str, optional)  : human-readable task label.
        queue               (str, optional)  : target queue.
        threads             (int, optional)  : task-spooler slot count.
        prioritize          (bool, optional) : move to front of the queue.
    """

    # Command
    raw_command = str(json_input["command"]).strip()
    _validate_command(raw_command)

    # Container (name or ID of an already-running container)
    container = json_input.get("container")
    if not container:
        raise ValueError("'container' is required for docker exec mode")
    _validate_container_name(container)

    # Extra docker exec parameters (client-supplied; e.g. -e VAR=val, -w /path, -u user)
    docker_extra_params = json_input.get("docker_extra_params", "")
    if docker_extra_params:
        _validate_docker_extra_params(docker_extra_params)
        docker_extra_params = _sanitize_docker_extra_params(docker_extra_params)

    # Extra command prefix
    command_prefix = json_input.get("command_prefix", "")
    if command_prefix:
        _validate_command(command_prefix)
        command_prefix = command_prefix.strip() + " " # add space after prefix

    # Extra command postfix
    command_postfix = json_input.get("command_postfix", "")
    if command_postfix:
        _validate_command(command_postfix)
        command_postfix = " " + command_postfix.strip() # add space before postfix

    # Task context (also sets json_input["threads"])
    ctx = _setup_task_context(json_input)

    prioritize = json_input.get("prioritize", False)

    # Write the JSON input for traceability
    with open(ctx.analysis_file, "w") as f:
        f.write(json.dumps(json_input))

    # Build docker exec command: docker exec [PARAMS] CONTAINER COMMAND...
    docker_cmd = ["docker", "exec"]
    if docker_extra_params:
        docker_cmd += _safe_split(docker_extra_params)
    docker_cmd.append(container)
    docker_cmd += _safe_split(command_prefix)
    docker_cmd += _safe_split(raw_command)
    docker_cmd += _safe_split(command_postfix)
    docker_cmd_str = " ".join(shlex.quote(tok) for tok in docker_cmd)

    ts_cmd = f"{ctx.ts_env}{ts} -N {ctx.task_slots} -L {ctx.analysis_id_name}" if ts else ""
    my_cmd = _build_ts_bash_wrapper_exec(
        ts_cmd, docker_cmd_str, ctx.analysis_info_file, ctx.analysis_output_file
    )

    _ = _queue_task(my_cmd, prioritize=prioritize, ts_cmd=f"{ctx.ts_env}{ts}")
    return ctx.analysis_id_name


def queue_command_docker(json_input: dict) -> str:
    """Queue a command in a new ephemeral Docker container (``docker run``).

    The command is read from the ``command`` key of json_input.
    Returns analysis_id_name or raises RuntimeError.

    JSON keys:
        command             (str, required)  : command to run inside the container.
        image               (str, required)  : Docker image to use.
        docker_extra_params (str, optional)  : extra ``docker run`` flags.
        use_stark_container_mount (bool, optional): mount STARK volumes.
        analysis_name       (str, optional)  : human-readable task label.
        queue               (str, optional)  : target queue.
        threads             (int, optional)  : task-spooler slot count / --cpus.
        memory              (str, optional)  : Docker memory limit, e.g. ``4G``.
        prioritize          (bool, optional) : move to front of the queue.
    """

    # Command
    command = json_input["command"]

    # Image
    image = json_input.get("image")
    if not image:
        raise ValueError("'image' is required for docker run mode")
    _validate_image(image)

    # Extra docker parameters (client-supplied; validated before embedding)
    docker_extra_params = json_input.get("docker_extra_params", "")
    if docker_extra_params:
        _validate_docker_extra_params(docker_extra_params)
        docker_extra_params = _sanitize_docker_extra_params(docker_extra_params)

    # Extra command prefix
    command_prefix = json_input.get("command_prefix", "")
    if command_prefix:
        _validate_command(command_prefix)
        command_prefix = command_prefix.strip() + " "  # add space after prefix

    # Extra command postfix
    command_postfix = json_input.get("command_postfix", "")
    if command_postfix:
        _validate_command(command_postfix)
        command_postfix = " " + command_postfix.strip()  # add space before postfix

    # Task context (also sets json_input["threads"])
    ctx = _setup_task_context(json_input)

    # Docker parameters
    docker_name = f"--name {ctx.analysis_id_name}"
    docker_parameters = f"--rm {docker_name} {docker_extra_params}"
    if json_input.get("use_stark_container_mount", False):
        docker_parameters += f" {docker_stark_container_mount}"

    prioritize = json_input.get("prioritize", False)

    # Write the JSON input for the container
    with open(ctx.analysis_file, "w") as f:
        f.write(json.dumps(json_input))

    docker_parameters += _build_docker_resource_flags(ctx.threads, json_input)

    # Build docker command
    docker_cmd = ["docker", "run"]
    docker_cmd += _safe_split(docker_parameters)
    docker_cmd += [image] 
    docker_cmd += _safe_split(command_prefix) 
    docker_cmd += _safe_split(command) 
    docker_cmd += _safe_split(command_postfix)
    docker_cmd_str = " ".join(shlex.quote(tok) for tok in docker_cmd)

    ts_cmd = f"{ctx.ts_env}{ts} -N {ctx.task_slots} -L {ctx.analysis_id_name}" if ts else ""
    my_cmd = _build_ts_bash_wrapper(
        ts_cmd, docker_cmd_str, ctx.analysis_id_name, ctx.analysis_info_file, ctx.analysis_output_file
    )

    _ = _queue_task(my_cmd, prioritize=prioritize, ts_cmd=f"{ctx.ts_env}{ts}")
    return ctx.analysis_id_name


def queue_command_docker_compose(json_input: dict) -> str:
    """Queue a command via ``docker-compose run``.

    The command is read from the ``command`` key of json_input.
    Returns analysis_id_name or raises RuntimeError.

    JSON keys:
        command             (str, required)  : command to run in the service container.
        service             (str, required)  : docker-compose service name.
        docker_compose_file (str, required)  : path to the docker-compose YAML file.
        docker_extra_params (str, optional)  : extra ``docker-compose run`` flags.
        use_stark_container_mount (bool, optional): mount STARK volumes.
        analysis_name       (str, optional)  : human-readable task label.
        queue               (str, optional)  : target queue.
        threads             (int, optional)  : task-spooler slot count.
        prioritize          (bool, optional) : move to front of the queue.
    """

    # Command
    command = json_input["command"]

    # Service
    service = json_input.get("service")
    if not service:
        raise ValueError("'service' is required for docker-compose mode")

    # Docker Compose file
    docker_compose_file = json_input.get("docker_compose_file")
    if not docker_compose_file:
        raise ValueError(
            "'docker_compose_file' is required when using 'command_docker_compose'"
        )
    if not os.path.isfile(docker_compose_file):
        raise ValueError(f"'docker_compose_file' does not exist: {docker_compose_file}")

    # Extra docker parameters (client-supplied; validated before embedding)
    docker_extra_params = json_input.get("docker_extra_params", "")
    if docker_extra_params:
        _validate_docker_extra_params(docker_extra_params)
        docker_extra_params = _sanitize_docker_extra_params(docker_extra_params)

    # Extra command prefix
    command_prefix = json_input.get("command_prefix", "")
    if command_prefix:
        _validate_command(command_prefix)
        command_prefix = command_prefix.strip() + " "  # add space after prefix

    # Extra command postfix
    command_postfix = json_input.get("command_postfix", "")
    if command_postfix:
        _validate_command(command_postfix)
        command_postfix = " " + command_postfix.strip()  # add space before postfix

    # Task context (also sets json_input["threads"])
    ctx = _setup_task_context(json_input)

    # Docker parameters
    docker_name = f"--name {ctx.analysis_id_name}"
    docker_parameters = f"--rm {docker_name} {docker_extra_params}"
    if json_input.get("use_stark_container_mount", False):
        docker_parameters += f" {docker_stark_container_mount}"

    # Prioritize
    prioritize = json_input.get("prioritize", False)

    # Touch output
    with open(ctx.analysis_output_file, "w") as f:
        f.write("")

    # Write the JSON input
    with open(ctx.analysis_file, "w") as f:
        f.write(json.dumps(json_input))

    # Build docker-compose command
    docker_cmd = ["docker-compose", "-f", docker_compose_file, "run"]
    docker_cmd += _safe_split(docker_parameters)
    docker_cmd.append(service)
    docker_cmd += _safe_split(command_prefix)
    docker_cmd += _safe_split(command)
    docker_cmd += _safe_split(command_postfix)
    docker_cmd_str = " ".join(shlex.quote(tok) for tok in docker_cmd)

    ts_cmd = f"{ctx.ts_env}{ts} -N {ctx.task_slots} -L {ctx.analysis_id_name}" if ts else ""
    my_cmd = _build_ts_bash_wrapper(
        ts_cmd, docker_cmd_str, ctx.analysis_id_name, ctx.analysis_info_file, ctx.analysis_output_file
    )

    _ = _queue_task(my_cmd, prioritize=prioritize, ts_cmd=f"{ctx.ts_env}{ts}")
    return ctx.analysis_id_name
