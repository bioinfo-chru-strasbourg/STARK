#!/usr/bin/env python3
"""STARK task launcher.

Replaces the bash wrappers (_build_ts_bash_wrapper / _build_ts_bash_wrapper_exec)
used to submit jobs to task-spooler.

Runs a command, redirects stdout+stderr to a .output file,
writes "finished" / "failed" / "killed" to a .info file,
and handles SIGTERM/SIGINT appropriately per mode.

Modes:
  docker-run      docker run ... IMAGE COMMAND
                  Signal: docker stop <name>

  docker-exec     docker exec [params] CONTAINER COMMAND
                  Signal: kill <PID>  (container is left untouched)

  docker-compose  docker-compose -f FILE run ... SERVICE COMMAND
                  Signal: docker stop <name>

  endpoint        curl --config <file> --data @<payload>
                  Signal: kill <PID>

Usage:
  launch.py --mode docker-run      --name <id> --output <path> --info <path> --cmd-file <path>
  launch.py --mode docker-exec              --output <path> --info <path> --cmd-file <path>
  launch.py --mode docker-compose  --name <id> --output <path> --info <path> --cmd-file <path>
  launch.py --mode endpoint                 --output <path> --info <path> \\
                  --curl-config <path> --curl-data <path> [--remove]

Exit codes:
  0   command succeeded (+ no JSON-RPC error in endpoint mode)
  1   command failed / JSON-RPC error
  2   argument error
  N   curl exit code (endpoint mode, network/HTTP failure)
"""

import argparse
import json
import os
import signal
import subprocess
import sys


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_info(info_file: str, status: str) -> None:
    """Write status string to the .info file (best-effort)."""
    try:
        with open(info_file, "w") as f:
            f.write(status)
    except OSError:
        pass


# ---------------------------------------------------------------------------
# Mode implementations
# ---------------------------------------------------------------------------


def _run_with_docker_stop(cmd: list, name: str, output_file: str, info_file: str) -> int:
    """Run cmd; on SIGTERM/SIGINT, stop the named docker container.

    Equivalent to _build_ts_bash_wrapper.
    """
    proc = None

    def _on_signal(signum, frame):
        if proc is not None:
            subprocess.run(["docker", "stop", name], capture_output=True)
        _write_info(info_file, "failed")
        sys.exit(1)

    signal.signal(signal.SIGTERM, _on_signal)
    signal.signal(signal.SIGINT, _on_signal)

    with open(output_file, "w") as out:
        proc = subprocess.Popen(cmd, stdout=out, stderr=out)
        proc.wait()

    return proc.returncode


def _run_with_pid_kill(cmd: list, output_file: str, info_file: str) -> int:
    """Run cmd; on SIGTERM/SIGINT, kill the process by PID (container untouched).

    Equivalent to _build_ts_bash_wrapper_exec.
    """
    proc = None

    def _on_signal(signum, frame):
        if proc is not None:
            try:
                proc.kill()
            except ProcessLookupError:
                pass
        _write_info(info_file, "killed")
        sys.exit(1)

    signal.signal(signal.SIGTERM, _on_signal)
    signal.signal(signal.SIGINT, _on_signal)

    with open(output_file, "w") as out:
        proc = subprocess.Popen(cmd, stdout=out, stderr=out)
        proc.wait()

    return proc.returncode


def _run_endpoint(
    curl_config: str,
    curl_data: str,
    remove: bool,
    output_file: str,
    info_file: str,
) -> int:
    """POST via curl using a config file; check the JSON-RPC response for errors.

    Integrates the logic of curl.py with proper .output/.info lifecycle.
    Signal: kill the curl process by PID (no docker container involved).
    """
    proc = None

    def _on_signal(signum, frame):
        if proc is not None:
            try:
                proc.kill()
            except ProcessLookupError:
                pass
        _write_info(info_file, "killed")
        sys.exit(1)

    signal.signal(signal.SIGTERM, _on_signal)
    signal.signal(signal.SIGINT, _on_signal)

    cmd = [
        "curl", "--silent", "--show-error", "--fail",
        "--config", curl_config,
        "--data", f"@{curl_data}",
    ]

    with open(output_file, "w") as out:
        proc = subprocess.Popen(cmd, stdout=out, stderr=out)
        proc.wait()

    curl_rc = proc.returncode

    # Remove config / payload files if requested
    if remove:
        for path in (curl_config, curl_data):
            try:
                os.remove(path)
            except FileNotFoundError:
                pass

    # Network / HTTP error
    if curl_rc != 0:
        return curl_rc

    # JSON-RPC error: read back the output to inspect the response
    try:
        with open(output_file) as f:
            response = f.read()
        data = json.loads(response)
        rpc_rc = (data.get("error") or {}).get("code", 0) or 0
    except (OSError, json.JSONDecodeError, AttributeError):
        rpc_rc = 0

    return 1 if rpc_rc != 0 else 0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="STARK task launcher — wraps task-spooler jobs with output/info lifecycle management."
    )
    parser.add_argument(
        "--mode", "-mode",
        required=True,
        choices=["docker-run", "docker-exec", "docker-compose", "endpoint"],
        help="Launch mode",
    )
    parser.add_argument(
        "--output", "-output",
        required=True,
        metavar="PATH",
        help="Path to the .output file (stdout+stderr of the command)",
    )
    parser.add_argument(
        "--info", "-info",
        required=True,
        metavar="PATH",
        help="Path to the .info file (finished / failed / killed)",
    )
    parser.add_argument(
        "--name", "-name",
        default="",
        metavar="NAME",
        help="Container name — used for 'docker stop' on signal (docker-run / docker-compose)",
    )
    # Docker modes
    parser.add_argument(
        "--cmd-file", "-cmd-file",
        default=None,
        metavar="PATH",
        help="Path to a JSON file containing the command as an array, e.g. [\"docker\",\"run\",\"--rm\",...]'",
    )
    parser.add_argument(
        "--remove", "-remove",
        action="store_true",
        default=False,
        help="Remove cmd-file / curl config and payload files after the call",
    )
    # Endpoint mode
    parser.add_argument(
        "--curl-config", "-curl-config",
        default=None,
        metavar="PATH",
        help="Path to the curl config file (endpoint mode)",
    )
    parser.add_argument(
        "--curl-data", "-curl-data",
        default=None,
        metavar="PATH",
        help="Path to the JSON payload file (endpoint mode)",
    )

    args = parser.parse_args()

    # --- Dispatch ---

    if args.mode == "endpoint":
        if not args.curl_config or not args.curl_data:
            print(
                "Error: --curl-config and --curl-data are required for endpoint mode",
                file=sys.stderr,
            )
            sys.exit(2)
        rc = _run_endpoint(
            args.curl_config, args.curl_data, args.remove, args.output, args.info
        )

    else:
        if not args.cmd_file:
            print("Error: --cmd-file is required for docker modes", file=sys.stderr)
            sys.exit(2)
        try:
            with open(args.cmd_file) as f:
                cmd = json.load(f)
            if not isinstance(cmd, list) or not cmd:
                raise ValueError("cmd-file must contain a non-empty JSON array")
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            print(f"Error: invalid --cmd-file: {exc}", file=sys.stderr)
            sys.exit(2)

        if args.remove:
            try:
                os.remove(args.cmd_file)
            except FileNotFoundError:
                pass

        if args.mode in ("docker-run", "docker-compose"):
            rc = _run_with_docker_stop(cmd, args.name, args.output, args.info)
        else:  # docker-exec
            rc = _run_with_pid_kill(cmd, args.output, args.info)

    # Write final .info status and exit
    _write_info(args.info, "finished" if rc == 0 else "failed")
    sys.exit(rc)


if __name__ == "__main__":
    main()
