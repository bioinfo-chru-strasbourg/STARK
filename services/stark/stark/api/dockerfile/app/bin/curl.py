#!/usr/bin/env python3
"""RPC curl wrapper: executes a curl call from a config file and checks the JSON-RPC response.

Usage:
    curl.py --config <path> [--remove]

Exit codes:
    0   success (curl OK + no JSON-RPC error)
    1   JSON-RPC error (.error.code != 0)
    N   curl error (network / HTTP failure, N = curl exit code)
"""

import argparse
import json
import os
import subprocess
import sys


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Execute a curl call from a config file and check the JSON-RPC response."
    )
    parser.add_argument(
        "--config", "-config",
        required=True,
        metavar="PATH",
        help="Path to the curl config file",
    )
    parser.add_argument(
        "--data",
        "-data",
        required=True,
        metavar="PATH",
        help="Path to the curl data file",
    )
    parser.add_argument(
        "--remove", "-remove",
        action="store_true",
        default=False,
        help="Remove the config file after the call",
    )
    args = parser.parse_args()

    # 1. RPC call
    result = subprocess.run(
        ["curl", "--silent", "--show-error", "--fail", "--config", args.config, "--data", f"@{args.data}"],
        capture_output=True,
        text=True,
    )
    response = result.stdout
    curl_rc = result.returncode

    # 1.1 Remove file if requested
    if args.remove:
        try:
            os.remove(args.config)
            os.remove(args.data)
        except FileNotFoundError:
            pass

    # 1.2 Print response
    print(response, end="")

    # 2. Network/curl error
    if curl_rc != 0:
        print(f"Curl error: {curl_rc}", file=sys.stderr)
        sys.exit(curl_rc)

    # 3. RPC error
    try:
        data = json.loads(response)
        rpc_rc = (data.get("error") or {}).get("code", 0) or 0
    except (json.JSONDecodeError, AttributeError):
        rpc_rc = 0

    if rpc_rc != 0:
        print(f"RPC error: {rpc_rc}", file=sys.stderr)
        sys.exit(1)

    # 4. Success
    sys.exit(0)


if __name__ == "__main__":
    main()
