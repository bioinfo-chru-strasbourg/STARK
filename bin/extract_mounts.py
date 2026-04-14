#!/usr/bin/env python
import argparse
import json
import subprocess
import sys


def run_cmd(cmd):
    """Run a shell command and return stdout or exit on failure."""
    try:
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        return output.decode().strip()
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {' '.join(cmd)}")
        print(e.output.decode())
        sys.exit(1)


def get_container_id():
    """Read the hostname which contains the short container ID."""
    with open("/etc/hostname") as f:
        return f.read().strip()


def main(docker):

    # Container ID
    container_id = get_container_id()
    # print(f"[INFO] Detected container hostname (short ID): {container_id}")

    # Inspect JSON
    try:
        raw_json = run_cmd([docker, "inspect", container_id])
        data = json.loads(raw_json)[0]
    except Exception as e:
        print("[ERROR] Failed to parse docker inspect JSON:", e)
        sys.exit(1)

    # print("\n=== Mounts detected ===")

    mounts_array = []
    mounts_array_for_check = []

    # 1) MOUNTS (runtime mounts)
    mounts = data.get("Mounts", [])
    if mounts:
        # print("\n# Mounts[]:")
        for m in mounts:
            src = m.get("Source")
            dst = m.get("Destination")
            mode = m.get("Mode", "")
            rw = m.get("RW", True)

            opt = mode if mode else ("rw" if rw else "ro")

            mount = f"{src}:{dst}:{opt}"
            #print(mount.split(":")[0:2])
            if mount.split(":")[0:2] not in mounts_array_for_check:
                # print(f"-v {mount}")
                mounts_array.append(mount)
                mounts_array_for_check.append(mount.split(":")[0:2])
            # print(f"-v {src}:{dst}:{opt}")
    # else:
    #     print("\n# No Mounts[] found")

    # 2) HostConfig.binds
    binds = data.get("HostConfig", {}).get("Binds", [])
    if binds:
        # print("\n# HostConfig.Binds:")
        for b in binds:
            mount = f"{b}"
            #print(mount.split(":")[0:2])
            if mount.split(":")[0:2] not in mounts_array_for_check:
                # print(f"-v {mount}")
                mounts_array.append(mount)
                mounts_array_for_check.append(mount.split(":")[0:2])
    # else:
    #     print("\n# No HostConfig.Binds found")

    if not mounts_array:
        print("")
    else:
        print(" -v " + " -v ".join(mounts_array))

    # print("\n=== End ===")


if __name__ == "__main__":

    # Arguments
    # add parameter for the docker binary
    parser = argparse.ArgumentParser(description="Extract mounts from within a Docker container")
    parser.add_argument("--docker", default="docker", help="Path to the docker binary (default: docker)")
    args = parser.parse_args()

    # catch arguments
    # docker
    docker = args.docker
    # check if docker is available
    try:
        run_cmd([docker, "--version"])
    except Exception as e:
        print(f"[ERROR] Docker binary '{docker}' is not available or not working: {e}")
        sys.exit(1)

    main(docker=docker)
