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


def main(docker, volumes_from_enable=False) -> str:

    result = ""

    # Container ID
    container_id = get_container_id()

    # If volumes_from_enable is True, print the --volumes-from option and exit
    if volumes_from_enable:
        result = f" --volumes-from {container_id}"

    # Retrieve the mounts from the container using docker inspect
    else:

        # Inspect JSON
        try:
            raw_json = run_cmd([docker, "inspect", container_id])
            data = json.loads(raw_json)[0]
        except Exception as e:
            print("[ERROR] Failed to parse docker inspect JSON:", e)
            sys.exit(1)

        # Init
        mounts_array = []
        mounts_array_for_check = []

        # 1) MOUNTS (runtime mounts)
        mounts = data.get("Mounts", [])
        if mounts:

            # For each mount, construct the mount string and add it to the mounts_array if it's not already present
            for m in mounts:

                # Extract source, destination, mode, and rw from the mount dictionary
                src = m.get("Source")
                dst = m.get("Destination")
                mode = m.get("Mode", "")
                rw = m.get("RW", True)

                # Determine the mount option based on mode and rw
                opt = mode if mode else ("rw" if rw else "ro")

                # Construct the mount string in the format "source:destination:options"
                mount = f"{src}:{dst}:{opt}"
                if mount.split(":")[0:2] not in mounts_array_for_check:
                    mounts_array.append(mount)
                    mounts_array_for_check.append(mount.split(":")[0:2])

        # 2) HostConfig.binds
        binds = data.get("HostConfig", {}).get("Binds", [])
        if binds:
            # For each bind, construct the mount string and add it to the mounts_array if it's not already present
            for b in binds:

                # Extract source, destination, and options from the bind string
                mount = f"{b}"

                # Check if the mount (source and destination) is already in mounts_array_for_check
                if mount.split(":")[0:2] not in mounts_array_for_check:
                    mounts_array.append(mount)
                    mounts_array_for_check.append(mount.split(":")[0:2])

        # Print the mounts in the required format
        if not mounts_array:
            print("")
        else:
            result = " -v " + " -v ".join(mounts_array)

    return result


if __name__ == "__main__":

    # Arguments
    # add parameter for the docker binary
    parser = argparse.ArgumentParser(description="Extract mounts from within a Docker container")
    parser.add_argument("--docker", default="docker", help="Path to the docker binary (default: docker)")
    parser.add_argument("--volumes_from_enable", action="store_true", help="Enable --volumes-from option (default: False)")
    args = parser.parse_args()

    # catch arguments
    # docker
    docker = args.docker
    # volumes_from_enable
    volumes_from_enable = args.volumes_from_enable

    # check if docker is available
    try:
        run_cmd([docker, "--version"])
    except Exception as e:
        print(f"[ERROR] Docker binary '{docker}' is not available or not working: {e}")
        sys.exit(1)

    output = main(docker=docker, volumes_from_enable=volumes_from_enable)
    print(output)
