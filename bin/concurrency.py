#!/usr/bin/env python3
# Author refactor: Samuel Nicaise, Antony Le Béchec
# Refactor NFS-safe concurrency

import argparse
import errno
import glob
import os
import subprocess
import time
import random

# from os.path import join as osj


# -------------------------------
#  LOCKFILE UTILS (NFS-SAFE)
# -------------------------------


def make_unique_lockfile_path(prefix, target):
    """Create a unique lockfile name including PID to avoid collisions."""
    base = os.path.basename(target)
    return f"{prefix}{base}.pid{os.getpid()}"


def atomic_create_lock(lockfile_path):
    """
    Create lockfile atomically using O_CREAT|O_EXCL.
    Returns True if the file was created, False if it already exists.
    """
    try:
        fd = os.open(lockfile_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.close(fd)
        return True
    except OSError as e:
        if e.errno == errno.EEXIST:
            return False
        raise  # real error → rethrow


def delete_lock(lockfile_path):
    """Remove a lockfile if it exists."""
    try:
        os.remove(lockfile_path)
    except FileNotFoundError:
        pass


def count_current_jobs(lockfile_prefix):
    """Count number of lockfiles in the directory."""
    return len(glob.glob(lockfile_prefix + "*"))


# -------------------------------
#   CONCURRENCY CONTROL
# -------------------------------


def is_rule_startable(lockfile_prefix, target, max_jobs, shell_mode=False):
    """
    NFS-safe rule start control:
    1. Create own lockfile atomically
    2. Count lockfiles
    3. Rollback if > max_jobs
    """
    lockfile = make_unique_lockfile_path(lockfile_prefix, target)

    # 1. Create lockfile atomically
    created = atomic_create_lock(lockfile)
    if not created:
        # Target is already started by same PID or improbable race
        if shell_mode:
            print("0")
        return False

    # 2. Count total concurrency
    total_jobs = count_current_jobs(lockfile_prefix)

    if total_jobs > max_jobs:
        # 3. Rollback: remove OUR lock only
        delete_lock(lockfile)
        if shell_mode:
            print("0")
        return False

    # OK
    if shell_mode:
        print("1")
    return True


def rule_finished(lockfile_prefix, target):
    """
    Clean lockfiles belonging to the current PID and target.
    This prevents deleting the lockfile of another job.
    """
    pattern = f"{lockfile_prefix}{os.path.basename(target)}.pid{os.getpid()}"
    for lock in glob.glob(pattern + "*"):
        delete_lock(lock)


def randomized_sleep(target):
    """Spread random time based on name hash."""
    seed = sum(ord(x) for x in target)
    random.seed(seed)
    wait = random.random() * 15
    time.sleep(wait)


def launch_when_possible(
    cmd, lockfile_prefix, target, max_jobs, current_dir, timeout_hours
):
    """
    Main loop: retry until allowed by semaphore.
    """
    randomized_sleep(target)
    starting = time.time()
    timeout = timeout_hours * 3600

    while True:
        try:
            if is_rule_startable(lockfile_prefix, target, max_jobs):
                # Launch command
                print("Launching:", cmd)

                # Subprocess safely, without shell injection
                ret = subprocess.call(cmd, shell=True, cwd=current_dir)

                # Job finished
                rule_finished(lockfile_prefix, target)
                return ret

        finally:
            # Clean our own lock if any
            rule_finished(lockfile_prefix, target)

        if time.time() - starting > timeout:
            raise RuntimeError(f"Timeout: couldn't create target: {target}")

        time.sleep(10)


# -------------------------------
#   CLI
# -------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    parser_launch = subparsers.add_parser("launch")
    parser_startable = subparsers.add_parser("is_rule_startable")
    parser_finished = subparsers.add_parser("rule_finished")

    # Common args
    for p in (parser_launch, parser_startable, parser_finished):
        p.add_argument("-l", "--lockfile_prefix", required=True)
        p.add_argument("-t", "--target", required=True)

    # Additional args
    for p in (parser_launch, parser_startable):
        p.add_argument("-m", "--max_jobs", type=int, required=True)
        p.add_argument("--shell_mode", action="store_true", default=False)

    parser_launch.add_argument("-c", "--cmd", required=True)
    parser_launch.add_argument("-d", "--current_dir", default=os.getcwd())
    parser_launch.add_argument("-z", "--timeout", type=int, default=48)

    args = parser.parse_args()

    if args.subcommand == "launch":
        launch_when_possible(
            args.cmd,
            args.lockfile_prefix,
            args.target,
            args.max_jobs,
            args.current_dir,
            args.timeout
        )

    elif args.subcommand == "is_rule_startable":
        is_rule_startable(
            args.lockfile_prefix,
            args.target,
            args.max_jobs,
            args.shell_mode
        )

    elif args.subcommand == "rule_finished":
        rule_finished(
            args.lockfile_prefix,
            args.target
        )
