#!/usr/bin/env python

from fastapi import APIRouter  # pyright: ignore[reportMissingImports]

from peers import get_local_metrics, get_self_hostname, get_self_url, load_peers

router = APIRouter()


@router.get("/whoami")
def whoami():
    """Return this node's hostname identifier."""
    return {"id": get_self_hostname()}


@router.get("/metrics")
def metrics():
    """Return slot metrics for all queues on this node.

    Response format:
    {
        "node":     "container-abc123",
        "self_url": "http://192.168.1.10:8001",   // null if unknown
        "queues": {
            "stark":  {"configured": 1, "running": 0, "queued": 0},
            "light":  {"configured": 4, "running": 2, "queued": 1}
        }
    }
    """
    return {
        "node": get_self_hostname(),
        "self_url": get_self_url(),
        "queues": get_local_metrics(),
    }


@router.get("/peers")
def peers():
    """Return the list of configured peers."""
    return {"peers": load_peers()}
