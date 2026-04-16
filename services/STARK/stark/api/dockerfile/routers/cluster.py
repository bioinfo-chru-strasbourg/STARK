#!/usr/bin/env python

import json as _json

from fastapi import APIRouter  # pyright: ignore[reportMissingImports]

from peers import (
    get_local_metrics,
    get_peers_metrics,
    get_peers_tasks,
    get_self_hostname,
    get_self_url,
    load_peers,
)

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


@router.get("/cluster/summary")
async def cluster_summary():
    """Aggregate queue slot metrics across all nodes (self + peers).

    Response format:
    {
        "nodes": {
            "node1": {
                "url":    "http://192.168.1.10:4200",
                "status": "online",
                "queues": {
                    "stark": {"configured": 1, "running": 1, "queued": 0, "available": 0},
                    "light": {"configured": 4, "running": 2, "queued": 1, "available": 1}
                }
            },
            "node2": {"url": "http://192.168.1.11:4200", "status": "offline"}
        },
        "totals": {
            "stark": {"configured": 2, "running": 1, "queued": 0, "available": 1},
            "light": {"configured": 8, "running": 2, "queued": 1, "available": 5}
        }
    }

    The local node is always included. Offline peers appear with status "offline"
    and no queue detail.
    """
    peers_metrics = await get_peers_metrics()  # {url: queues_dict}

    peer_name_map = {
        p["url"]: p.get("name", p["url"]) for p in load_peers() if p.get("url")
    }

    self_url = get_self_url()
    self_name = (
        peer_name_map.get(self_url, get_self_hostname())
        if self_url
        else get_self_hostname()
    )
    local_queues_raw = get_local_metrics()

    nodes: dict = {}

    # Self node
    nodes[self_name] = {
        "url": self_url,
        "status": "online",
        "queues": {
            q: {**v, "available": max(0, v["configured"] - v["running"] - v["queued"])}
            for q, v in local_queues_raw.items()
        },
    }

    # Peers
    responding_urls = set(peers_metrics.keys())
    for p in load_peers():
        url = p.get("url", "")
        name = p.get("name", url)
        if not url:
            continue
        # Skip self — already added as the local node above
        if self_url and url == self_url:
            continue
        if url in responding_urls:
            raw = peers_metrics[url]
            nodes[name] = {
                "url": url,
                "status": "online",
                "queues": {
                    q: {
                        **v,
                        "available": max(
                            0, v["configured"] - v["running"] - v["queued"]
                        ),
                    }
                    for q, v in raw.items()
                },
            }
        else:
            nodes[name] = {"url": url, "status": "offline"}

    # Totals across all online nodes
    totals: dict = {}
    for node_data in nodes.values():
        if node_data.get("status") != "online":
            continue
        for q, v in node_data.get("queues", {}).items():
            if q not in totals:
                totals[q] = {"configured": 0, "running": 0, "queued": 0, "available": 0}
            totals[q]["configured"] += v.get("configured", 0)
            totals[q]["running"] += v.get("running", 0)
            totals[q]["queued"] += v.get("queued", 0)
            totals[q]["available"] += v.get("available", 0)

    return {"nodes": nodes, "totals": totals}


@router.get("/cluster/tasks")
async def cluster_tasks():
    """Return the consolidated task list from all nodes (self + peers).

    Each task carries two extra fields compared to GET /list:
        "node":     "node1"                      human-readable name from peers.json
        "node_url": "http://192.168.1.10:4200"   URL of the originating node

    Tasks are sorted by state priority (running → queued → waiting → finished)
    then by task id within each state.
    """
    from routers.queue import list_task  # local import to avoid circular deps

    peer_name_map = {
        p["url"]: p.get("name", p["url"]) for p in load_peers() if p.get("url")
    }
    self_url = get_self_url()
    self_name = (
        peer_name_map.get(self_url, get_self_hostname())
        if self_url
        else get_self_hostname()
    )

    # Local tasks
    local_response = await list_task()
    local_tasks = _json.loads(local_response.body)
    for t in local_tasks:
        t["node"] = self_name
        t["node_url"] = self_url

    # Peers tasks
    peers_tasks = await get_peers_tasks()  # {url: {"name": str, "tasks": list}}
    all_tasks = list(local_tasks)
    for url, data in peers_tasks.items():
        for t in data.get("tasks", []):
            t["node"] = data.get("name", url)
            t["node_url"] = url
            all_tasks.append(t)

    # Sort: running first, then queued, then others, then finished
    _state_order = {"running": 0, "queued": 1, "waiting": 2, "finished": 3}
    all_tasks.sort(
        key=lambda t: (
            _state_order.get(t.get("state", "").lower(), 4),
            t.get("id") or 0,
        )
    )

    return all_tasks
