#!/usr/bin/env python

import json as _json
from typing import Optional, Union

import httpx  # pyright: ignore[reportMissingImports]
from fastapi import (
    APIRouter,
    Depends,
    HTTPException,
    Query,
    Response,
    status,
)  # pyright: ignore[reportMissingImports]

from authentication import get_current_user_or_service
from config import STARK_API_KEY, local_port
from models import User
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

    if not self_url:
        self_url = "http://" + get_self_hostname() + f":{local_port}"

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


@router.get("/queues")
def list_queues():
    """Return the list of configured queue names."""
    from queues import load_queues

    return {"queues": list(load_queues().keys())}


@router.get("/modules")
def list_modules():
    """Return the list of configured modules with their public metadata.

    The Docker image name is intentionally excluded from the response to
    avoid leaking internal infrastructure details to clients.
    """
    from modules import load_modules

    return {
        "modules": [
            {
                "name": name,
                "description": cfg.get("description", ""),
                "defaults": cfg.get("defaults", {}),
            }
            for name, cfg in load_modules().items()
        ]
    }


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

    if not self_url:
        self_url = "http://" + get_self_hostname() + f":{local_port}"

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

    # # Sort: running first, then queued, then others, then finished
    # _state_order = {"running": 0, "queued": 1, "waiting": 2, "finished": 3}
    # all_tasks.sort(
    #     key=lambda t: (
    #         _state_order.get(t.get("state", "").lower(), 4),
    #         t.get("id") or 0,
    #     )
    # )

    return all_tasks


@router.get("/cluster/proxy/queue")
async def cluster_proxy_queue(
    node_url: str = Query(..., description="Target node base URL"),
    action: str = Query(...),
    id: Optional[str] = Query(None),
    queue: Optional[str] = Query(None),
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """Proxy a /queue action to a specific cluster node.

    Used by the cluster tasks UI to forward Info / Log / Analysis / Kill /
    Prioritize / Remove actions to the node that owns the task.
    Destructive actions require admin or service key on this node.
    """
    restricted = {"kill", "prioritize", "remove", "swap"}
    if action in restricted and authorized != "service":
        if not isinstance(authorized, User) or "admin" not in authorized.groups:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin group required for this action",
            )

    # get URL
    if not node_url or node_url is None or node_url == "null":
        node_url = "http://" + get_self_hostname() + f":{local_port}"

    params = f"action={action}"
    if id:
        params += f"&id={id}"
    if queue:
        params += f"&queue={queue}"
    url = f"{node_url}/queue?{params}"

    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            resp = await client.get(url, headers={"X-API-Key": STARK_API_KEY})
        return Response(
            content=resp.content,
            status_code=resp.status_code,
            media_type=resp.headers.get("content-type", "text/plain"),
        )
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"Node unreachable: {exc}") from exc


@router.post("/cluster/proxy/relaunch/{ts_id}")
async def cluster_proxy_relaunch(
    ts_id: str,
    node_url: str = Query(..., description="Target node base URL"),
    queue: Optional[str] = Query(None),
    prioritize: Optional[bool] = Query(None),
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """Proxy a /relaunch request to a specific cluster node."""
    if authorized != "service":
        if not isinstance(authorized, User) or "admin" not in authorized.groups:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin group required for relaunch",
            )

    # get URL
    if not node_url or node_url is None or node_url == "null":
        node_url = "http://" + get_self_hostname() + f":{local_port}"

    q = f"?queue={queue}" if queue else ""
    if prioritize is not None:
        q += f"&prioritize={prioritize}" if q else f"?prioritize={prioritize}"

    url = f"{node_url}/relaunch/{ts_id}{q}"

    try:
        async with httpx.AsyncClient(timeout=15.0) as client:
            resp = await client.post(url, headers={"X-API-Key": STARK_API_KEY})
        return Response(
            content=resp.content,
            status_code=resp.status_code,
            media_type=resp.headers.get("content-type", "text/plain"),
        )
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"Node unreachable: {exc}") from exc


@router.get("/cluster/proxy/archive/{analysis_id_name}/log")
async def cluster_proxy_archive_log(
    analysis_id_name: str,
    node_url: str = Query(..., description="Target node base URL"),
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """Proxy GET /archive/{id}/log to the node that owns the archive files."""
    if not node_url or node_url == "null":
        node_url = "http://" + get_self_hostname() + f":{local_port}"
    url = f"{node_url}/archive/{analysis_id_name}/log"
    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            resp = await client.get(url, headers={"X-API-Key": STARK_API_KEY})
        return Response(
            content=resp.content,
            status_code=resp.status_code,
            media_type=resp.headers.get("content-type", "text/plain"),
        )
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"Node unreachable: {exc}") from exc


@router.get("/cluster/proxy/archive/{analysis_id_name}/json")
async def cluster_proxy_archive_json(
    analysis_id_name: str,
    node_url: str = Query(..., description="Target node base URL"),
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """Proxy GET /archive/{id}/json to the node that owns the archive files."""
    if not node_url or node_url == "null":
        node_url = "http://" + get_self_hostname() + f":{local_port}"
    url = f"{node_url}/archive/{analysis_id_name}/json"
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            resp = await client.get(url, headers={"X-API-Key": STARK_API_KEY})
        return Response(
            content=resp.content,
            status_code=resp.status_code,
            media_type=resp.headers.get("content-type", "text/plain"),
        )
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"Node unreachable: {exc}") from exc


@router.post("/cluster/proxy/archive/{analysis_id_name}/relaunch")
async def cluster_proxy_archive_relaunch(
    analysis_id_name: str,
    node_url: str = Query(..., description="Target node base URL"),
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """Proxy POST /archive/{id}/relaunch to the node that owns the archive files."""
    if authorized != "service":
        if not isinstance(authorized, User) or "admin" not in authorized.groups:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin group required for relaunch",
            )
    if not node_url or node_url == "null":
        node_url = "http://" + get_self_hostname() + f":{local_port}"
    url = f"{node_url}/archive/{analysis_id_name}/relaunch"
    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            resp = await client.post(url, headers={"X-API-Key": STARK_API_KEY})
        return Response(
            content=resp.content,
            status_code=resp.status_code,
            media_type=resp.headers.get("content-type", "text/plain"),
        )
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"Node unreachable: {exc}") from exc


@router.delete("/cluster/proxy/archive/{analysis_id_name}")
async def cluster_proxy_archive_delete(
    analysis_id_name: str,
    node_url: str = Query(..., description="Target node base URL"),
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """Proxy DELETE /archive/{id} to the node that owns the archive files."""
    if authorized != "service":
        if not isinstance(authorized, User) or "admin" not in authorized.groups:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin group required for delete",
            )
    if not node_url or node_url == "null":
        node_url = "http://" + get_self_hostname() + f":{local_port}"
    url = f"{node_url}/archive/{analysis_id_name}"
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            resp = await client.delete(url, headers={"X-API-Key": STARK_API_KEY})
        return Response(
            content=resp.content,
            status_code=resp.status_code,
            media_type=resp.headers.get("content-type", "text/plain"),
        )
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"Node unreachable: {exc}") from exc


@router.get("/cluster/archives")
async def cluster_archives():
    """Return the aggregated archives from all nodes (self + peers).

    Each entry is enriched with 'node' and 'node_url' fields.
    Results are sorted by mtime descending (newest first).
    """
    import asyncio

    from routers.queue import list_archives  # local import to avoid circular deps

    peer_name_map = {
        p["url"]: p.get("name", p["url"]) for p in load_peers() if p.get("url")
    }
    self_url = get_self_url()
    self_name = (
        peer_name_map.get(self_url, get_self_hostname())
        if self_url
        else get_self_hostname()
    )

    if not self_url:
        self_url = "http://" + get_self_hostname() + f":{local_port}"

    # Local archives
    local_response = await list_archives()
    local_archives = _json.loads(local_response.body)
    for a in local_archives:
        a["node"] = self_name
        a["node_url"] = self_url

    # Peer archives
    peers = load_peers()
    peer_entries = [
        (p.get("name", p.get("url", "")), p["url"])
        for p in peers
        if p.get("url") and (not self_url or p["url"] != self_url)
    ]

    all_archives = list(local_archives)

    if peer_entries:
        async with httpx.AsyncClient(timeout=10.0) as client:
            responses = await asyncio.gather(
                *[
                    client.get(f"{url}/archives", headers={"X-API-Key": STARK_API_KEY})
                    for _, url in peer_entries
                ],
                return_exceptions=True,
            )
        for (name, url), resp in zip(peer_entries, responses):
            if isinstance(resp, Exception):
                continue
            if resp.status_code == 200:
                try:
                    archives = resp.json()
                    for a in archives:
                        a["node"] = name
                        a["node_url"] = url
                    all_archives.extend(archives)
                except Exception:
                    pass

    all_archives.sort(key=lambda a: a.get("mtime", 0), reverse=True)
    return all_archives
