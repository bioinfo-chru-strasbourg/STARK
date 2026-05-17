#!/usr/bin/env python
"""
nodes.py — Cluster / mini-orchestrator logic.

Responsibilities:
  - Load nodes from config/nodes.json (auto-created empty on first start).
  - Resolve this node's own URL (env var → /whoami scan → None = local-only).
  - Compute local task-spooler metrics per queue.
  - Collect metrics from all nodes asynchronously (httpx).
  - Score nodes and select the best one for a given queue.
  - Forward an /analysis request transparently to a target node.

If no nodes are configured, or if self-URL cannot be determined, all analyses
are run locally — the cluster logic is entirely transparent / opt-in.

nodes.json format (auto-created in config/ on first start):
{
  "nodes": [
    {"name": "node1", "url": "http://192.168.1.10:8001"},
    {"name": "node2", "url": "http://192.168.1.10:8002"},
    {"name": "node3", "url": "http://192.168.1.11:8001"}
  ]
}
"""

import asyncio
import json
import logging
import os
import socket
import subprocess
from typing import Optional

import httpx  # pyright: ignore[reportMissingImports]

logger = logging.getLogger(__name__)

from config import NODES_FILE, STARK_API_KEY, STARK_API_SELF_URL, shell, ts
from queues import _queue_daemon_is_active, _resolve_queue_socket, load_queues
from tasks import _read_task_slots, _resolve_task_slots

# from routers.cluster import whoami  # pyright: ignore[reportMissingImports]

# ---------------------------------------------------------------------------
# Self-identity
# ---------------------------------------------------------------------------

_self_hostname: str = socket.gethostname()
_self_url_cache: Optional[str] = None  # populated lazily


def get_self_hostname() -> str:
    return _self_hostname


def get_self_url() -> Optional[str]:
    """Return the URL at which *this* node is reachable.

    Resolution order:
      1. STARK_API_SELF_URL environment variable (explicit, recommended).
      2. Lazy scan of configured nodes: POST /whoami to each node; if the
         returned hostname matches ours, that node entry is our own URL.
      3. None — self-URL unknown; analyses will always run locally.

    The result is cached after the first successful resolution.
    """
    global _self_url_cache

    # 1. Explicit env var
    if STARK_API_SELF_URL:
        return STARK_API_SELF_URL

    # Return cached value if already resolved
    if _self_url_cache is not None:
        return _self_url_cache

    # 2. Scan nodes to find which one is us
    nodes = load_nodes()
    for p in nodes:
        url = p.get("url", "")
        if not url:
            continue
        if not url.startswith("http"):
            url = f"http://{url}"
        try:
            r = httpx.get(f"{url}/whoami", timeout=1.0)
            if r.is_success and r.json().get("id") == _self_hostname:
                _self_url_cache = url
                return url
        except Exception:
            continue

    # 3. Unknown — local-only mode
    # return "http://localhost"
    return None


# ---------------------------------------------------------------------------
# Nodes config
# ---------------------------------------------------------------------------

def load_nodes() -> list:
    """Load nodes from config/nodes.json.

    If the file does not exist it is auto-created with an empty nodes list so
    the user has a ready-to-fill template. Returns [] if file is absent,
    malformed, or contains no nodes — in all those cases analyses run locally.

    Expected format:
    {
      "nodes": [
        {"name": "node1", "url": "http://192.168.1.10:8001"},
        {"name": "node2", "url": "http://192.168.1.11:8001"}
      ]
    }
    """
    default: dict = {"nodes": []}

    if not os.path.exists(NODES_FILE):
        try:
            os.makedirs(os.path.dirname(NODES_FILE), exist_ok=True)
            with open(NODES_FILE, "w") as f:
                json.dump(default, f, indent=2)
        except OSError:
            pass
        return []

    try:
        with open(NODES_FILE, "r") as f:
            data = json.load(f)
        return data.get("nodes", [])
    except (FileNotFoundError, json.JSONDecodeError):
        return []


# ---------------------------------------------------------------------------
# Local metrics
# ---------------------------------------------------------------------------

def get_local_metrics() -> dict:
    """Return task-spooler slot metrics for every configured queue.

    Format:
        {
            "stark":  {"configured": 1, "running": 0, "queued": 0},
            "light":  {"configured": 4, "running": 2, "queued": 1},
        }

    Counts are derived from 'ts -l' output.  If a queue daemon is inactive,
    running/queued default to 0 (the queue is considered free but offline).
    """
    queues = load_queues()
    default_name = next(iter(queues))
    metrics: dict = {}

    for queue_name, cfg in queues.items():
        is_default = queue_name == default_name
        configured = int(cfg.get("slots", 1))

        if not _queue_daemon_is_active(queue_name, cfg, is_default):
            metrics[queue_name] = {"configured": configured, "running": 0, "queued": 0}
            continue

        socket_path = _resolve_queue_socket(queue_name, cfg, is_default)
        socket_part = f"TS_SOCKET={socket_path} " if socket_path else ""
        q_env = f"{socket_part}TS_SAVELIST={cfg['savelist']} TS_SLOTS={cfg['slots']} "

        running = 0
        queued = 0
        try:
            result = subprocess.run(
                f"{q_env} {ts} -l",
                shell=True,
                executable=shell,
                capture_output=True,
                text=True,
                check=False,
                timeout=5,
            )
            for line in result.stdout.strip().split("\n")[1:]:
                parts = line.split()
                if len(parts) < 2:
                    continue
                state = parts[1].lower()
                # Count slot consumption, not task count:
                # a task submitted with -N k occupies k slots.
                slot_cost = _read_task_slots(line, configured)
                if slot_cost is None or slot_cost == 0:
                    slot_cost = 1  # fallback: unknown → assume 1 slot
                if state == "running":
                    running += slot_cost
                elif state == "queued":
                    queued += slot_cost
        except Exception:
            pass

        metrics[queue_name] = {
            "configured": configured,
            "running": running,
            "queued": queued,
        }

    return metrics


# ---------------------------------------------------------------------------
# Node metrics (async)
# ---------------------------------------------------------------------------

async def get_nodes_metrics() -> dict:
    """Collect /metrics from all configured nodes concurrently.

    Returns {node_url: queues_dict}, where queues_dict has the same structure
    as get_local_metrics().  Unreachable nodes are silently skipped.
    """
    nodes = load_nodes()
    if not nodes:
        return {}

    headers = {"X-API-Key": STARK_API_KEY}
    # Exclude self to avoid double-counting when this node is in nodes.json
    self_url = get_self_url()
    urls = [
        p["url"]
        for p in nodes
        if p.get("url") and (not self_url or p["url"] != self_url)
    ]

    logger.debug("Collecting metrics from nodes: %s", urls)

    async with httpx.AsyncClient(timeout=1.0) as client:
        responses = await asyncio.gather(
            *[client.get(f"{url}/metrics", headers=headers) for url in urls],
            return_exceptions=True,
        )

    result: dict = {}
    for url, resp in zip(urls, responses):
        if isinstance(resp, Exception):
            logger.debug("Node %s unreachable: %s", url, resp)
            continue
        if resp.status_code == 200:
            try:
                data = resp.json()
                result[url] = data.get("queues", {})
            except Exception:
                continue

    return result


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


def compute_best_node(
    queue_name: str, all_metrics: dict, json_input: dict
) -> Optional[str]:
    """Return the URL of the node with the most capacity for queue_name.

    Score = configured - running - queued  (higher → more slots available).
    Returns None if all_metrics is empty.
    """

    # If no metrics, return None to run locally
    if not all_metrics:
        return None

    # Extract requested slots from input, default to None (no specific requirement)
    requested_slots = json_input.get("threads", None)

    # Determine max_slots across cluster for this queue to use as a fallback when requested_slots is not specified or exceeds max capacity.
    max_slots = max(
        q.get(queue_name, {}).get("configured", 0) for q in all_metrics.values()
    )

    # fallback = max capacity across cluster
    if requested_slots is None or requested_slots > max_slots:
        requested_slots = max_slots

    # Track best node and its priority tuple for comparison
    best_url: Optional[str] = None
    best_priority: Optional[tuple] = None

    for url, queues in all_metrics.items():

        # Get metrics for the requested queue
        q = queues.get(queue_name, {})

        # If queue not configured on that node, consider it unavailable for that node (score = -inf)
        if not q or q.get("configured", 0) <= 0:
            continue

        # Extract metrics for easier readability
        configured = q.get("configured", 0)
        running = q.get("running", 0)
        queued = q.get("queued", 0)

        # Calculate score based on configured slots minus running and queued slots
        available = configured - running - queued

        # overload distance (can be negative)
        delta = available - requested_slots

        # PRIORITY DESIGN
        # 1. First classify nodes by capability: can this node satisfy the request (delta >= 0) or not (delta < 0)?
        # 2. Among capable nodes, prefer those that are not overloaded (delta < 0) over those that are (delta >= 0).
        # 3. Then minimize overload (negative delta) or waste (positive delta).
        # 4. Penalize overload more heavily than waste by making it a primary sorting key.

        is_capable = configured >= requested_slots
        is_overloaded = delta < 0
        overload_penalty = abs(delta) if is_overloaded else 0
        slack_penalty = delta if delta > 0 else 0

        # final score:
        priority = (
            not is_capable,  # True (good) > False (bad)
            is_overloaded,  # False (good) < True (bad)
            overload_penalty,  # minimize overload
            slack_penalty,  # minimize wasted slots
        )

        if best_priority is None or priority < best_priority:
            best_priority = priority
            best_url = url

    if best_url is None:
        # Choose first node with the queue configured, even if it has no free slots (best effort)
        for url, queues in all_metrics.items():
            q = queues.get(queue_name, {})
            if q and q.get("configured", 0) > 0:
                best_url = url
                break

    return best_url


# ---------------------------------------------------------------------------
# Forward
# ---------------------------------------------------------------------------

async def forward_request(target_url: str, request):
    """Transparently forward a request to target_url, preserving the original path.

    Adds X-STARK-Forwarded: 1 to prevent routing loops (a forwarded request is
    always executed locally on the receiving node).
    """
    from fastapi.responses import Response  # pyright: ignore[reportMissingImports]

    headers = {
        k: v
        for k, v in request.headers.items()
        if k.lower() not in ("host", "content-length")
    }
    headers["X-STARK-Forwarded"] = "1"

    async with httpx.AsyncClient(timeout=30.0) as client:
        resp = await client.post(
            f"{target_url}{request.url.path}" + ('?' + request.url.query if request.url.query else ''),
            content=await request.body(),
            headers=headers,
        )

    return Response(
        content=resp.content,
        status_code=resp.status_code,
        media_type=resp.headers.get("content-type", "text/plain"),
    )


# ---------------------------------------------------------------------------
# Node tasks (async)
# ---------------------------------------------------------------------------


async def get_nodes_tasks() -> dict:
    """Collect /list (task list) from all configured nodes concurrently.

    Returns {node_url: tasks_list}, where each tasks_list is the same
    structure returned by the local GET /list endpoint.
    Unreachable nodes are silently skipped.

    The name of each node (from nodes.json) is included so callers can
    map a URL back to a human-readable name.
    """
    nodes = load_nodes()
    if not nodes:
        return {}

    headers = {"X-API-Key": STARK_API_KEY}
    # Exclude self to avoid duplicate tasks when this node is in nodes.json
    self_url = get_self_url()
    node_entries = [
        (p.get("name", p.get("url", "")), p["url"])
        for p in nodes
        if p.get("url") and (not self_url or p["url"] != self_url)
    ]

    async with httpx.AsyncClient(timeout=5.0) as client:
        responses = await asyncio.gather(
            *[client.get(f"{url}/list", headers=headers) for _, url in node_entries],
            return_exceptions=True,
        )

    result: dict = {}
    for (name, url), resp in zip(node_entries, responses):
        if isinstance(resp, Exception):
            logger.debug("Node %s (%s) unreachable for tasks: %s", name, url, resp)
            continue
        if resp.status_code == 200:
            try:
                result[url] = {"name": name, "tasks": resp.json()}
            except Exception:
                continue

    return result
