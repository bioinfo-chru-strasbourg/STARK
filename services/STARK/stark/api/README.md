# STARK API

A **FastAPI**-based web service for launching and monitoring [STARK](https://github.com/bioinfo-chru-strasbourg/STARK) bioinformatics analyses via a task queue ([task-spooler](https://github.com/thomaspreece/task-spooler)) and Docker.

![STARK API dashboard](images/api.png)

---

## Features

- **Web UI** — browser-based dashboard to launch analyses, monitor the task queue, inspect logs and results
- **REST API** — JSON endpoints consumable by external services (e.g. `STARK.listener`)
- **Dual authentication** — JWT bearer tokens for human users + static API key (`X-API-Key`) for service-to-service calls
- **Role-based access** — `admin` group required for destructive actions (kill, remove, prioritize, relaunch) and for launching new analyses
- **Live queue** — auto-refreshing task table with state-aware action buttons
- **Inline detail panels** — Info / Log / Analysis JSON displayed inline, persisted across auto-refreshes, with one-click copy
- **Time tracking** — elapsed time for running tasks (computed from `Start time`), final duration for finished tasks (from `Time run`)
- **Relaunch** — re-queue a finished task from its original JSON parameters

---

## Stack

| Component | Role |
|---|---|
| FastAPI + Uvicorn | HTTP server |
| Jinja2 | HTML templating (index.html) |
| python-jose | JWT encoding/decoding |
| passlib[bcrypt] | Password hashing |
| task-spooler (`ts`) | Job queue |
| Docker | STARK analysis runtime |

---

## Directory layout

```
dockerfile/
├── app.py               # FastAPI application
├── requirements.txt
├── Dockerfile
├── config/
│   └── users.json       # User database (auto-created on first run)
├── static/
│   ├── script.js        # Frontend logic
│   └── style.css
└── templates/
    └── index.html       # Dashboard template
images/
    └── api.png
STARK.env                # Environment variable defaults
README.md
```

---

## Configuration

All settings are passed as **environment variables** (e.g. via `docker-compose` or a `.env` file).

| Variable | Default | Description |
|---|---|---|
| `STARK_API_KEY` | `a_default_super_secret_api_key` | Static API key for service-to-service auth (`X-API-Key` header) |
| `STARK_API_REFRESH_INTERVAL` | `10` | Queue auto-refresh interval **in seconds** |
| `DOCKER_STARK_IMAGE` | `stark` | Docker image used to run analyses |
| `TS` | _(empty)_ | Path to the `ts` binary |
| `TS_SAVELIST` | `/ts-tmp` | task-spooler save directory |
| `TS_SLOTS` | `1` | Number of parallel task-spooler slots |
| `SHELL` | `/bin/ash` | Shell used to run queued commands |
| `DOCKER_STARK_SERVICE_STARK_API_CONTAINER_MOUNT` | _(empty)_ | Extra Docker volume mounts for the STARK container |
| `DOCKER_STARK_SERVICE_STARK_API_LOG_FOLDER` | `/STARK/services/stark/stark/api` | Folder where `.json`, `.info` and `.output` files are written |
| `DOCKER_STARK_SERVICE_STARK_API_RUNS_FOLDER` | `/STARK/input/runs` | Default folder scanned for run directories |

Example `STARK.env`:
```env
STARK_API_KEY=my_secret_key
STARK_API_REFRESH_INTERVAL=60
```

---

## User management

Users are stored in `config/users.json` (auto-created with default accounts on first start).

```json
{
  "stark": {
    "password": "password",
    "groups": ["admin"]
  },
  "guest": {
    "password": "guest",
    "groups": []
  }
}
```

- Users in the `admin` group can **launch**, **relaunch**, **kill**, **prioritize**, and **remove** tasks.
- Users without `admin` can only **view** the queue, and display **Info / Log / Analysis** details.

---

## REST API

### Authentication

**Obtain a JWT token:**
```bash
TOKEN=$(curl -s -X POST "http://localhost:8000/token" \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=stark&password=password" \
  | python3 -c "import sys,json; print(json.load(sys.stdin)['access_token'])")
```

All subsequent requests use:
```
Authorization: Bearer $TOKEN
```

Service-to-service calls use instead:
```
X-API-Key: <STARK_API_KEY>
```

---

### Endpoints

#### `GET /`
Returns the web dashboard (HTML).

---

#### `POST /token`
Obtain a JWT access token.

| Field | Type | Description |
|---|---|---|
| `username` | form | Username |
| `password` | form | Password |

**Response:** `{ "access_token": "...", "token_type": "bearer" }`

---

#### `GET /me`
Returns the current user's username and groups.

**Response:** `{ "username": "stark", "groups": ["admin"] }`

---

#### `POST /analysis`
Launch a new STARK analysis. **Admin only.**

**Body (JSON):**
```json
{ "run": "MY_RUN_NAME" }
```

**Response:** `STARK.<ID>.<analysisIDNAME>` (plain text, 200) or error (400/403).

**Example:**
```bash
curl -s -X POST "http://localhost:8000/analysis" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"run": "MY_RUN"}'
```

---

#### `GET /queue`
Query or act on the task queue.

| `action` | Admin required | Description |
|---|---|---|
| `list` | No | Returns the full queue as a JSON array |
| `info` | No | `ts -i <id>`: raw task info |
| `log` | No | Contents of the `.output` log file |
| `analysis` | No | Contents of the original `.json` parameters file (pretty-printed) |
| `kill` | Yes | Stop the Docker container + `ts -k <id>` |
| `prioritize` | Yes | `ts -u <id>`: move task to front of queue |
| `remove` | Yes | `ts -r <id>`: remove task from queue |
| `swap` | Yes | `ts -U <id>`: swap with next task |

**List tasks:**
```bash
curl -s "http://localhost:8000/queue?action=list" \
  -H "Authorization: Bearer $TOKEN" | python3 -m json.tool
```

**Get task log:**
```bash
curl -s "http://localhost:8000/queue?action=log&id=3" \
  -H "Authorization: Bearer $TOKEN"
```

**Kill a task:**
```bash
curl -s "http://localhost:8000/queue?action=kill&id=3" \
  -H "Authorization: Bearer $TOKEN"
```

**List response format:**
```json
[
  {
    "id": "2",
    "state": "running",
    "output": "/ts-tmp/...",
    "elevel": "",
    "times": "142.3",
    "run_name": "MY_RUN"
  }
]
```

---

#### `POST /relaunch/{ts_id}`
Re-queue a finished task using its original JSON parameters. **Admin only.**

```bash
curl -s -X POST "http://localhost:8000/relaunch/3" \
  -H "Authorization: Bearer $TOKEN"
```

**Response:** `STARK.<new_ID>.<analysisIDNAME>` (plain text, 200).

---

## Web UI

The dashboard is accessible at `http://localhost:8000/`.

| Section | Description |
|---|---|
| **Launch Analysis** | Run name field + Advanced JSON textarea. Visible to admins only. |
| **Task Queue** | Auto-refreshing table. Columns: ID, State, E-Level, Time, Run Name, Actions. |

**Action buttons per state:**

| State | Buttons |
|---|---|
| `running` | Info, Log, Analysis, **Kill** |
| `queued` | Info, Log, Analysis, **Prioritize**, **Remove** |
| `finished` | Info, Log, Analysis, **Relaunch** |

Red buttons (Kill, Prioritize, Remove, Relaunch) are only visible to `admin` users.

Clicking **Info**, **Log**, or **Analysis** opens an inline detail panel below the row. The panel persists across auto-refreshes and includes a **Copy** button.
