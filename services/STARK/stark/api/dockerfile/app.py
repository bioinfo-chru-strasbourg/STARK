#!/usr/bin/env python

from fastapi import (
    FastAPI,
    Request,
    Depends,
    HTTPException,
    status,
    Query,
    Form,
    Security,
)
from fastapi.responses import HTMLResponse, PlainTextResponse, JSONResponse
from fastapi.templating import Jinja2Templates
from fastapi.security import (
    OAuth2PasswordBearer,
    OAuth2PasswordRequestForm,
    APIKeyHeader,
)
from fastapi.staticfiles import StaticFiles
from jose import JWTError, jwt
from passlib.context import CryptContext
from pydantic import BaseModel
import os
import subprocess
import json
import random
import string
import time
import datetime
import re
import concurrent.futures
from typing import Optional, Union, List

# --- Configuration ---
SECRET_KEY = "a_very_secret_key_that_should_be_in_a_config_file"
ALGORITHM = "HS256"
# USERS_FILE = "users.json"
USERS_FILE = "config/users.json"
STARK_API_KEY = os.environ.get("STARK_API_KEY", "a_default_super_secret_api_key")

pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")

app = FastAPI()

app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates")


# --- Models ---
class Token(BaseModel):
    access_token: str
    token_type: str


class User(BaseModel):
    username: str
    groups: List[str] = []


# --- Authentication ---

api_key_header = APIKeyHeader(name="X-API-Key", auto_error=False)


def get_api_key(api_key: str = Security(api_key_header)):
    if api_key == STARK_API_KEY:
        return api_key
    else:
        return None


def load_users():
    """Loads users from the JSON file.
    
    Expected format:
    {
        "stark": {"password": "secret", "groups": ["admin", "users"]},
        "john":  {"password": "pass",   "groups": ["users"]}
    }
    Legacy format (plain password string) is also supported.
    """
    #default = {"stark": {"password": "password", "groups": ["admin"]}}
    default = {
        "stark": {"password": "password", "groups": ["admin"]},
        "guest": {"password": "guest", "groups": []},
    }
    if not os.path.exists(USERS_FILE):
        os.makedirs(os.path.dirname(USERS_FILE), exist_ok=True)
        with open(USERS_FILE, "w") as f:
            json.dump(default, f, indent=2)
    try:
        with open(USERS_FILE, "r") as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return default


def authenticate_user(username: str, password: str):
    """Authenticates a user. Supports both legacy (plain string) and new (dict) formats."""
    users_db = load_users()
    if username not in users_db:
        return False
    entry = users_db[username]
    # New format: {"password": "...", "groups": [...]}
    if isinstance(entry, dict):
        if entry.get("password") == password:
            return User(username=username, groups=entry.get("groups", []))
    # Legacy format: plain password string
    elif entry == password:
        return User(username=username, groups=[])
    return False


def create_access_token(data: dict):
    to_encode = data.copy()
    encoded_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)
    return encoded_jwt


oauth2_scheme = OAuth2PasswordBearer(tokenUrl="token", auto_error=False)


async def get_current_user(token: str = Depends(oauth2_scheme)):
    if not token:
        return None
    credentials_exception = HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Could not validate credentials",
        headers={"WWW-Authenticate": "Bearer"},
    )
    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
        username: str = payload.get("sub")
        if username is None:
            raise credentials_exception
    except JWTError:
        raise credentials_exception

    users_db = load_users()
    if username not in users_db:
        raise credentials_exception
    entry = users_db[username]
    groups = entry.get("groups", []) if isinstance(entry, dict) else []
    return User(username=username, groups=groups)


async def get_current_user_or_service(
    api_key: str = Security(api_key_header), user: User = Depends(get_current_user)
):
    if api_key and api_key == STARK_API_KEY:
        return "service"
    if user:
        return user
    raise HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED, detail="Not authenticated"
    )


# --- Routes ---


@app.post("/token", response_model=Token)
async def login_for_access_token(form_data: OAuth2PasswordRequestForm = Depends()):
    user = authenticate_user(form_data.username, form_data.password)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect username or password",
            headers={"WWW-Authenticate": "Bearer"},
        )
    access_token = create_access_token(data={"sub": user.username, "groups": user.groups})
    return {"access_token": access_token, "token_type": "bearer"}


@app.get("/me")
async def get_me(user: User = Depends(get_current_user)):
    """Returns the current user's info (username and groups)."""
    if not user:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Not authenticated")
    return {"username": user.username, "groups": user.groups}


@app.get("/", response_class=HTMLResponse)
async def read_root(request: Request):
    import os

    version = str(int(os.path.getmtime("static/script.js")))
    return templates.TemplateResponse(
        "index.html",
        {
            "request": request,
            "version": version,
            "refresh_interval_ms": refresh_interval_ms,
        },
    )


# --- STARK API LOGIC ---

# Docker STARK image
docker_stark = os.environ.get("DOCKER_STARK_IMAGE", "stark")
ts = os.environ.get("TS", "")
ts_savelist = os.environ.get("TS_SAVELIST", "/ts-tmp")
ts_slots = os.environ.get("TS_SLOTS", "1")
shell = os.environ.get("SHELL", "/bin/ash")
ts_env = f"TS_SAVELIST={ts_savelist} TS_SLOTS={ts_slots} "
docker_stark_container_mount = os.environ.get(
    "DOCKER_STARK_SERVICE_STARK_API_CONTAINER_MOUNT", ""
)
docker_stark_api_log_folder = os.environ.get(
    "DOCKER_STARK_SERVICE_STARK_API_LOG_FOLDER", "/STARK/services/stark/stark/api"
)
docker_stark_api_runs_folder = os.environ.get(
    "DOCKER_STARK_SERVICE_STARK_API_RUNS_FOLDER", "/STARK/input/runs"
)
refresh_interval_ms = int(os.environ.get("STARK_API_REFRESH_INTERVAL", "10")) * 1000


def randomStringDigits(stringLength=6):
    """Generate a random string of letters and digits"""
    lettersAndDigits = string.ascii_letters + string.digits
    return "".join(random.choice(lettersAndDigits) for i in range(stringLength))


def queue_analysis(json_input: dict) -> str:
    """Build and queue a STARK Docker analysis from a JSON input dict.
    Returns the analysis IDNAME string on success, or raises RuntimeError."""
    json_dump = json.dumps(json_input)

    analysesID = randomStringDigits(12)
    runID = "UNKNOWN"
    runMD5 = randomStringDigits(41)

    if "run" in json_input:
        name = json_input["run"].split(":")[0]
        runID = os.path.basename(name)

        runFolder = ""
        if os.path.isdir(name):
            runFolder = name
        elif os.path.isdir(os.path.join(docker_stark_api_runs_folder, name)):
            runFolder = os.path.join(docker_stark_api_runs_folder, name)

        if runFolder:
            myCmd = f"find {runFolder} -maxdepth 1 -type f -print0 | xargs -0 sha1sum | cut -b-40 | sha1sum | awk '{{print $1}}'"
        else:
            myCmd = f"echo {name} | sha1sum | awk '{{print $1}}'"

        runMD5 = (
            subprocess.run(myCmd, shell=True, stdout=subprocess.PIPE)
            .stdout.decode("utf-8")
            .strip()
        )

        if "analysis_name" in json_input:
            runID = json_input["analysis_name"]

    elif "analysis_name" in json_input:
        runID = json_input["analysis_name"]
        runMD5 = randomStringDigits(41)

    analysesRUNNAME = runID
    analysesNAME = f"ID-{runMD5}-NAME-{runID}"
    analysisIDNAME = f"STARK.{analysesID}.{analysesNAME}"

    docker_name = f" --name {analysisIDNAME} "
    docker_parameters = f" --rm {docker_stark_container_mount} {docker_name} "

    analysisFOLDER = docker_stark_api_log_folder
    analysisFILE = os.path.join(analysisFOLDER, f"{analysisIDNAME}.json")
    analysisINFOFILE = os.path.join(analysisFOLDER, f"{analysisIDNAME}.info")
    analysisOUTPUTFILE = os.path.join(analysisFOLDER, f"{analysisIDNAME}.output")

    with open(analysisFILE, "w") as f:
        f.write(json_dump)

    ts_cmd = f"{ts_env}{ts} -L {analysisIDNAME}" if ts else ""
    myCmd = (
        f'{ts_cmd} sh -c "'
        f"docker run {docker_parameters} {docker_stark} "
        f"--analysis_name={analysesRUNNAME} --analysis={analysisFILE} "
        f"> {analysisOUTPUTFILE} 2>&1 "
        f"&& (echo 'done' > {analysisINFOFILE} && exit 0) "
        f"|| (echo 'failed' > {analysisINFOFILE} && exit 1)\""
    )

    task_id = (
        subprocess.run(myCmd, shell=True, stdout=subprocess.PIPE)
        .stdout.decode("utf-8")
        .strip()
    )
    if not task_id:
        raise RuntimeError("task-spooler returned no task ID")
    return analysisIDNAME


@app.post("/analysis")
async def stark_launch(
    request: Request,
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    # Restrict launch to admin users (service calls are always allowed)
    if authorized != "service":
        if not isinstance(authorized, User) or "admin" not in authorized.groups:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin group required to launch analyses",
            )
    json_input = await request.json()
    try:
        analysisIDNAME = queue_analysis(json_input)
        return PlainTextResponse(content=analysisIDNAME, status_code=200)
    except Exception as e:
        return PlainTextResponse(content=f"KO: {e}", status_code=400)


@app.get("/queue")
async def queue(
    action: str = Query(
        "list",
        enum=[
            "list",
            "info",
            "log",
            "analysis",
            "kill",
            "prioritize",
            "remove",
            "swap",
        ],
    ),
    id: Optional[str] = Query(None),
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """
    Get the task spooler queue or perform an action on a task.
    """
    if action == "list":
        command = f"{ts_env} {ts} -l"  # Use -l to get the list of tasks
        try:
            result = subprocess.run(
                command,
                shell=True,
                executable=shell,
                capture_output=True,
                text=True,
                check=False,
            )

            if result.returncode != 0 and (
                "tasks" in result.stderr.lower() or not result.stdout
            ):
                return JSONResponse(content=[])

            lines = result.stdout.strip().split("\n")
            tasks = []
            if len(lines) > 1:
                # The header line is not always consistent, so we parse based on column content.
                # This regex is designed to be more flexible.
                line_regex = re.compile(
                    r"^(?P<id>\d+)\s+"
                    r"(?P<state>\w+)\s+"
                    r"(?P<output>\S+)\s+"
                    r"(?P<elevel>\S*)\s*"  # Optional E-Level
                    r"(?P<times>[\d./\s-]*)\s+"  # Optional Times
                    r"(?P<command>.*)$"
                )

                for line in lines[1:]:
                    match = line_regex.match(line)
                    if match:
                        data = match.groupdict()
                        run_name_match = re.search(r"-NAME-([^\s\]]+)", data["command"])
                        run_name = run_name_match.group(1) if run_name_match else "N/A"

                        tasks.append(
                            {
                                "id": data["id"].strip(),
                                "state": data["state"].strip(),
                                "output": data["output"].strip(),
                                "elevel": (
                                    "SUCCESS"
                                    if data["elevel"].strip() == "0"
                                    else (
                                        "FAILED: " + data["elevel"].strip()
                                        if data["elevel"].isdigit()
                                        else ""
                                    )
                                ),
                                "times": "",
                                "run_name": run_name,
                            }
                        )

            # Enrich running and finished tasks with Time run from ts -i (parallel)
            def fetch_time(task):
                state = task["state"].lower()
                if state not in ("running", "finished"):
                    return
                try:
                    info_result = subprocess.run(
                        f"{ts_env} {ts} -i {task['id']}",
                        shell=True,
                        executable=shell,
                        capture_output=True,
                        text=True,
                        check=False,
                        timeout=5,
                    )
                    stdout = info_result.stdout
                    if state == "finished":
                        # ts -i reports "Time run: X.XXXs" once the task is done
                        time_match = re.search(r"Time run:\s*([\d.]+)s", stdout)
                        task["times"] = time_match.group(1) if time_match else ""
                    else:
                        # For running tasks, compute elapsed = now - Start time
                        start_match = re.search(r"Start time:\s*(.+)", stdout)
                        if start_match:
                            start_str = start_match.group(1).strip()
                            try:
                                # Parse "Mon Mar 30 18:24:26 2026" style
                                start_dt = datetime.datetime.strptime(
                                    start_str, "%a %b %d %H:%M:%S %Y"
                                )
                                elapsed = time.time() - start_dt.timestamp()
                                task["times"] = f"{elapsed:.1f}"
                            except ValueError:
                                task["times"] = ""
                        else:
                            task["times"] = ""
                except Exception:
                    task["times"] = ""

            tasks_to_enrich = [
                t for t in tasks if t["state"].lower() in ("running", "finished")
            ]
            if tasks_to_enrich:
                with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
                    executor.map(fetch_time, tasks_to_enrich)

            return JSONResponse(content=tasks)
        except FileNotFoundError:
            raise HTTPException(status_code=500, detail=f"Command not found: {ts}")
    else:
        # Restrict destructive actions to admin group
        restricted_actions = {"kill", "prioritize", "remove", "swap"}
        if action in restricted_actions and authorized != "service":
            if not isinstance(authorized, User) or "admin" not in authorized.groups:
                raise HTTPException(
                    status_code=status.HTTP_403_FORBIDDEN,
                    detail="Admin group required for this action",
                )

        # Handling for other actions like info, log, kill
        if not id:
            raise HTTPException(
                status_code=400, detail="Task ID is required for this action"
            )

        action_map = {
            "info": "-i",
            "kill": "-k",
            "prioritize": "-u",
            "remove": "-r",
            "swap": "-U",
        }

        # For kill, also stop the Docker container if it is running
        if action == "kill":
            info_cmd = f"{ts_env} {ts} -i {id}"
            info_result = subprocess.run(
                info_cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=False,
                timeout=10,
            )
            container_match = re.search(r"--name\s+(\S+)", info_result.stdout)
            if container_match:
                container_name = container_match.group(1)
                subprocess.run(
                    f"docker stop {container_name}",
                    shell=True,
                    capture_output=True,
                    text=True,
                    check=False,
                    timeout=30,
                )

        # For analysis: return the original JSON used to launch the task
        elif action == "analysis":
            info_cmd = f"{ts_env} {ts} -i {id}"
            info_result = subprocess.run(
                info_cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=False,
                timeout=10,
            )
            name_match = re.search(r"--name\s+(\S+)", info_result.stdout)
            if not name_match:
                return PlainTextResponse(
                    "Could not determine task name from ts info.", status_code=404
                )
            analysis_name = name_match.group(1)
            json_file = os.path.join(
                docker_stark_api_log_folder, f"{analysis_name}.json"
            )
            if not os.path.isfile(json_file):
                return PlainTextResponse(
                    f"JSON file not found: {json_file}", status_code=404
                )
            try:
                with open(json_file, "r", errors="replace") as f:
                    content = f.read()
                # Pretty-print if valid JSON
                try:
                    content = json.dumps(
                        json.loads(content), indent=2, ensure_ascii=False
                    )
                except json.JSONDecodeError:
                    pass
                return PlainTextResponse(content)
            except OSError as e:
                return PlainTextResponse(
                    f"Error reading JSON file: {e}", status_code=500
                )

        # For log: try docker logs for running containers, fallback to ts output file
        if action == "log":
            # Get task info to extract the container/analysis name
            info_cmd = f"{ts_env} {ts} -i {id}"
            info_result = subprocess.run(
                info_cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=False,
                timeout=10,
            )
            container_match = re.search(r"--name\s+(\S+)", info_result.stdout)
            if not container_match:
                return PlainTextResponse(
                    "Could not determine task name from ts info.", status_code=404
                )

            analysis_name = container_match.group(1)
            output_file = os.path.join(
                docker_stark_api_log_folder, f"{analysis_name}.output"
            )

            if not os.path.isfile(output_file):
                return PlainTextResponse(
                    f"Output file not found: {output_file}", status_code=404
                )
            try:
                with open(output_file, "r", errors="replace") as f:
                    return PlainTextResponse(f.read())
            except OSError as e:
                return PlainTextResponse(
                    f"Error reading output file: {e}", status_code=500
                )

        ts_action = action_map.get(action)
        if not ts_action:
            raise HTTPException(status_code=400, detail=f"Invalid action: {action}")

        command = f"{ts_env} {ts} {ts_action} {id}"
        try:
            result = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                check=False,
                timeout=10,
            )
            response_text = result.stdout if result.stdout else result.stderr
            return PlainTextResponse(response_text)
        except subprocess.TimeoutExpired:
            return PlainTextResponse(f"Command timed out.", status_code=504)
        except FileNotFoundError:
            return PlainTextResponse(f"Command not found: {ts}", status_code=500)


@app.post("/relaunch/{ts_id}")
async def relaunch_task(
    ts_id: str,
    authorized: Union[User, str] = Depends(get_current_user_or_service),
):
    """Re-queue a finished task using its original JSON file."""
    if authorized != "service":
        if not isinstance(authorized, User) or "admin" not in authorized.groups:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin group required to relaunch analyses",
            )

    # Retrieve the analysis name from ts -i
    info_result = subprocess.run(
        f"{ts_env} {ts} -i {ts_id}",
        shell=True,
        capture_output=True,
        text=True,
        check=False,
        timeout=10,
    )
    name_match = re.search(r"--name\s+(\S+)", info_result.stdout)
    if not name_match:
        raise HTTPException(
            status_code=404, detail="Could not determine task name from ts info."
        )
    analysis_name = name_match.group(1)

    json_file = os.path.join(docker_stark_api_log_folder, f"{analysis_name}.json")
    if not os.path.isfile(json_file):
        raise HTTPException(status_code=404, detail=f"JSON file not found: {json_file}")

    try:
        with open(json_file, "r") as f:
            json_input = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        raise HTTPException(status_code=500, detail=f"Error reading JSON file: {e}")

    try:
        analysisIDNAME = queue_analysis(json_input)
        return PlainTextResponse(content=analysisIDNAME, status_code=200)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Relaunch failed: {e}")


if __name__ == '__main__':
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
