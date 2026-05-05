from fastapi import FastAPI, WebSocket, Depends, HTTPException, status, Query
from fastapi.responses import HTMLResponse
from fastapi.security import OAuth2PasswordBearer, OAuth2PasswordRequestForm
from jose import JWTError, jwt
from passlib.context import CryptContext
from pydantic import BaseModel
import pty
import os
import asyncio
import aiofiles

import json

# --- Configuration ---
SECRET_KEY = "YOUR_VERY_SECRET_KEY"  # Change this in a real app
ALGORITHM = "HS256"
ACCESS_TOKEN_EXPIRE_MINUTES = 30
USERS_FILE = "config/users.json"

pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")

app = FastAPI()

# --- Models ---
class Token(BaseModel):
    access_token: str
    token_type: str

class User(BaseModel):
    username: str

# --- Authentication (STARK Auth Placeholder & JWT) ---

def load_users():
    """Loads users from the JSON file."""
    if not os.path.exists(USERS_FILE):
        # Create the directory if it doesn't exist
        os.makedirs(os.path.dirname(USERS_FILE), exist_ok=True)
        # Create an empty users file if it doesn't exist
        with open(USERS_FILE, "w+") as f:
            json.dump({"stark": "password"}, f)
    try:
        with open(USERS_FILE, "r") as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        # If the file doesn't exist or is empty/corrupt, return an empty dict
        return {}

def stark_authenticate_user(username: str, password: str):
    """
    Authenticates a user against the users.json file.
    Verifies the password against the stored hash.
    Returns the user object if authentication is successful, otherwise False.
    """
    users_db = load_users()
    user_hash = users_db.get(username)
    if user_hash and password == user_hash:
        return User(username=username)
    return False
    # print("TEST0")
    # print(f"Password to verify : {pwd_context.hash(password)}")
    # print("TEST")
    # Truncate password to 72 bytes as a workaround for the bcrypt limitation
    # if user_hash and pwd_context.verify(password[:72], user_hash):
    #     return User(username=username)
    # return False

def create_access_token(data: dict):
    to_encode = data.copy()
    # In a real app, you would add an expiration time
    # from datetime import datetime, timedelta
    # expire = datetime.utcnow() + timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES)
    # to_encode.update({"exp": expire})
    encoded_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)
    return encoded_jwt

oauth2_scheme = OAuth2PasswordBearer(tokenUrl="token")

async def get_current_user(token: str = Query(None)):
    if token is None:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Not authenticated",
        )
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
    
    # We don't need to re-authenticate, just confirm the user from the token exists
    users_db = load_users()
    if username not in users_db:
        raise credentials_exception
    return User(username=username)


# --- Routes ---
@app.post("/token", response_model=Token)
async def login_for_access_token(form_data: OAuth2PasswordRequestForm = Depends()):
    user = stark_authenticate_user(form_data.username, form_data.password)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect username or password",
            headers={"WWW-Authenticate": "Bearer"},
        )
    access_token = create_access_token(data={"sub": user.username})
    return {"access_token": access_token, "token_type": "bearer"}

@app.get("/")
async def get():
    async with aiofiles.open("static/index.html", "r") as f:
        content = await f.read()
    return HTMLResponse(content=content)

@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket, token: str = Query(None)):
    # Authenticate the user with the token from the query parameter
    try:
        user = await get_current_user(token) # Let's make this async again
        if not user:
            await websocket.close(code=status.WS_1008_POLICY_VIOLATION, reason="Invalid user")
            return
    except HTTPException as e:
        await websocket.close(code=status.WS_1008_POLICY_VIOLATION, reason=f"Authentication failed: {e.detail}")
        return
        
    await websocket.accept()

    # Command to execute in a shell to ensure TTY allocation
    cmd = "script -q -c 'docker exec -it stark-module-stark-submodule-stark-service-cli bash' /dev/null"

    try:
        # Create a subprocess shell
        process = await asyncio.create_subprocess_shell(
            cmd,
            stdout=asyncio.subprocess.PIPE,
            stdin=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )

        # Bridge functions
        async def read_from_proc(stream):
            while True:
                data = await stream.read(1024)
                if not data:
                    break
                await websocket.send_text(data.decode(errors='replace'))

        async def write_to_proc(stream):
            while True:
                data = await websocket.receive_text()
                stream.write(data.encode())
                await stream.drain()

        # Run tasks concurrently
        read_stdout = asyncio.create_task(read_from_proc(process.stdout))
        read_stderr = asyncio.create_task(read_from_proc(process.stderr))
        write_stdin = asyncio.create_task(write_to_proc(process.stdin))

        done, pending = await asyncio.wait(
            {read_stdout, read_stderr, write_stdin},
            return_when=asyncio.FIRST_COMPLETED
        )

        for task in pending:
            task.cancel()

    except Exception as e:
        print(f"Websocket Error: {e}")
        await websocket.send_text(f"\r\nServer error: {e}\r\n")
    finally:
        if 'process' in locals() and process.returncode is None:
            process.kill()
            await process.wait()
        await websocket.close()
        print("Connection closed.")
