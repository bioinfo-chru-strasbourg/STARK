#!/usr/bin/env python

import json
import os
from typing import Union

from fastapi import Depends, HTTPException, Security, status  # pyright: ignore[reportMissingImports]
from fastapi.security import (  # pyright: ignore[reportMissingImports]
    APIKeyHeader,
    OAuth2PasswordBearer,
)
from jose import JWTError, jwt  # pyright: ignore[reportMissingModuleSource, reportMissingImports]
from passlib.context import CryptContext  # pyright: ignore[reportMissingModuleSource, reportMissingImports]

from config import ALGORITHM, SECRET_KEY, STARK_API_KEY, USERS_FILE
from models import User

pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")
api_key_header = APIKeyHeader(name="X-API-Key", auto_error=False)
oauth2_scheme = OAuth2PasswordBearer(tokenUrl="token", auto_error=False)


def get_api_key(api_key: str = Security(api_key_header)):
    if api_key == STARK_API_KEY:
        return api_key
    return None


def load_users() -> dict:
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
    users_db = load_users()
    if username not in users_db:
        return False
    entry = users_db[username]
    if isinstance(entry, dict):
        if entry.get("password") == password:
            return User(username=username, groups=entry.get("groups", []))
    elif entry == password:
        return User(username=username, groups=[])
    return False


def create_access_token(data: dict) -> str:
    to_encode = data.copy()
    encoded_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)
    return encoded_jwt


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
    api_key: str = Security(api_key_header),
    user: User = Depends(get_current_user),
) -> Union[User, str]:
    if api_key and api_key == STARK_API_KEY:
        return "service"
    if user:
        return user
    raise HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED, detail="Not authenticated"
    )
