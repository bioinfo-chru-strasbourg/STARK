#!/usr/bin/env python

from pydantic import BaseModel  # pyright: ignore[reportMissingImports]
from typing import List


class Token(BaseModel):
    access_token: str
    token_type: str


class User(BaseModel):
    username: str
    groups: List[str] = []
