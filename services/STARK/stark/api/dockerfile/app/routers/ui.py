#!/usr/bin/env python

import os

from fastapi import APIRouter, Request  # pyright: ignore[reportMissingImports]
from fastapi.responses import HTMLResponse  # pyright: ignore[reportMissingImports]
from fastapi.templating import Jinja2Templates  # pyright: ignore[reportMissingImports]

from config import refresh_interval_ms

router = APIRouter()
templates = Jinja2Templates(directory="templates")


@router.get("/", response_class=HTMLResponse)
async def read_root(request: Request):
    version = str(int(os.path.getmtime("static/script.js")))
    return templates.TemplateResponse(
        "index.html",
        {
            "request": request,
            "version": version,
            "refresh_interval_ms": refresh_interval_ms,
        },
    )
