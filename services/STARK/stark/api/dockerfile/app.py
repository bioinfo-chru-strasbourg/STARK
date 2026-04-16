#!/usr/bin/env python

from fastapi import FastAPI  # pyright: ignore[reportMissingImports]
from fastapi.staticfiles import StaticFiles  # pyright: ignore[reportMissingImports]

from routers import auth, ui, analysis, cluster
from routers import queue as queue_router

app = FastAPI()
app.mount("/static", StaticFiles(directory="static"), name="static")

app.include_router(auth.router)
app.include_router(ui.router)
app.include_router(cluster.router)
app.include_router(analysis.router)
app.include_router(queue_router.router)


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
