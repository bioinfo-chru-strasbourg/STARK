#!/usr/bin/env python

from fastapi import FastAPI  # pyright: ignore[reportMissingImports]
from fastapi.staticfiles import StaticFiles  # pyright: ignore[reportMissingImports]
from fastapi.responses import FileResponse  # pyright: ignore[reportMissingImports]

from routers import auth, ui, analysis, cluster, archive
from routers import queue as queue_router
from config import local_port

app = FastAPI()
app.mount("/static", StaticFiles(directory="static"), name="static")


# Temporary route for favicon
@app.get("/favicon.ico", include_in_schema=False)
async def favicon():
    return FileResponse("static/favicon.ico")


app.include_router(auth.router)
app.include_router(ui.router)
app.include_router(cluster.router)
app.include_router(analysis.router)
app.include_router(queue_router.router)
app.include_router(archive.router)


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=int(local_port))
