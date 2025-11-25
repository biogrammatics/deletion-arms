"""
FastAPI application for the Deletion Arms Designer
"""

import os
from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.middleware.cors import CORSMiddleware

from web.routes import router

# Get the directory containing this file
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Create FastAPI app
app = FastAPI(
    title="Deletion Arms Designer",
    description="Design gene knockout constructs using half-site restriction enzyme strategy",
    version="1.0.0"
)

# Add CORS middleware for future Rails integration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure appropriately for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files
app.mount("/static", StaticFiles(directory=os.path.join(BASE_DIR, "static")), name="static")

# Setup templates
templates = Jinja2Templates(directory=os.path.join(BASE_DIR, "templates"))

# Include API routes
app.include_router(router)


@app.get("/")
async def home(request: Request):
    """Serve the main web UI"""
    return templates.TemplateResponse("index.html", {"request": request})


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
