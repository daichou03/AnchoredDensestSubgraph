from fastapi import FastAPI, Request, HTTPException, Form, Body
from fastapi.responses import HTMLResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from neo4j import GraphDatabase, basic_auth
from fastapi.middleware.cors import CORSMiddleware
import httpx


HTTP_CREDS_USERNAME = "neo4j"
HTTP_CREDS_PASSWORD = "12345678"

app = FastAPI()

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # You can specify domains for security
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files
app.mount("/static", StaticFiles(directory="app/static"), name="static")

# Set up templates
templates = Jinja2Templates(directory="templates")

@app.get("/")
async def serve_home(request: Request):
    return templates.TemplateResponse({"request": request}, "index.html")

@app.get("/")
async def get_index():
    return FileResponse('templates/index.html')


@app.post("/get-node-count/")
async def get_node_count(database_name: str = Form(None)):
    if not database_name:
        raise HTTPException(status_code=400, detail="Database name must not be empty")
    print(database_name)
    uri = "bolt://localhost:7687"  # Adjust as necessary
    username = HTTP_CREDS_USERNAME  # Default username, replace with your credentials
    password = HTTP_CREDS_PASSWORD  # Replace with your password
    driver = GraphDatabase.driver(uri, auth=basic_auth(username, password))

    try:
        with driver.session(database=database_name) as session:
            query = "MATCH (n) RETURN count(n) as node_count"
            result = session.run(query)
            node_count = result.single()["node_count"]
        return {"database": database_name, "node_count": node_count}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
    finally:
        driver.close()


@app.post("/load-graph/")
async def load_graph(filename: str = Body(...), dir: str = Body(...)):
    try:
        response = httpx.post(
            "http://localhost:8080",
            json={"action": "load-graph", "filename": filename, "dir": dir}
        )
        response.raise_for_status()
        return {"message": response.text}
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/solve-gads/")
async def solve_gads(R: list = Body(...), weight_map: list = Body(...)):
    try:
        response = httpx.post("http://localhost:8080", json={"action": "solve-gads", "R": R, "weight_map": weight_map})
        response.raise_for_status()
        return response.json()
    except httpx.HTTPStatusError as e:
        raise HTTPException(status_code=e.response.status_code, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

