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
templates = Jinja2Templates(directory="app/templates")

@app.get("/")
async def serve_home(request: Request):
    return templates.TemplateResponse({"request": request}, "index.html")

@app.get("/")
async def get_index():
    return FileResponse('templates/index.html')

@app.get("/view-graph/")
async def view_graph(request: Request):
    return templates.TemplateResponse("graph.html", {"request": request})

@app.get("/specimen-1/")
async def test_html(request: Request):
    return templates.TemplateResponse("specimen-1.html", {"request": request})

# @app.get("/test/")
# async def test_html(request: Request):
#     return templates.TemplateResponse("test.html", {"request": request})


@app.post("/execute-query/")
async def execute_query(database_name: str = Form(...), query: str = Form(...)):
    if not database_name:
        raise HTTPException(status_code=400, detail="Database name must not be empty")
    if not query:
        raise HTTPException(status_code=400, detail="Query must not be empty")
    uri = "bolt://localhost:7687"
    username = HTTP_CREDS_USERNAME
    password = HTTP_CREDS_PASSWORD
    driver = GraphDatabase.driver(uri, auth=basic_auth(username, password))

    try:
        with driver.session(database=database_name) as session:
            result = session.run(query)
            return result.data()
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
    finally:
        driver.close()
    

@app.post("/get-node-count/")
async def get_node_count(database_name: str = Form(...)):
    data = await execute_query(database_name, "MATCH (n) RETURN count(n) as node_count")
    print(data)
    return {"database": database_name, "node_count": data[0]["node_count"]}


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

