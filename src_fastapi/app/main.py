from fastapi import FastAPI, Request, HTTPException, Form
from fastapi.responses import HTMLResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from neo4j import GraphDatabase, basic_auth
from fastapi.middleware.cors import CORSMiddleware


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
app.mount("/static", StaticFiles(directory="static"), name="static")

# Set up templates
templates = Jinja2Templates(directory="templates")

@app.get("/")
async def serve_home(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})

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
