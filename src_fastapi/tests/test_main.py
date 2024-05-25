from fastapi.testclient import TestClient
from ..app.main import app

client = TestClient(app)

# Sanity
def test_get():
    response = client.get("/")
    assert response.status_code == 200

# FastAPI -> Julia server
def test_load_graph_sync():
    response = client.post("/load-graph-sync/", json="got")
    assert response.status_code == 200
    assert response.json() == {"message": "Synchronized data got loaded successfully"}

def test_load_graph():
    response = client.post("/load-graph/", json={"filename": "pincer.in", "dir": "../Example_small/"})
    assert response.status_code == 200
    assert response.json() == {"message": "Graph loaded successfully from: ../Example_small/pincer.in"}

def test_solve_gads():
    response = client.post("/solve-gads/", json={"R": [1,2], "weight_map": [2,2,0,0,0,0,-1]})
    assert response.status_code == 200
    assert response.json() == [1,2,3,4,5]
