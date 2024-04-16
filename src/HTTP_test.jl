using HTTP
using JSON3

function fetch_graph_data()
    url = "http://localhost:7687/db/data/transaction/commit"  # Adjust port if needed
    auth = base64("username:password")  # Replace with your Neo4j credentials
    headers = [
        "Content-Type" => "application/json",
        "Authorization" => "Basic $auth"
    ]
    body = JSON3.write({
        "statements" => [{
            "statement" => "MATCH (n)-[r]->(m) RETURN id(n), id(m)",
            "resultDataContents" => ["row"]
        }]
    })

    response = HTTP.post(url, headers, body)
    if response.status == 200
        return JSON3.read(response.body)["results"][1]["data"]
    else
        error("Failed to fetch data: $(response.status) $(String(response.body))")
    end
end