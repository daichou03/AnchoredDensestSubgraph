using PyCall
using SparseArrays
include("HTTP_creds.jl")

neo4j = pyimport("neo4j")

function GetNeo4jSession(database_name, username=HTTP_CREDS_USERNAME, password=HTTP_CREDS_PASSWORD)
    uri = "bolt://localhost:7687"

    driver = neo4j.GraphDatabase.driver(uri, auth=(username, password))
    return driver.session(database = database_name)
end

# query = "MATCH (n)-[r]->(m) RETURN id(n), id(m)"
# run_query(driver, query)


function run_query(session, query)
    result = session.run(query)
    for record in result
        println(record)
    end
end


# Assuming all nodes in neo4j has a key attribute called id that is 1-indexed
# As unweighted, undirected
function fetch_graph_data_as_CSC(session)
    # This query should return pairs of node IDs representing edges
    result = session.run("MATCH (n)-[r]->(m) RETURN n.id as id1, m.id as id2")
    edges = []
    for record in result
        push!(edges, record)  # +1 because Julia is 1-indexed
    end
    return edges_to_sparse_matrix(edges)
end

function edges_to_sparse_matrix(edges)
    N = maximum([max(edge...) for edge in edges])
    M = length(edges)
    ei = zeros(Int64, M*2)
    ej = zeros(Int64, M*2)
    count = 0
    for (v1, v2) in edges
        if v1 == v2
            continue
        elseif v1 > v2
            v1, v2 = v2, v1
        end
        count += 1
        ei[count*2-1] = v1
        ej[count*2-1] = v2
        ei[count*2] = v2
        ej[count*2] = v1
    end
    A = sparse(ei[1:count*2], ej[1:count*2], ones(Float64, count*2), N, N)
    return A
end

# To test
function send_ids_to_neo4j(session, ids)
    # Define a Cypher query that uses the list of IDs
    query = """
    UNWIND $ids AS id
    MATCH (n) WHERE n.id = id
    RETURN n
    """
    # Execute the query passing the list of IDs as a parameter
    result = session.run(query, Dict("ids" => ids))
end

# function fetch_nodes_by_ids(session, ids)
#     # Define a Cypher query that fetches nodes based on a list of IDs
#     query = """
#     UNWIND $ids AS node_id
#     MATCH (n)
#     WHERE id(n) = node_id
#     RETURN id(n) AS node_id
#     """
#     # Execute the query passing the list of IDs as a parameter
#     result = session.run(query, Dict("ids" => ids))
#     fetched_ids = Int64[]  # Create an empty array to store node IDs

#     # Iterate over the results and populate the array
#     for record in result
#         push!(fetched_ids, record["node_id"])
#     end
#     return fetched_ids
# end