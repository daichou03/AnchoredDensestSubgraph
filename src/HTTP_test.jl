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
    # Function to execute a query and print results
    try
        result = session.run(query)
        for record in result
            println(record)
        end
    finally
        #session.close()
    end
end

# As unweighted, undirected
function fetch_graph_data_as_CSC(session, attrib="id")
    # This query should return pairs of node IDs representing edges
    result = session.run("MATCH (n)-[r]->(m) RETURN n.$(attrib) as id1, m.$(attrib) as id2")
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