var config = {
    container_id: "viz",
    server_url: "bolt://localhost:7687",
    server_user: "neo4j",
    server_password: "12345678",
    initial_cypher: "MATCH (n)-[r]->(m) RETURN n, r, m",
};

var viz = new NeoVis.default(config);
viz.render();