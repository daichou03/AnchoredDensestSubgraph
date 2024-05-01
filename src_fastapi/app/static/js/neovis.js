document.addEventListener('DOMContentLoaded', function() {
    var config = {
        container_id: "viz",
        server_url: "bolt://localhost:7687",
        //server_user: "neovisuser",
        //server_password: "12345678",
        serverDatabase: "zachary",
        initial_cypher: "MATCH (n)-[r]->(m) RETURN n, r, m limit 100",
        labels: {
            // Define labels and their properties for visualization
        },
        relationships: {
            // Define relationship types and their properties
        },
    };

    var viz = new NeoVis.default(config);
    viz.render();
});