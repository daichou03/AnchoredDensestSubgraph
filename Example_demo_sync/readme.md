Data graphs under this folder is also in *.in format.
However, they are also used by the demo project in a "synchronized" mode that corresponds to a graph database in neo4j.
That is, you need to make a *.in file and modify the graph database in neo4j such that:
- Like other *.in, the *.in file (only) has information for nodes and edges, and are undirected.
- Nodes in the corresponding neo4j database needs to have a "julia_id" property that aligns with the "id" of the nodes for this *.in, 1-indexed. Of course, edges in neo4j database and in *.in file must match.

*** Setup julia_id property
For example, you have julia_id.csv like this:
```csv
name,julia_id
Addam-Marbrand,1
Aegon-Frey-(son-of-Stevron),2
Aegon-I-Targaryen,3
Aegon-Targaryen-(son-of-Rhaegar),4
Aegon-V-Targaryen,5
```

Copy this file to import folder of the neo4j database, and run this query:
```cypher
LOAD CSV WITH HEADERS FROM 'file:///julia_id.csv' AS row
MATCH (n)
WHERE n.name = row.name
SET n.julia_id = toInteger(row.julia_id)
```

