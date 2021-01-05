README - Last update: 20201223, daichou03

Local Densest Subgraph

(Currently as a fork of HypergraphFlowClustering)

------

Documentation

Install Julia

Package to install:
MatrixNetworks
TODO: Need to check.

Read graph file:
For loading unweighted, undirected graph:

A = readIN("lobster.in")

- The first line is the number of vertices and the number of edges respectively
- The remaining lines should be the edge list, vertices are 1-indexed.

Example:
5 8
1 2
1 3
1 4
2 3
2 4
3 4
3 5
4 5

For loading weighted, directed graph:

B = readSMAT("graph.smat")

Find (global) densest subgraph and its density:

GlobalMaximumDensity(A)

TODO: which algorithm?
The subgraph is a maximum densest subgraph.

Find local densest subgraph and its density with reference vertices:

R = vec([1 2])
LocalMaximumDensity(A, R)

Note R is a vector of indices of vertices in A, 1-indexed.

------
Acknowledgement (for code)

- Fork of https://github.com/nveldt/HypergraphFlowClustering
- Currently using maxflow.jl, will use more
- May use/modify HyperLocal.jl
