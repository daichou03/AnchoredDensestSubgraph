README - Last update: 20201223, daichou03

Local Densest Subgraph

(Currently as a fork of HypergraphFlowClustering)

------

# Getting Started

## Prerequisites
(Linux) Install HDF5 if you don't already have.

> Case 1 - Ubuntu: HDF5 not configured properly at the very start, uninstall, install a lower version to force update dependency, uninstall, then install the newest version again can fix.
> `apt-get -u install hdf5-tools`
> Still getting `UndefVarError: libhdf5 not defined`
> In Julia:
> using Pkg
> Pkg.rm("HDF5")
> Pkg.add(Pkg.PackageSpec(;name="HDF5",version="0.11.1"))
> Pkg.build("HDF5")
> Pkg.rm("HDF5")
> Pkg.add("HDF5")

## Install 64-bit Julia

## Install packages
In Julia:
```julia
using Pkg
Pkg.add("MatrixNetworks")
Pkg.add("MAT")
Pkg.add("StatsBase")
```
* MatrixNetworks
* TODO: Need to check.

## Testing
Under HypergraphFlowClustering/src, enter julia.
```julia
include("Test_yd.jl")
```

Read graph file:
For loading unweighted, undirected graph:

A = readIN("lobster.in", "../Example_small")

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
GlobalAnchoredDensestSubgraph(A, R)

Note R is a vector of indices of vertices in A, 1-indexed.

------
Acknowledgement (for code)

- Fork of https://github.com/nveldt/HypergraphFlowClustering
- Currently using maxflow.jl, will use more
- May use/modify HyperLocal.jl
