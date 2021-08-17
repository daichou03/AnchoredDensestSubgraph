README - Last update: 20210704

# Anchored Densest Subgraph

------

## Getting Started

### Prerequisites
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

### Install 64-bit Julia

### Packages required
In Julia:
```julia
using Pkg
Pkg.add("MatrixNetworks")
Pkg.add("MAT")
Pkg.add("StatsBase")
```

### Testing
Under ./src, enter julia.
```julia
include("Query_test_yd.jl")
```

Read graph file:  
Some toy data graphs are in /Example_small/:

```julia
A = readIN("lobster.in", "../Example_small")
```

Some small, preprocessed real-world data graphs are in /Example_SCC/:

```julia
A = readIN("eucore.in", "../Example_SCC/")
```

#### Accepted input graph
This algorithm can work on any data graphs you downloaded.  
The graph needs to be preprocessed as unweighted, undirected if it isn't.  
`readIN()` accepts the following format:

- The first line is the number of vertices and the number of edges respectively.
- The remaining lines are the list of edges (only 1 direction is needed), vertices are 1-indexed.

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

Find densest subgraph (by global, common definition) of `A` and its density:

```julia
GlobalDensestSubgraph(A)
```

Find local densest subgraph of `A` and its density with reference vertices:

```julia
R = [1,2]
LocalAnchoredDensestSubgraph(A, R)
```
Note `R` is a 1-indexed vector of indices of vertices in `A`.

------
Acknowledgement (for code)

- Fork of https://github.com/nveldt/HypergraphFlowClustering
- Using maxflow.jl
