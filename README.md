README - Last update: 20210818

# Anchored Densest Subgraph

------

## Getting Started

### Prerequisites
(Linux) Install HDF5 if you don't already have.

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
Under `./src`, enter `julia`.

Load all packages needed for testing:
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

The loaded data graph (`A`) is in `SparseMatrixCSC` format.

#### Accepted input graph by `readIN()`
This algorithm can work on other data graphs you downloaded.  
All our experimental data graphs can be found at [SNAP](https://snap.stanford.edu/data/index.html), [KONECT](http://konect.cc/) and [Network Repository](https://networkrepository.com/networks.php).

Once you downloaded  (for example, [uk2007](http://konect.cc/networks/dimacs10-uk-2007-05/) (large!)) and extracted the data,  
`readIN()` accepts file of the following format:

- The first line is the number of vertices and the number of edges respectively.
- The remaining lines are the list of edges (unweighted, undirected, unidirectional), vertices are 1-indexed.

Example:  
```
5 8  
1 2  
1 3  
1 4  
2 3  
2 4  
3 4  
3 5  
4 5  
```

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

- Originally fork of https://github.com/nveldt/HypergraphFlowClustering, using `maxflow.jl` with modifications.

------
## Troubleshooting

Case 1 - In Linux, when installing Julia packages, found HDF5 related errors:  
This is due to either HDF5 is not configured properly. Uninstall, install a lower version to force update dependency, uninstall, then install the newest version again can fix.  
In Linux:
```
apt-get -u install hdf5-tools
```

If still getting `UndefVarError: libhdf5 not defined`:  
In Julia:  
```
using Pkg  
Pkg.rm("HDF5")  
Pkg.add(Pkg.PackageSpec(;name="HDF5",version="0.11.1"))  
Pkg.build("HDF5")  
Pkg.rm("HDF5")  
Pkg.add("HDF5") 
```
