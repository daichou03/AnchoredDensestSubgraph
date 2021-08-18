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
This algorithm can work on other data graphs you downloaded, for example, [uk2007](http://konect.cc/networks/dimacs10-uk-2007-05/) (large!).  
All our experimental data graphs can be found at [SNAP](https://snap.stanford.edu/data/index.html), [KONECT](http://konect.cc/) and [Network Repository](https://networkrepository.com/networks.php).

In our experiment, we preprocessed all the data graphs we used so that they can be loaded by `readIN()` (see above for syntax):  
- The first line is the number of vertices and the number of edges respectively.  
- The remaining lines are the list of edges (unweighted, undirected), vertices are 1-indexed.  
- A line starting with `#` and `%` is ignored.

You may need to preprocess the raw data first if it does not meet the above conditions.

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

Alternatively, there is a `readRaw()` that does the same as `readIN()` but assumes the file does not have a header line for # vertices and # edges.  
You pass the number of vertices and the number of edges as the parameters:

```julia
A = readRaw("zebra.txt", 27, 111, "../Example_raw")
```

This may come in handy when you downloaded a file with an unweighted, undirected edge list without the need of writing out a new file or modify the original file.

#### Find densest Graph

Find anchored densest subgraph of `A` with reference vertices `R` and its anchored density:

```julia
R = [1,2]
LocalAnchoredDensestSubgraph(A, R) # The strongly-local implementation
```
Note `R` is a 1-indexed vector of indices of vertices in `A`.

------
#### Acknowledgement (for code)

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
