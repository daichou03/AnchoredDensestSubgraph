README - Last update: 20240218

# Generalized Anchored Densest Subgraph

------

## Getting Started

### Prerequisites
(Linux) Install HDF5 if you don't already have it.

### Install 64-bit Julia

### Packages required
In Julia:
```julia
using Pkg
Pkg.add("MatrixNetworks")
Pkg.add("MAT")
Pkg.add("StatsBase")
```

### Run 1 instance
Under `./src`, enter `julia`.

Import this to load the core algorithm:
```julia
include("LP_algorithm.jl")
```

Read graph file:
Some toy data graphs are in /Example_small/, for example:

```julia
A = readIN("lobster.in", "../Example_small/")
```

Loading this graph is equivalent to:
```julia
A = sparse([1,1,1,2,2,3,3,4,2,3,4,3,4,4,5,5], [2,3,4,3,4,4,5,5,1,1,1,2,2,3,3,4], ones(Float64, 16), 5, 5)
```

Run GADS under graph `A` with `R = [1,2]`, x $(\omega_{12}) = 1$, y $(\omega_{24}) = -1$ (see Section 1.4):
```julia
R = [1,2]
x = 1
y = -1
# Each number in this array corresponds to a specific edge weight.
# Due to the limitation of the LP algorithm, changing edge weights other than x and y does not (always) guarantee correctness.
weight = [2,x,0,0,0,0,y]
SolveLPAnchoredDensestSubgraphGeneric(A, R, weight)
```
You have just run the LP algorithm introduced in the original paper!


Some small (compared to other real-world data graphs), preprocessed real-world data graphs are in /Example_SCC/. For example:

```julia
A = readIN("eucore.in", "../Example_SCC/")
```

The loaded data graph (`A`) is in `SparseMatrixCSC` format.

To look for codes that mass (re)produce the experimental results:
```julia
include("LP_compare_test.jl")
```

#### Accepted input graph by `readIN()`
This algorithm can work on other data graphs you downloaded, for example, [uk2007](http://konect.cc/networks/dimacs10-uk-2007-05/) (large!).  
All our experimental data graphs can be found at [SNAP](https://snap.stanford.edu/data/index.html), [KONECT](http://konect.cc/) and [Network Repository](https://networkrepository.com/networks.php).

In our experiment, we preprocessed all the data graphs we used so that they can be loaded by `readIN()` (see above for syntax):  
- The first line contains two integers - the number of vertices (**n**) and the number of edges (**m**), respectively.
- For the next **m** lines, each line contains two integers representing an (unweighted, undirected) edge; all vertices are 1-indexed.
- Any contents after these **m** lines will not be processed; Any line starts with `#` and `%` is ignored. 

If it does not meet the above conditions, you may need to preprocess the raw data first.

Example of a legitimate .in file:  
```
% "Lobster" toy graph
# 4-clique:
5 8  
1 2  
1 3  
1 4  
2 3  
2 4  
3 4  
# + An extra extruding triangle:
3 5  
4 5  

This line is also a comment.
```

Alternatively, a `readRaw()` does the same as `readIN()` but assumes the file does not have a header line for # vertices and # edges.  
Instead, you pass the number of vertices and the number of edges as the parameters:

```julia
A = readRaw("zebra.txt", 27, 111, "../Example_raw")
```

------
#### Acknowledgments

This project is a fork of [HypergraphFlowClustering](https://github.com/nveldt/HypergraphFlowClustering) by [Nate Veldt](https://github.com/nveldt). We are grateful for [Nate Veldt]'s contributions, such as:
- `maxflow.jl` with modifications.
