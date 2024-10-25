README - Last update: 20241025

THIS IS THE RELEASE BRANCH FOR PAPER:
# On Density-based Local Community Search

Depending on your objective, follow steps as below:
```mermaid
flowchart LR
    A>"Getting Started"]
    B>"(Optional) Change LP Solver"]
    C1>"Run One Instance"]
    C2>"Run Case Study"]
    C3>"Run Experiment"]
    A --> B
    B --> C1
    B --> C2
    B --> C3
```

---
```mermaid
flowchart LR
    A["Getting Started"]
```

### Prerequisites
Julia

### Packages required
Installing packages in Julia for example:
```julia
using Pkg
Pkg.add("MatrixNetworks")
Pkg.add("MAT")
Pkg.add("StatsBase")
Pkg.add("JuMP")
Pkg.add("HiGHS")  # Skip if using other LP Solver
```

### Working folder
For all tasks, open cmd/bash, navigate to the repository's root directory, then enter the `./src` folder and run `julia`

---
```mermaid
flowchart LR
    B["(Optional) Change LP Solver"]
```
By default, the LP (Linear-Programming) solver `HiGHS` would be (installed and) used.
If you want to use another LP solver, take `CPLEX` for example, you need to refer to:
- Install [IBM ILOG CPLEX Optimization Studio](https://www.ibm.com/products/ilog-cplex-optimization-studio). You will need to have/obtain a license
- Install [CPLEX for Julia](https://www.ibm.com/products/ilog-cplex-optimization-studio). Including `Pkg.add("CPLEX")`.
- In the file `.\src\LP_load_solver.jl`, comment the blocks using HiGHS and uncomment the blocks using CPLEX.

---
```mermaid
flowchart LR
    C1["Run One Instance"]
```

```julia
include("LP_algorithm.jl")
```

### Read/create a graph
Load a toy graph in `./Example_small/` (there are some other toy graphs in the same folder for exploration):

```julia
A = readIN("lobster.in", "../Example_small/")
```

This graph is same as:
```julia
A = sparse([1,1,1,2,2,3,3,4,2,3,4,3,4,4,5,5], [2,3,4,3,4,4,5,5,1,1,1,2,2,3,3,4], ones(Float64, 16), 5, 5)
```

Run GADS under graph `A` with $R = [1,2], x (\omega_{12}) = 1, y (\omega_{24}) = -1$ (see Section 1.4):
```julia
R = [1,2]
x = 1
y = -1
# 7 Numbers in weight corresponds to the Weight Configuration with 10 edge weights as Definition 1 (in paper)
# minus 3 edge weights excluded by Definition 4-C2 (edge weights of edges between V_3 and V_4).
# Note: GADS algorithm only guarantee to work on 0 <= x <= 2, y <= 0 and all other weights same as below.
weight = [2,x,0,0,0,0,y]  # Weight Configuration Ω
SolveLPAnchoredDensestSubgraphGeneric(A, R, weight)
```

Say the output is `(densestSubgraph(1.2, [1, 2, 3, 4, 5]), 0.001)`, it means
- $S^*_{\Omega, R} = [1, 2, 3, 4, 5]$, the Local Densest Graph of G under weight configuration $\Omega$ and seed set $R$;
- $\rho^*_{\Omega, R} = 1.2$, the local density of the above graph;
- Runtime is $0.001$ seconds.


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

This code is a fork of [HypergraphFlowClustering](https://github.com/nveldt/HypergraphFlowClustering) by [Nate Veldt](https://github.com/nveldt). We are grateful for [Nate Veldt]'s contributions, such as:
- `maxflow.jl` with modifications.
