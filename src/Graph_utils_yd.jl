using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase # TODO: To install
using Random
using Base
include("maxflow.jl")
include("Helper_io.jl")

#------------
# Graph Utils
#------------

# ----------------------------------
# Laplacians package functions start
# ----------------------------------
# Credit:
# https://github.com/danspielman/Laplacians.jl/blob/master/src/graphUtils.jl
# Too heavy to load the entire package, copy the specific function here instead.

# Unweighted degree.
GetDegree(mat::SparseMatrixCSC{Tv,Ti}, v::Ti) where {Tv,Ti} = mat.colptr[v+1]-mat.colptr[v]

# --------------------------------
# Laplacians package functions end
# --------------------------------


# Convert unweighted graph to transition graph
function toTransitionGraph(B::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
    P = copy(B)
    N = size(P, 1)
    for v = 1:N
        vdeg = GetDegree(P, v)
        if vdeg > 0
            for pt = P.colptr[v] : P.colptr[v+1] - 1
                P.nzval[pt] = 1 / vdeg
            end
        end
    end
    return P
end


# YD 20200201: https://github.com/JuliaLang/julia/blob/master/stdlib/SparseArrays/src/sparsevector.jl
function GetAdjacency(B::SparseMatrixCSC, V::Int64, Self::Bool=true)
    L = SparseArrays.nonzeroinds(B[:,V])
    if Self
        L = insert_and_dedup!(L, V) # Keep this list sorted
    end
    return L
end

# Imitating a user try to submit a list of people including themselves and some of their friends (neigbhours), of #Size people in total (including themselves).
# If the user don't have that many friends, submit all their friends + self instead.
function GetRandomAdjacency(B::SparseMatrixCSC, V::Int64, Size::Int64)
    L = GetAdjacency(B, V, false)
    if length(L) >= Size
        L = StatsBase.sample(L, Size - 1, replace=false, ordered=true)
    end
    L = insert_and_dedup!(L, V)
    return L
end

function GetRandomAdjacency(B::SparseMatrixCSC, Size::Int64)
    return GetRandomAdjacency(B, rand(1:size(B, 1)), Size)
end

# Use Set to rewrite this.
# YD 20210114: Depends on Size, get faster than V1 with larger size.
# Compared with V1:
# 120% time on |S| = 20 - 200
# 91% time on |S| = 2000
function GetSampleUntilSizeV2(B::SparseMatrixCSC, V::Int64, Size::Int64)
    r = [V]
    adj_set = Set(GetAdjacency(B, V, false))
    size = 1
    while size < Size && length(adj_set) > 0
        next = rand(adj_set)
        append!(r, next)
        adj_set = setdiff(union(adj_set, Set(GetAdjacency(B, next, true))), Set([next]))
        size += 1
    end
    return r
end

# Connected components
# The implementations below are not performant.
# Use Laplacians package - refer:
# https://github.com/danspielman/Laplacians.jl/blob/bbc21749131b5c6452e2fa5b4733f84129d4ab09/src/graphAlgs.jl#L157
# Note: loading the package itself takes some time.

function ExtractConnectedComponent(B::SparseMatrixCSC)
    explored = GetAdjacency(B,1,true)
    adj = setdiff(explored, [1])
    while length(adj) > 0
        adj = setdiff(GetComponentAdjacency(B, adj, false), explored)
        explored = union(explored, adj)
    end
    return explored
end

# Returns the number of connected components.
# If ReturnLargestCCIndex, returns the largest connected component's index INSTEAD (I know this is bad code, not necessary to make it better for now).
# Also if want to know the largest CC's index, definitely can stop early.
function DetectConnectedComponents(B::SparseMatrixCSC, ReturnLargestCCIndex::Bool=false, ShowLengthOfComponents::Bool=false)
    remaining = copy(B)
    components = 0
    largestCCIndex = 1
    largestCC = 0
    while size(remaining, 1) > 0
        nextCC = ExtractConnectedComponent(remaining)
        remaining_components = collect(setdiff(Set(1:size(remaining, 1)), nextCC))
        remaining = remaining[remaining_components, remaining_components]
        components += 1
        if components == 1 || length(largestCC) < length(nextCC)
            largestCC = nextCC
            largestCCIndex = components
        end
        if ShowLengthOfComponents
            println(string("Length of connected component #", components, ": ", length(nextCC)))
        end
    end
    if ReturnLargestCCIndex
        return largestCCIndex
    else
        return components
    end
end

# Laplacians.biggestComp
function RetrieveLargestConnectedComponent(B::SparseMatrixCSC)
    if @isdefined biggestComp
        biggestComp(B) # using Laplacians
    else
        YDRetrieveLargestConnectedComponent(B)
    end
end

function YDRetrieveLargestConnectedComponent(B::SparseMatrixCSC)
    println("WARNING: This LCC algorithm is not optimised and takes forever for large graphs! Use Laplacians.biggestComp() instead.")
    largestCCIndex = DetectConnectedComponents(B, true, false)
    remaining = copy(B)
    for i = 1:largestCCIndex-1
        nextCC = ExtractConnectedComponent(remaining)
        remaining_components = setdiff(1:size(remaining, 1), nextCC)
        remaining = remaining[remaining_components, remaining_components]
    end
    nextCC = ExtractConnectedComponent(remaining)
    return remaining[nextCC,nextCC]
end

function DetectConnectedComponents(B::SparseMatrixCSC, ShowLengthOfComponents::Bool=false)
    remaining = copy(B)
    components = 0
    while size(remaining, 1) > 0
        explored = Set(GetAdjacency(remaining,1,true))
        adj_set = setdiff(explored, Set([1]))
        while length(adj_set) > 0
            adj_set = setdiff(SetGetComponentAdjacency(remaining, collect(adj_set), false), explored)
            explored = union(explored, adj_set)
        end
        remaining_components = collect(setdiff(Set(1:size(remaining, 1)), explored))
        remaining = remaining[remaining_components, remaining_components]
        components += 1
        if ShowLengthOfComponents
            println(string("Length of connected component #", components, ": ", length(explored)))
        end
    end
    return components
end

# Note there are |S| convertions from array to set, and 1 conversion from set to array.
# 20210219: Slow on large sets.
function GetComponentAdjacency(B::SparseMatrixCSC, S::Vector{Int64}, Self::Bool=true)
    return collect(SetGetComponentAdjacency(B,S,Self))
end

function SetGetComponentAdjacency(B::SparseMatrixCSC, S::Vector{Int64}, Self::Bool=true)
    N = size(B,1)
    L = reduce(union, map(x->Set(GetAdjacency(B,x,false)), S))
    if Self
        L = union(L, Set(S))
    else
        L = setdiff(L, Set(S))
    end
    return L
end

# Returns whether the connected subgraph S is in a connected component at least of size Size.
# Basically finding the connected component (in an inefficient way) but can return quickly if Size is small.
function ConnectedComponentSizeAtLeast(B::SparseMatrixCSC, S::Vector{Int64}, Size::Integer)
    if length(S) >= Size
        return true
    end
    N = size(B,1)
    L = copy(S)
    size_old = length(L)
    while true
        for V in L
            L = union(L, GetAdjacency(B,V,true))
            if length(L) >= Size
                return true
            end
        end
        if size_old == length(L)
            return false
        else
            size_old = length(L)
        end
    end
end

function GetVolume(B::SparseMatrixCSC, S::Vector{Int64})
    sum(map(v->GetDegree(B,v), S))
end

function GetVolume(B::SparseMatrixCSC)
    N = size(B,1)
    sum(map(v->GetDegree(B,v), 1:N))
end

function GetAllDegrees(B::SparseMatrixCSC)
    N = size(B,1)
    collect(zip(1:N, map(v->GetDegree(B,v), 1:N)))
end

function GetAllDegreesFile(B::SparseMatrixCSC, FileName::String)
    io = open(FileName, "w")
    for i = 1:size(B,1)
        deg = GetDegree(B,i)
        wtf = string(deg, "\n")
        println(wtf)
        write(io, wtf)
    end
    close(io)
end

function GetInducedVolume(B::SparseMatrixCSC, S::Vector{Int64})
    N = size(S,1)
    sum(map(v->GetDegree(B[S,S],v), 1:N))
end

function PopSourceForFlowNetworkResult(S::Vector{Int64})
    S_ret = copy(S)
    popfirst!(S_ret)
    return map(x->x-1, S_ret)
end

function GetOrderByDegreeGraphIndices(B::SparseMatrixCSC)
    globalDegreeZip = collect(zip(1:size(B,1), map(x -> GetDegree(B,x), 1:size(B,1))))
    orderedIndices = sort(globalDegreeZip, by=x->x[2])
end

# Set an upper bound for node to be taken by degree.
# c_max: the maximum degree in C
# N: #nodes of the entire graph
# We cap the degree of nodes to be taken in R by c_max * X,
# Where X is smaller if c_max is close to N, and larger if c_max is small, with a bound (min_scale, max_scale).
# Otherwise, X = log(N / c_max) / log(log_scale).
mutable struct rNodeDegreeCap
    min_scale::Float64
    log_scale::Float64
    max_scale::Float64
end

DEFAULT_R_NODE_DEGREE_CAP = rNodeDegreeCap(1.0, 2.0, 10.0)
DEFAULT_R_NODE_DEGREE_CAP = rNodeDegreeCap(8.0, 2.0, 8.0)
NULL_R_NODE_DEGREE_CAP = rNodeDegreeCap(2.0^32, 2.0, 2.0^32)

function GetRNodeDegreeCap(c_max::Int64, N::Int64, RNodeDegreeCap::rNodeDegreeCap)
    scale = min(RNodeDegreeCap.max_scale, (max(RNodeDegreeCap.min_scale, log(RNodeDegreeCap.log_scale, N / c_max))))
    return floor(Int64, scale * c_max)
end