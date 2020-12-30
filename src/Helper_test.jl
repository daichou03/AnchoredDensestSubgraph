using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
include("maxflow.jl")
include("Helper_io.jl")
include("Helper_yd.jl")

mutable struct rSeed
    seed::Int64
    adj::Vector{Int64}
    degree::Float64
    volume::Float64
    induced_volume::Float64
    induced_maximum_density::Float64
    induced_maximum_densest_graph_size::Int64
    local_density::Float64
    local_densest_graph_size::Int64
end

function GetDegree(B::SparseMatrixCSC, V::Int64)
    sum(B[V,:])
end

function GetAdjacency(B::SparseMatrixCSC, V::Int64, Self::Bool)
    N = size(B,1)
    L = map(z->z[1], filter(a->a[2]>0, collect(zip(1:N,B[V,:]))))
    if Self
        L = prepend!(L, V)
    end
    return L
end

function GetAdjacency(B::SparseMatrixCSC, V::Int64)
    GetAdjacency(B, V, true)
end

function GetLeaveHighestDegAdjacency(B::SparseMatrixCSC, V::Int64)
    L = GetAdjacency(B, V, false)
    highest_index = findmax(map(z->GetDegree(B,z), L))[2]
    L[highest_index] = V # Remove the one with highest index and replace it with the seed
    return L
end

function GetVolume(B::SparseMatrixCSC, S::Vector{Int64})
    sum(map(v->GetDegree(B,v), S))
end

function GetAllDegrees(B::SparseMatrixCSC)
    N = size(B,1)
    collect(zip(1:N, map(v->GetDegree(B,v), 1:N)))
end

function GetInducedVolume(B::SparseMatrixCSC, S::Vector{Int64})
    sum(B[S,S])
end

function GetSeedReport(B::SparseMatrixCSC, V::Int64)
    adj = GetAdjacency(B, V)
    inducedMD = GlobalMaximumDensity(B[adj,adj])
    localMD = LocalMaximumDensity(B, adj)
    rSeed(V, adj, GetDegree(B, V), GetVolume(B, adj), GetInducedVolume(B, adj), inducedMD.alpha_star, length(inducedMD.source_nodes)-1, localMD.alpha_star, length(localMD.source_nodes)-1)
end

function GetLeaveHighestDegSeedReport(B::SparseMatrixCSC, V::Int64)
    adj = GetLeaveHighestDegAdjacency(B, V)
    inducedMD = GlobalMaximumDensity(B[adj,adj])
    localMD = LocalMaximumDensity(B, adj)
    rSeed(V, adj, GetDegree(B, V)-1, GetVolume(B, adj), GetInducedVolume(B, adj), inducedMD.alpha_star, length(inducedMD.source_nodes)-1, localMD.alpha_star, length(localMD.source_nodes)-1)
end

function SearchForNonDegeneratingSeed(B::SparseMatrixCSC)
    r = Vector{Int64}()
    for i = 1:size(B,1)
        rep = GetSeedReport(B,i)
        if rep.local_density - rep.induced_maximum_density > 1e-6
            println(string(i, " is a non-degenerate seed."))
            push!(r,i)
        end
    end
    r
end

function SearchForNonDegeneratingLeaveHighestDegSeed(B::SparseMatrixCSC, init::Int64)
    r = Vector{Int64}()
    for i = init:size(B,1)
        rep = GetLeaveHighestDegSeedReport(B,i)
        if rep.local_density - rep.induced_maximum_density > 1e-6
            println(string(i, " is a non-degenerate seed."))
            push!(r,i)
        end
    end
    r
end

function SearchForNonDegeneratingLeaveHighestDegSeed(B::SparseMatrixCSC)
    SearchForNonDegeneratingLeaveHighestDegSeed(B, 1)
end

function GetOneHopHighestDegreeNeighbourList(B::SparseMatrixCSC)
    N = size(B,1)
    List = repeat([0], N)
    for i = 1:N
        L = GetAdjacency(B, i, false)
        List[i] = L[findmax(map(z->GetDegree(B,z), L))[2]]
    end
    return List
end

function SearchWithOneHopHighestDegreeNeighbour(B::SparseMatrixCSC, NeighbourList::Vector{Int64})
    N = size(B,1)
    UniqueList = sort(unique(NeighbourList))
    Len = length(UniqueList)
    println(string("# Unique highest degree neighbours: ",Len))
    r = Vector{Int64}()
    for i = 1:Len
        V = UniqueList[i]
        rep = GetSeedReport(B,V)
        if rep.local_density - rep.induced_maximum_density > 1e-6
            println(string(V, " is a non-degenerate seed."))
            push!(r,V)
        end
    end
    r
end

fbgov = readIN("../Example/fbgov.in")
