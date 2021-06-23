using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("maxflow.jl")
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Core_algorithm_yd.jl")
include("Test_utils_yd.jl")
include("Utils.jl")
include("CP_GreedyL.jl")
include("CS_Amazon.jl")

function GetDensity(B::SparseMatrixCSC, S::Vector{Int64})
    return GetVolume(B[S,S]) / length(S)
end

function GetAnchoredDensity(B::SparseMatrixCSC, R::Vector{Int64}, S::Vector{Int64})
    S_in_R = intersect(S, R)
    S_out_R = setdiff(S, S_in_R)
    return (GetVolume(B[S,S]) - GetVolume(B, S_out_R)) / length(S)
end

function GetConductance(B::SparseMatrixCSC, S::Vector{Int64})
    return 1 - GetVolume(B[S,S]) / GetVolume(B, S)
end

function GetLocalConductance(B::SparseMatrixCSC, R::Vector{Int64}, S::Vector{Int64}, epsilon::Float64=1.0)
    O_R = GetVolume(B, R) - epsilon * GetVolume(B, setdiff(S, R)) - GetVolume(B, setdiff(R, S))
    return O_R > 0 ? (GetVolume(B, S) - GetVolume(B[S,S])) / O_R : Inf
end

function GetLScore(B::SparseMatrixCSC, S::Vector{Int64})
    return LScore(B, S)
end

function GetPropRinS(R::Vector{Int64}, S::Vector{Int64})
    return 1 - length(setdiff(R, S)) / length(R)
end

function GetPropSoutR(R::Vector{Int64}, S::Vector{Int64})
    return length(setdiff(S, R)) / length(S)
end

function ReportCommunity(B::SparseMatrixCSC, R::Vector{Int64}, S::Vector{Int64})
    return join([length(S),
        GetDensity(B,S),
        GetAnchoredDensity(B,R,S),
        GetConductance(B,S),
        GetLocalConductance(B,R,S),
        GetLScore(B,S),
        GetPropRinS(R,S),
        GetPropSoutR(R,S)], "|")
end

function BulkReportCommunity(B::SparseMatrixCSC, Rs::Any, Ss::Any, Name::String)
    folder = string(CS_AMAZON_FOLDER, "EV-", Name, "/")
    mkpath(folder)
    for i in 1:length(Rs)
        io = open(string(folder,i,".txt"), "w")
        for j in 1:length(Rs[i])
            write(io, string(ReportCommunity(B, Rs[i][j], Ss[i][j]), "\n"))
        end
        close(io)
    end
end

# For Graph Editor

function ExportGraphEditorR(R, Ss, Name::String, Folder::String=string(CS_AMAZON_FOLDER, "GraphEditor/"))
    RUnion = copy(R)
    for s in Ss
        RUnion = union(RUnion, s)
    end
    sort!(RUnion)
    RUnionN = sort(GetComponentAdjacency(B, RUnion, true))

    # RIndices, edgelist, R, S1, S2, ...
    RUnionSubsetInds = orderedSubsetIndices(RUnionN, RUnion)
    RsubsetInds = orderedSubsetIndices(RUnionN, sort(R))
    SsubsetIndss = Any[]
    for i in 1:length(Ss)
        append!(SsubsetIndss, 0)
        SsubsetIndss[i] = orderedSubsetIndices(RUnionN, sort(Ss[i]))
    end
    Bsubset = B[RUnionN, RUnionN]

    mkpath(string(Folder,Name))
    io_edgelist = open(string(Folder,Name,"/edgelist.csv"), "w")
    for v1 in RUnionSubsetInds
        v1N = GetAdjacency(Bsubset, v1, false)
        for v2 in v1N
            if v1 < v2 || !(v2 in RUnionSubsetInds)
                write(io_edgelist, string(v1, ",", v2, "\n"))
            end
        end
    end
    close(io_edgelist)

    io_inds = open(string(Folder,Name,"/indices.csv"), "w")
    for v in RUnionN
        write(io_inds, string(v, "\n"))
    end
    close(io_inds)

    for i in 1:length(SsubsetIndss)
        io_inds = open(string(Folder,Name,"/S-", i,".csv"), "w")
        write(io_inds, string("Id,Label,Node,Color,Size", "\n"))
        for v in 1:length(RUnionN)
            color = string("#", ((v in SsubsetIndss[i]) ? lpad(string(204*256^(i-1),base=16), 6, "0") : "FFFFFF")) # Note only work if i <= 3
            size = ((v in RsubsetInds) ? 1 : 0)
            write(io_inds, string(join([v, v, v, color, size], ","), "\n"))
        end
        close(io_inds)
    end
end
