using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("maxflow.jl")
include("Utils_io.jl")
include("Utils_graph.jl")
include("Utils_warmup.jl")
include("Utils.jl")

#############
# For Gephi #
#############
# Amazon: Folder=folderString(CS_AMAZON_FOLDER, "Single", "GraphEditor")
# DBLP: Folder=folderString(CS_DBLP_FOLDER, "Single", "GraphEditor")

function ExportGraphEditor(B::SparseMatrixCSC, R, Ss, Name::String, Folder::String)
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
        push!(SsubsetIndss, orderedSubsetIndices(RUnionN, sort(Ss[i])))
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
            color = string("#", ((v in SsubsetIndss[i]) ? lpad(string(204*256^(3-i),base=16), 6, "0") : "FFFFFF")) # Note only work if i <= 3
            size = ((v in RsubsetInds) ? 1 : 0)
            write(io_inds, string(join([v, v, v, color, size], ","), "\n"))
        end
        close(io_inds)
    end
end

# Changes:
# Output all 3 colours simultaneously for filtering
# Does not take neighbour nodes anymore. Instead, having degree and residual degree (from untaken neighbours) as labels.

function ExportGraphEditorForDBLP(B::SparseMatrixCSC, V::Int64, R::Vector{Int64}, Ss, Name::String, Folder::String)
    RUnion = copy(R)
    for s in Ss
        RUnion = union(RUnion, s)
    end
    sort!(RUnion)

    # RIndices, edgelist, R, S1, S2, ...
    RsubsetInds = orderedSubsetIndices(RUnion, sort(R))
    SsubsetIndss = Any[]
    for i in 1:length(Ss)
        push!(SsubsetIndss, orderedSubsetIndices(RUnion, sort(Ss[i])))
    end
    Bsubset = B[RUnion, RUnion]
    folderName = folderString(Folder,Name)

    mkpath(folderName)
    io_edgelist = open(string(folderName,"edgelist.csv"), "w")
    for v1 in RsubsetInds
        v1N = GetAdjacency(Bsubset, v1, false)
        for v2 in v1N
            if v1 < v2
                write(io_edgelist, string(v1, ",", v2, "\n"))
            end
        end
    end
    close(io_edgelist)

    io_inds = open(string(folderName,"indices.csv"), "w")
    for v in RUnion
        write(io_inds, string(v, "\n"))
    end
    close(io_inds)

    io_inds = open(string(folderName,"S.csv"), "w")
    write(io_inds, string("Id,Label,Node,Color,Polygon,Deg,IsKey,LogDeg,ResDeg,OriginInd,OriginName", "\n"))
    for v in 1:length(RUnion)
        color = "#"
        for i in 1:length(SsubsetIndss)
            color = string(color, v in SsubsetIndss[i] ? "FF" : "00")
        end
        polygon = ((v in RsubsetInds) ? 0 : 3)
        deg = GetDegree(B, RUnion[v])
        iskey = (v == orderedSubsetIndices(RUnion, [V])[1]) ? 2 : ((v in RsubsetInds) ? 1 : 0)
        logDeg = log(deg)
        resDeg = deg - GetDegree(Bsubset, v)
        originInd = RUnion[v]
        originName = replace(allNames[originInd], "_"=>" ")
        label = string(deg, " - ", originName)
        write(io_inds, string(join([v, label, v, color, polygon, deg, iskey, logDeg, resDeg, originInd, originName], ","), "\n"))
    end
    close(io_inds)
end