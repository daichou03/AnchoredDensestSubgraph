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

# function GetPropSoutR(R::Vector{Int64}, S::Vector{Int64})
#     return length(setdiff(S, R)) / length(S)
# end

function GetPropSinR(R::Vector{Int64}, S::Vector{Int64})
    return 1 - length(setdiff(S, R)) / length(S)
end

function ReportCommunity(B::SparseMatrixCSC, R::Vector{Int64}, S::Vector{Int64})
    return join([length(S),
        GetDensity(B,S),
        GetAnchoredDensity(B,R,S),
        GetConductance(B,S),
        GetLocalConductance(B,R,S),
        GetLScore(B,S),
        GetPropRinS(R,S),
        GetPropSinR(R,S)], "|")
end

function BulkReportCommunity(B::SparseMatrixCSC, Rs::Any, Ss::Any, TestName::String, AlgName::String)
    folder = string(CS_AMAZON_FOLDER, "Report/", TestName, "/EV-", AlgName, "/")
    mkpath(folder)
    for i in 1:length(Rs)
        io = open(string(folder,i,".txt"), "w")
        for j in 1:length(Rs[i])
            write(io, string(ReportCommunity(B, Rs[i][j], Ss[i][j]), "\n"))
        end
        close(io)
    end
    BulkReportRSize(Rs, TestName)
end

function BulkReportRSize(Rs::Any, TestName::String)
    folder = string(CS_AMAZON_FOLDER, "Report/", TestName, "/R/")
    mkpath(folder)
    for i in 1:length(Rs)
        io = open(string(folder,i,".txt"), "w")
        for j in 1:length(Rs[i])
            write(io, string(length(Rs[i][j]), "\n"))
        end
        close(io)
    end
end

ALG_REPORT_NAMES = ["EV-LA", "EV-GL", "EV-FS"]
REPORT_METRICS = ["Length", "Density", "R-Subgraph Density", "Conductance", "Local Conductance", "L", "% R in S", "% S in R"]
REPORT_METRIC_FOLDER_NAME = ["length", "density", "rsdensity", "conductance", "lconductance", "lscore", "rins", "sinr"]

IMPUTE_VALUES = [0.0, 0.0, 0.0, 1.0, 999999.0, 0.0, 0.0, 0.0]

function ImputeNaNs(Values::Vector{Float64}, ImputeValues = IMPUTE_VALUES)
    ret = copy(Values)
    for i in 1:length(Values)
        if isnan(ret[i]) || isinf(ret[i])
            ret[i] = ImputeValues[i]
        end
        
    end
    return ret
end

function IntegrateReport(TestName::String, NumReports::Int64=41)
    statsAlgs = []
    for i_alg in 1:length(ALG_REPORT_NAMES)
        folder = string(CS_AMAZON_FOLDER, "Report/", TestName, "/", ALG_REPORT_NAMES[i_alg], "/")
        statsDegs = []
        for i_deg in 1:NumReports
            io = open(string(folder,i_deg,".txt"))
            stats = Array{Float64}(undef, length(REPORT_METRIC_FOLDER_NAME))
            count = 0
            while !eof(io)
                stat = ImputeNaNs(map(x->parse(Float64, x), split(readline(io), "|")))
                count += 1
                for j in 1:length(REPORT_METRIC_FOLDER_NAME)
                    stats[j] += stat[j]
                end
            end
            close(io)
            append!(statsDegs, 0)
            statsDegs[i_deg] = map(x->x/count, stats)
        end
        append!(statsAlgs, 0)
        statsAlgs[i_alg] = statsDegs
    end
    # R length
    folder = string(CS_AMAZON_FOLDER, "Report/", TestName, "/R/")
    rLengths = []
    for i_deg in 1:NumReports
        io = open(string(folder,i_deg,".txt"))
        rLength = 0
        count = 0
        while !eof(io)
            rLength += parse(Int64, readline(io))
            count += 1
        end
        close(io)
        append!(rLengths, 0)
        rLengths[i_deg] = map(x->x/count, rLength)
    end
    # For each metrics output data
    for i_metric in 1:length(REPORT_METRIC_FOLDER_NAME)
        output_folder = string(CS_AMAZON_FOLDER, "ReportIntegrated/", TestName, "/")
        mkpath(output_folder)
        io = open(string(output_folder, string(REPORT_METRIC_FOLDER_NAME[i_metric], ".txt")), "w")
        for i_deg in 1:NumReports
            line = []
            for i_alg in 1:length(ALG_REPORT_NAMES)
                append!(line, 0)
                line[i_alg] = statsAlgs[i_alg][i_deg][i_metric]
            end
            # Append length of R
            if i_metric == 1
                append!(line, 0)
                line[4] = rLengths[i_deg]
            end
            write(io, string(join(line, " "), "\n"))
        end
        close(io)
    end
end

# For Gephi
function ExportGraphEditor(R, Ss, Name::String, Folder::String=string(CS_AMAZON_FOLDER, "GraphEditor/"))
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
