using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("Helper_io.jl")
include("Utils.jl")
include("Graph_utils_yd.jl")

# This part of the code should be made safe to be included by both our algorithm and any external baseline algorithms that may not be compatiable with our algorithms.

CS_SIMPLE_FOLDER = "../CaseStudy/Simple/"

##################################
# For simple effectiveness tests #
##################################

SIMPLE_TEST_DATA_NAMES = ["amazon", "dblp", "youtube", "skitter", "livejournal", "orkut"]

# I/O R and V

function ExportSimpleRs(Vs::Vector{Int64}, Rs::Any, DataName::String)
    folder = folderString(CS_SIMPLE_FOLDER, DataName, "R")
    mkpath(folder)
    io_v = open(string(folder,"v.txt"), "w")
    io_r = open(string(folder,"r.txt"), "w")
    for i in 1:length(Vs)
        write(io_v, string(Vs[i], "\n"))
        write(io_r, string(join(Rs[i], ","), "\n"))
    end
    close(io_v)
    close(io_r)
end

function ImportSimpleRs(DataName::String)
    folder = folderString(CS_SIMPLE_FOLDER, DataName, "R")
    vs = Int64[]
    rs = Any[]

    io_v = open(string(folder,"v.txt"))
    while !eof(io_v)
        push!(vs, parse(Int64, readline(io_v)))
    end
    close(io_v)

    io_r = open(string(folder,"r.txt"))
    while !eof(io_r)
        push!(rs, map(x->parse(Int64, x), split(readline(io_r), ",")))
    end
        
    println(string(length(vs), " Rs imported."))
    return vs, rs
end

# I/O result set and time

function ExportSimpleResults(Ss::Any, Times::Any, DataName::String, AlgName::String)
    folder_s = folderString(CS_SIMPLE_FOLDER, DataName, "Result")
    folder_time = folderString(CS_SIMPLE_FOLDER, DataName, "Time")
    mkpath(folder_s)
    mkpath(folder_time)
    io_s = open(string(folder_s, AlgName, ".txt"), "w")
    io_time = open(string(folder_time, AlgName, ".txt"), "w")
    for j in 1:length(Ss)
        write(io_s, string(join(Ss[j], ","), "\n"))
        write(io_time, string(Times[j], "\n"))
    end
    close(io_s)
    close(io_time)
end

function ImportSimpleResults(DataName::String, AlgName::String)
    ss = Any[]
    io_s = open(string(folderString(CS_SIMPLE_FOLDER, DataName, "Result"), AlgName, ".txt"))
    while !eof(io_s)
        push!(ss, map(x->parse(Int64, x), split(readline(io_s), ",")))
    end
    close(io_s)

    times = Any[]
    io_time = open(string(folderString(CS_SIMPLE_FOLDER, DataName, "Time"), AlgName, ".txt"))
    while !eof(io_time)
        push!(times, parse(Float64, readline(io_time)))
    end
    close(io_time)

    println(string(length(ss), " results imported."))
    return ss, times
end
