using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase # TODO: To install
using Random
include("maxflow.jl")
include("Helper_io.jl")
include("Core_algorithm_yd.jl")
include("Graph_utils_yd.jl")
include("Utils.jl")

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

#--------------------------------------------------------
# Degeneracy test on R based on seed + neighbouring nodes
#--------------------------------------------------------

# Sampling by:
# chosen vertex and all its neighbours

function GetSeedAndNeighbours(B::SparseMatrixCSC, V::Int64)
    adj = GetAdjacency(B, V)
    GetGenericSeedReport(B,V,adj)
end

function SearchForNonDegeneratingSeedAndNeighbours(B::SparseMatrixCSC)
    r = Vector{Int64}()
    for i = 1:size(B,1)
        rep = GetSeedAndNeighbours(B,i)
        if rep.local_density - rep.induced_maximum_density > 1e-6
            println(string(i, " is a non-degenerate seed."))
            push!(r,i)
        end
    end
    r
end

# Test if the seed vertices with higher degree have higher chance to be non-degenerating

function GetOneHopHighestDegreeNeighbourList(B::SparseMatrixCSC)
    N = size(B,1)
    List = repeat([0], N)
    for i = 1:N
        L = GetAdjacency(B, i, false)
        List[i] = L[findmax(map(z->GetDegree(B,z), L))[2]]
    end
    return List
end

# Same as SearchForNonDegeneratingSeedAndNeighbours, but only seeds that are highest degree neighbour of some nodes
function SearchFOrNonDegeneratingSeedAndNeighboursOneHopHighestDegreeOnly(B::SparseMatrixCSC, NeighbourList::Vector{Int64})
    N = size(B,1)
    UniqueList = sort(unique(NeighbourList))
    Len = length(UniqueList)
    println(string("# Unique highest degree neighbours: ",Len))
    r = Vector{Int64}()
    for i = 1:Len
        V = UniqueList[i]
        rep = GetSeedAndNeighbours(B,V)
        if rep.local_density - rep.induced_maximum_density > 1e-6
            println(string(V, " is a non-degenerate seed."))
            push!(r,V)
        end
    end
    r
end

# Sampling by:
# chosen vertex and all its neighbours but the one with highest degree

function GetLeaveHighestDegAdjacency(B::SparseMatrixCSC, V::Int64)
    L = GetAdjacency(B, V, false)
    highest_index = findmax(map(z->GetDegree(B,z), L))[2]
    L[highest_index] = V # Remove the one with highest index and replace it with the seed
    return L
end

function GetLeaveHighestDegSeedReport(B::SparseMatrixCSC, V::Int64)
    adj = GetLeaveHighestDegAdjacency(B, V)
    GetGenericSeedReport(B,V,adj)
end

function SearchForNonDegeneratingLeaveHighestDegSeed(B::SparseMatrixCSC, init::Int64=1)
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

# Sampling by:
# chosen all its neighbours of a vertex (excluding itself).

function GetSeedExcludingSelfReport(B::SparseMatrixCSC, V::Int64)
    adj = GetAdjacency(B, V, false)
    GetGenericSeedReport(B,V,adj)
end

function SearchForNonDegeneratingSeedExcludingSelf(B::SparseMatrixCSC, init::Int64=1)
    r = Vector{Int64}()
    for i = init:size(B,1)
        rep = GetSeedExcludingSelfReport(B,i)
        if rep.local_density - rep.induced_maximum_density > 1e-6
            println(string(i, " is a non-degenerate seed."))
            push!(r,i)
        end
    end
    r
end

#---------------------------------------------
# Degeneracy test on R based on random walking
#---------------------------------------------

# Sampling by:
# starting with chosen/random vertex, sample a random connected component with fixed size

# Starting with a chosen vertex, include a random neighbour of a random vertex in the set to the set until the size of the set reaches Size.
function GetRandomWalkUntilSize(B::SparseMatrixCSC, V::Int64, Size::Int64)
    r = [V]
    while length(r) < Size
        next = 0
        while next == 0 || next in r
            seed = rand(r)
            next = rand(GetAdjacency(B, seed, false))
        end
        append!(r, next)
    end
    return r
end

function TestDegeneracyOnRandomWalkUntilSize(B::SparseMatrixCSC, Size::Int64, Tests::Int64, ShowSeed::Bool=false)
    N = size(B,1)
    nonDegCount = 0
    for i = 1:Tests
        seed = rand(1:N)
        sample = GetRandomWalkUntilSize(B,seed,Size)
        rep = GetGenericSeedReport(B,seed,sample)
        nonDeg = rep.local_density - rep.induced_maximum_density > 1e-6
        nonDegCount += (nonDeg ? 1 : 0)
        text = string("Test ", i, ": ", GetGenericSeedReport(B,seed,sample))
        if ShowSeed
            if nonDeg
                print_rgb(255,64,128,text)
                println()
            else
                print_rgb(255,255,255,text)
                println()
            end
        end
    end
    print_rgb(128,128,255,string("Non-degenerating R count: ", nonDegCount))
    println("")
    return nonDegCount
end

function TestDegeneracyOnRandomWalkUntilDifferentSize(B::SparseMatrixCSC, SizeFrom::Int64, SizeUntil::Int64, SizeInterval::Int64, Tests::Int64)
    size = SizeFrom
    while size <= SizeUntil
        print_rgb(255,255,128,string("Size = ", size, ": "))
        TestDegeneracyOnRandomWalkUntilSize(B,size,Tests,false)
        size += SizeInterval
    end
end

function TestDegeneracyOnRandomWalkUntilDifferentSize(B::SparseMatrixCSC, SizeFrom::Int64, SizeUntil::Int64, Tests::Int64)
    TestDegeneracyOnRandomWalkUntilDifferentSize(B, SizeFrom, SizeUntil, 1, Tests)
end

function TestDegeneracyOnRandomWalkUntilDifferentSize(B::SparseMatrixCSC, SizeUntil::Int64, Tests::Int64)
    TestDegeneracyOnRandomWalkUntilDifferentSize(B, 2, SizeUntil, Tests)
end

# Sampling by:
# starting with GetSampleUntilSize, then randomly removing some high degree nodes. That may lead to disjoint sets of nodes.
# Note that Size = final |R|, so if Size = 200 and Removes = 100, will get an R with |R| = 300 first then remove 100 nodes from it.
function GetRandomWalkUntilSizeThenRemoveHighDensity(B::SparseMatrixCSC, V::Int64, Size::Int64, Removes::Int64, DensityWeightFactor::Union{Int64,Float64})
    r = GetRandomWalkUntilSize(B, V, Size + Removes)
    weights = map(x->x[2]^DensityWeightFactor, GetAllDegrees(B[r,r]))
    removes = sample(1:length(r), Weights(weights), Removes, replace=false)
    deleteat!(r, sort(removes))
    return r
end

function BulkRandomWalkUntilSizeThenRemoveHighDensity(B::SparseMatrixCSC, Size::Int64, Removes::Int64, DensityWeightFactor::Union{Int64,Float64}, Tests::Int64)
    N = size(B, 1)
    samples = zeros(Int64, (Tests, Size))
    for row in eachrow(samples)
        ret = GetRandomWalkUntilSizeThenRemoveHighDensity(B,rand(1:N),Size,Removes,DensityWeightFactor)
        for i in 1:Size
            row[i] = ret[i]
        end
    end
    return samples
end

# for i = 1:size(samples, 1)
#     println(DetectConnectedComponents(B[samples[i,:], samples[i,:]]))
# end

function TestDegeneracyOnRandomWalkUntilSizeThenRemoveHighDensity(B::SparseMatrixCSC, Size::Int64, Removes::Int64, DensityWeightFactor::Union{Int64,Float64}, Tests::Int64, ShowSeed::Bool=false)
    N = size(B,1)
    nonDegCount = 0
    components = 0.0
    totalComponents = 0.0
    samples = BulkRandomWalkUntilSizeThenRemoveHighDensity(B,Size,Removes,DensityWeightFactor,Tests)
    for i = 1:Tests
        seed = rand(1:N)
        sample = samples[i,:]
        rep = GetGenericSeedReport(B,seed,sample)
        nonDeg = rep.local_density - rep.induced_maximum_density > 1e-6
        nonDegCount += (nonDeg ? 1 : 0)
        components = DetectConnectedComponents(B[sample,sample])
        totalComponents += components
        text = string("Test ", i, ": ", GetGenericSeedReport(B,seed,sample))
        if ShowSeed
            if nonDeg
                print_rgb(255,64,128,text)
                println()
            else
                print_rgb(255,255,255,text)
                println()
            end
            println(string("Number of components: ", components))
        end
    end
    print_rgb(128,128,255,string("Non-degenerating R count: ", nonDegCount))
    print_rgb(128,128,255,string(", average number of connected components: ", totalComponents / Tests))
    println("")
    return nonDegCount
end

function TestDegeneracyOnRandomWalkUntilSizeThenRemoveHighDensityDifferentRemoves(B::SparseMatrixCSC, Size::Int64, RemovesFrom::Int64, RemovesStep::Int64, RemovesTo::Int64, DensityWeightFactor::Union{Int64,Float64}, Tests::Int64)
    removes = RemovesFrom
    while removes <= RemovesTo
        print_rgb(255,255,128,string("Removes = ", removes, ": "))
        TestDegeneracyOnRandomWalkUntilSizeThenRemoveHighDensity(B,Size,removes,DensityWeightFactor,Tests,false)
        removes += RemovesStep
    end
end

# Sampling by:
# starting with chosen/random vertex, sample a random connected component until a minimum induced volume is reached.

function GetRandomWalkWithMinimumInducedVolume(B::SparseMatrixCSC, V::Int64, Volume::Int64)
    r = [V]
    while GetInducedVolume(B, r) < Volume
        next = 0
        while next == 0 || next in r
            seed = rand(r)
            next = rand(GetAdjacency(B, seed, false))
        end
        append!(r, next)
    end
    return r
end

function TestDegeneracyOnRandomWalkWithMinimumInducedVolume(B::SparseMatrixCSC, Volume::Int64, Tests::Int64, PrintRep::Bool=false)
    N = size(B,1)
    nonDegCount = 0
    for i = 1:Tests
        seed = rand(1:N)
        sample = GetRandomWalkWithMinimumInducedVolume(B,seed,Volume)
        rep = GetGenericSeedReport(B,seed,sample)
        nonDeg = rep.local_density - rep.induced_maximum_density > 1e-6
        nonDegCount += (nonDeg ? 1 : 0)
        text = string("Test ", i, ": ", GetGenericSeedReport(B,seed,sample))
        if PrintRep
            if nonDeg
                print_rgb(255,64,128,text)
                println()
            else
                print_rgb(255,255,255,text)
                println()
            end
        end
    end
    print_rgb(128,128,255,string("Non-degenerating R count: ", nonDegCount))
    println("")
    return nonDegCount
end

# Start with 2.
function TestDegeneracyOnRandomWalkWithDifferentMinimumInducedVolume(B::SparseMatrixCSC, VolumeUntil::Int64, Step::Int64, Tests::Int64)
    volume = 2
    while volume <= VolumeUntil
        print_rgb(255,255,128,string("InducedVolume >= ", volume, ": "))
        TestDegeneracyOnRandomWalkWithMinimumInducedVolume(B,volume,Tests,false)
        volume += Step
    end
end

# Sampling by:
# R with a fixed Size and must meet a minimum volume.
# Doing this by starting with GetSampleUntilSize (random walk), then using local search to change its components until its volume reaches (greater than or equal to) the required value.
# For now, this function returns an R with |R| = Size and induced_volume(R) >= MinVolume.
function GetFixedSizeWithMinimumDensity(B::SparseMatrixCSC, V::Int64, Size::Int64, MinVolume::Int64, MaxReseed::Int64=-1)
    if MinVolume > Size * (Size - 1)
        error(string("MinVolume is ", MinVolume, " which is greater than a graph with Size ", Size, " can possibly have."))
    end
    while MaxReseed != 0
        R = GetRandomWalkUntilSize(B, V, Size)
        retry = Size * 2 # TODO: number of retries reasonable? To see.
        while retry > 0 && GetInducedVolume(B, R) < MinVolume
            residualR = map(x -> Size - 1 - GetDegree(B[R,R], x), 1:Size) # For each vertex in r, the density it could possibly improve
            v = sample(1:Size, Weights(residualR))
            oldDegree = GetDegree(B[R,R], v)
            if oldDegree > 0
                soleNeighbour = rand(R[GetAdjacency(B[R,R], v, false)])
            else
                soleNeighbour = rand(R) # Any in R
            end
            Neighbours = shuffle(GetAdjacency(B, soleNeighbour))
            RMinus = copy(R)
            splice!(RMinus, v)
            for i = Neighbours
                newR = vcat(i, RMinus)
                newDegree = GetDegree(B[newR, newR], 1)
                if newDegree > oldDegree
                    retry += 1 # retry decreases only if can't find an improvement in a loop
                    R = newR
                    break
                end
            end
            retry -= 1
            # println(string("retry remaining: ", retry, ", current induced volume: ", GetInducedVolume(B, R)))
        end
        if GetInducedVolume(B, R) >= MinVolume
            return R
        end
        MaxReseed -= 1
    end
    error(string("Failed to find a subgraph with size ", Size, " and volume at least ", MinVolume, " in the given times of attempts. Either such subgraph does not exist or too rare for the purpose of sampling."))
end

function GetFixedSizeWithMinimumDensity(B::SparseMatrixCSC, Size::Int64, MinVolume::Int64)
    GetFixedSizeWithMinimumDensity(B, rand(1:size(B,1)), Size, MinVolume)
end

function TestDegeneracyOnFixedSizeWithMinimumDensity(B::SparseMatrixCSC, Size::Int64, MinVolume::Int64, Tests::Int64, ShowSeed::Bool=false)
    N = size(B,1)
    nonDegCount = 0
    for i = 1:Tests
        seed = rand(1:N)
        sample = GetFixedSizeWithMinimumDensity(B,seed,Size,MinVolume)
        rep = GetGenericSeedReport(B,seed,sample)
        nonDeg = rep.local_density - rep.induced_maximum_density > 1e-6
        nonDegCount += (nonDeg ? 1 : 0)
        text = string("Test ", i, ": ", GetGenericSeedReport(B,seed,sample))
        if ShowSeed
            if nonDeg
                print_rgb(255,64,128,text)
                println()
            else
                print_rgb(255,255,255,text)
                println()
            end
        end
    end
    print_rgb(128,128,255,string("Non-degenerating R count: ", nonDegCount))
    println("")
    return nonDegCount
end

function TestDegeneracyOnFixedSizeWithDifferentMinimumDensity(B::SparseMatrixCSC, Size::Int64, MinVolumeFrom::Int64, MinVolumeStep::Int64, MinVolumeTo::Int64, Tests::Int64, ShowSeed::Bool=false)
    volume = MinVolumeFrom
    while volume <= MinVolumeTo
        print_rgb(255,255,128,string("Min Volume >= ", volume, ": "))
        TestDegeneracyOnFixedSizeWithMinimumDensity(B,Size,volume,Tests)
        volume += MinVolumeStep
    end
end

# Relative density = 1 -> Clique only.
function TestDegeneracyOnDifferentSizeWithMinimumRelativeDensity(B::SparseMatrixCSC, SizeFrom::Int64, SizeUntil::Int64, SizeInterval::Int64, RelativeDensity::Float64, Tests::Int64)
    size = SizeFrom
    while size <= SizeUntil
        volume = Int64(floor(size * (size - 1) * RelativeDensity))
        print_rgb(255,255,128,string("Size = ", size, ": "))
        TestDegeneracyOnFixedSizeWithMinimumDensity(B,size,volume,Tests)
        size += SizeInterval
    end
end

# fbgov = readIN("../Example/fbgov.in")
eucore = RetrieveLargestConnectedComponent(readIN("../Example/email-Eu-core.in"))