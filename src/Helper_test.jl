using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase # TODO: To install
using Random
include("maxflow.jl")
include("Helper_io.jl")
include("Helper_yd.jl")
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

function GetDegree(B::SparseMatrixCSC, V::Int64)
    sum(B[V,:])
end

function GetAdjacency(B::SparseMatrixCSC, V::Int64, Self::Bool=true)
    N = size(B,1)
    L = map(z->z[1], filter(a->a[2]>0, collect(zip(1:N,B[V,:]))))
    if Self
        L = prepend!(L, V)
    end
    return L
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

# Not very efficient on large?
# Returns the Set of the connected component that #1 vertex is in.
function ExtractConnectedComponent(B::SparseMatrixCSC)
    explored = Set(GetAdjacency(B,1,true))
    adj_set = setdiff(explored, Set([1]))
    while length(adj_set) > 0
        adj_set = setdiff(SetGetComponentAdjacency(B, collect(adj_set), false), explored)
        explored = union(explored, adj_set)
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

function RetrieveLargestConnectedComponent(B::SparseMatrixCSC)
    largestCCIndex = DetectConnectedComponents(B, true, false)
    remaining = copy(B)
    for i = 1:largestCCIndex-1
        nextCC = ExtractConnectedComponent(remaining)
        remaining_components = collect(setdiff(Set(1:size(remaining, 1)), nextCC))
        remaining = remaining[remaining_components, remaining_components]
    end
    nextCC = collect(ExtractConnectedComponent(remaining))
    return B[nextCC,nextCC]
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
function GetComponentAdjacency(B::SparseMatrixCSC, S::Vector{Int64}, Self::Bool=true)
    return collect(SetGetComponentAdjacency(B,S,Self))
end

function SetGetComponentAdjacency(B::SparseMatrixCSC, S::Vector{Int64}, Self::Bool=true)
    N = size(B,1)
    L = reduce(union, map(x->Set(GetAdjacency(B,x,true)), S))
    if !Self
        L = setdiff(L, Set(S))
    end
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

function GetGenericSeedReport(B::SparseMatrixCSC, V::Int64, R::Vector{Int64})
    inducedMD = GlobalMaximumDensity(B[R,R])
    localMD = LocalMaximumDensity(B, R)
    rSeed(V, R, GetDegree(B, V), GetVolume(B, R), GetInducedVolume(B, R), inducedMD.alpha_star, length(inducedMD.source_nodes)-1, localMD.alpha_star, length(localMD.source_nodes)-1)
end

# Sampling by:
# chosen vertex and all its neighbours

function GetSeedReport(B::SparseMatrixCSC, V::Int64)
    adj = GetAdjacency(B, V)
    GetGenericSeedReport(B,V,adj)
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

function SearchForNonDegeneratingLeaveHighestDegSeed(B::SparseMatrixCSC)
    SearchForNonDegeneratingLeaveHighestDegSeed(B, 1)
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

# Sampling by:
# starting with chosen/random vertex, sample a random connected component with fixed size

# Starting with a chosen vertex, include a random neighbour of a random vertex in the set to the set until the size of the set reaches Size.
function GetSampleUntilSize(B::SparseMatrixCSC, V::Int64, Size::Int64)
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

function RandomSampleUntilSize(B::SparseMatrixCSC, Size::Int64, Tests::Int64, ShowSeed::Bool=false)
    N = size(B,1)
    nonDegCount = 0
    for i = 1:Tests
        seed = rand(1:N)
        sample = GetSampleUntilSize(B,seed,Size)
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

# Start with 2.
function RandomSampleDifferentSize(B::SparseMatrixCSC, SizeUntil::Int64, Tests::Int64)
    RandomSampleDifferentSize(B, 2, SizeUntil, Tests)
end

function RandomSampleDifferentSize(B::SparseMatrixCSC, SizeFrom::Int64, SizeUntil::Int64, Tests::Int64)
    RandomSampleDifferentSize(B, SizeFrom, SizeUntil, 1, Tests)
end

function RandomSampleDifferentSize(B::SparseMatrixCSC, SizeFrom::Int64, SizeUntil::Int64, SizeInterval::Int64, Tests::Int64)
    size = SizeFrom
    while size <= SizeUntil
        print_rgb(255,255,128,string("Size = ", size, ": "))
        RandomSampleUntilSize(B,size,Tests,false)
        size += SizeInterval
    end
end

# Sampling by:
# starting with GetSampleUntilSize, then randomly removing some high degree nodes. That may lead to disjoint sets of nodes.
function GetSampleUntilSizeThenRemoveHighDensity(B::SparseMatrixCSC, V::Int64, Size::Int64, Removes::Int64, DensityWeightFactor::Union{Int64,Float64})
    r = GetSampleUntilSize(B, V, Size + Removes)
    weights = map(x->x[2]^DensityWeightFactor, GetAllDegrees(B[r,r]))
    removes = sample(1:length(r), Weights(weights), Removes, replace=false)
    deleteat!(r, sort(removes))
    return r
end

function BulkSampleUntilSizeThenRemoveHighDensity(B::SparseMatrixCSC, Size::Int64, Removes::Int64, DensityWeightFactor::Union{Int64,Float64}, Tests::Int64)
    N = size(B, 1)
    samples = zeros(Int64, (Tests, Size))
    for row in eachrow(samples)
        ret = GetSampleUntilSizeThenRemoveHighDensity(B,rand(1:N),Size,Removes,DensityWeightFactor)
        for i in 1:Size
            row[i] = ret[i]
        end
    end
    return samples
end

# for i = 1:size(samples, 1)
#     println(DetectConnectedComponents(B[samples[i,:], samples[i,:]]))
# end

function RandomSampleUntilSizeThenRemoveHighDensity(B::SparseMatrixCSC, Size::Int64, Removes::Int64, DensityWeightFactor::Union{Int64,Float64}, Tests::Int64, ShowSeed::Bool=false)
    N = size(B,1)
    nonDegCount = 0
    components = 0.0
    totalComponents = 0.0
    samples = BulkSampleUntilSizeThenRemoveHighDensity(B,Size,Removes,DensityWeightFactor,Tests)
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

function RandomSampleUntilSizeThenRemoveHighDensityDifferentRemoves(B::SparseMatrixCSC, Size::Int64, RemovesFrom::Int64, RemovesStep::Int64, RemovesTo::Int64, DensityWeightFactor::Union{Int64,Float64}, Tests::Int64)
    removes = RemovesFrom
    while removes <= RemovesTo
        print_rgb(255,255,128,string("Removes = ", removes, ": "))
        RandomSampleUntilSizeThenRemoveHighDensity(B,Size,removes,DensityWeightFactor,Tests,false)
        removes += RemovesStep
    end
end

# Sampling by:
# starting with chosen/random vertex, sample a random connected component with a minimum induced volume

function GetSampleUntilInducedVolume(B::SparseMatrixCSC, V::Int64, Volume::Int64)
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

function RandomSampleUntilInducedVolume(B::SparseMatrixCSC, Volume::Int64, Tests::Int64, PrintRep::Bool=false)
    N = size(B,1)
    nonDegCount = 0
    for i = 1:Tests
        seed = rand(1:N)
        sample = GetSampleUntilInducedVolume(B,seed,Volume)
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
function RandomSampleDifferentInducedVolume(B::SparseMatrixCSC, VolumeUntil::Int64, Step::Int64, Tests::Int64)
    volume = 2
    while volume <= VolumeUntil
        print_rgb(255,255,128,string("InducedVolume >= ", volume, ": "))
        RandomSampleUntilInducedVolume(B,volume,Tests,false)
        volume += Step
    end
end

# Sampling by:
# Starting with GetSampleUntilSize, use random walking to change its components until its volume reaches a specified value.
# For now, this function returns an R with |R| = Size and induced_volume(R) >= MinVolume.
function GetSampleUntilDensity(B::SparseMatrixCSC, V::Int64, Size::Int64, MinVolume::Int64, MaxReseed::Int64=-1)
    if MinVolume > Size * (Size - 1)
        error(string("MinVolume is ", MinVolume, " which is greater than a graph with Size ", Size, " can possibly have."))
    end
    while MaxReseed != 0
        R = GetSampleUntilSize(B, V, Size)
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

function GetSampleUntilDensity(B::SparseMatrixCSC, Size::Int64, MinVolume::Int64)
    GetSampleUntilDensity(B, rand(1:size(B,1)), Size, MinVolume)
end

function RandomSampleUntilDensity(B::SparseMatrixCSC, Size::Int64, MinVolume::Int64, Tests::Int64, ShowSeed::Bool=false)
    N = size(B,1)
    nonDegCount = 0
    for i = 1:Tests
        seed = rand(1:N)
        sample = GetSampleUntilDensity(B,seed,Size,MinVolume)
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

function RandomSampleUntilDifferentDensity(B::SparseMatrixCSC, Size::Int64, MinVolumeFrom::Int64, MinVolumeStep::Int64, MinVolumeTo::Int64, Tests::Int64, ShowSeed::Bool=false)
    volume = MinVolumeFrom
    while volume <= MinVolumeTo
        print_rgb(255,255,128,string("Min Volume >= ", volume, ": "))
        RandomSampleUntilDensity(B,Size,volume,Tests)
        volume += MinVolumeStep
    end
end

# Relative density = 1 -> Clique only.
function RandomSampleDifferentSizeRelativeDensity(B::SparseMatrixCSC, SizeFrom::Int64, SizeUntil::Int64, SizeInterval::Int64, RelativeDensity::Float64, Tests::Int64)
    size = SizeFrom
    while size <= SizeUntil
        volume = Int64(floor(size * (size - 1) * RelativeDensity))
        print_rgb(255,255,128,string("Size = ", size, ": "))
        RandomSampleUntilDensity(B,size,volume,Tests)
        size += SizeInterval
    end
end

fbgov = readIN("../Example/fbgov.in")
