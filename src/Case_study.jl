using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase # TODO: To install
using Random
using Base
include("maxflow.jl")
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Core_algorithm_yd.jl")
include("Test_utils_yd.jl")
include("Utils.jl")
include("Test_yd.jl")
include("Query_test_yd.jl")

DBLP_CI_FILE = "dblpciv4.in"

println("Reading DBLP citation data...")
B = readIN(DBLP_CI_FILE, "../CaseStudy/IN/")

V = 95485
# As an example, say V = 95485, for getting info from the raw citation data given this index (note -1 for the raw file):
# grep -n "#index95484" DBLP-citation-Jan8.txt
# Say you get row number = 6809987, then:
# Raw$ awk 'FNR>=6809987 && FNR<=6810020' DBLP-citation-Jan8.txt

C = GenerateUserInputSet(B,V,2,4)
R = GenerateReferenceSetFixedWalks(B,C)
a = StronglyLocalMaximumDensity(B,R)