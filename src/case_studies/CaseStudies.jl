"""
CaseStudies.jl

Module for case study implementations that are universal to both ADS and GADS papers.
This module contains case study code for various datasets including Amazon, DBLP, 
and other real-world networks.

# Usage
```julia
using AnchoredDensestSubgraph.CaseStudies

# Load case study data
data = load_amazon_case_study()
data = load_dblp_case_study()

# Run case study evaluations
results = evaluate_case_study(data)
```
"""
module CaseStudies

# Core dependencies
using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
using CSV
using DataFrames

# Include utility functions
include("../Utils_io.jl")
include("../Utils_graph.jl")
include("../Utils.jl")

# Include case study implementations
include("amazon.jl")
include("amazon_la.jl")
include("amazon_stratified.jl")
include("dblp.jl")
include("dblp_la.jl")
include("evaluation_amazon.jl")
include("evaluation_simple.jl")
include("evaluation_single.jl")
include("simple.jl")
include("simple_la.jl")
include("generic.jl")
include("generic_la.jl")

# Include newer case study implementations
include("dblp_v2.jl")
include("dblp_generate.jl")
include("generic_v2.jl")
include("plot.jl")

# Export case study functions
export ExportSimpleRs, ImportSimpleRs, ProcessCaseStudy
export ProcessCaseStudyAmazon, ProcessCaseStudyDBLP, ProcessCaseStudySimple
export GenerateCaseStudyData, EvaluateCaseStudyResults

# Export evaluation functions
export ProcessEvaluationAmazon, ProcessEvaluationSimple, ProcessEvaluationSingle

# Export plotting functions
export PlotCaseStudyResults, GenerateCaseStudyPlots

end # module

