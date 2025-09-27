#!/usr/bin/env julia
# Test case studies module after reorganization

println("=== Testing Case Studies Module ===")

# Test 1: Check case studies directory structure
println("1. Testing case studies directory structure...")
case_study_files = [
    "src/case_studies/CaseStudies.jl",
    "src/case_studies/simple.jl",
    "src/case_studies/dblp.jl",
    "src/case_studies/amazon.jl",
    "src/case_studies/generic.jl",
    "src/case_studies/plot.jl"
]

all_case_study_files_exist = true
for file in case_study_files
    if isfile(file)
        println("   ✅ $file exists")
    else
        println("   ❌ $file missing")
        all_case_study_files_exist = false
    end
end

# Test 2: Test case studies module loading
println("\n2. Testing case studies module loading...")
try
    include("src/case_studies/CaseStudies.jl")
    println("   ✅ Case studies module loads successfully")
    
    # Test if key functions are available
    if @isdefined(CaseStudies)
        println("   ✅ Case studies module is defined")
    end
catch e
    println("   ❌ Case studies module failed to load: $e")
end

# Test 3: Test main module with case studies
println("\n3. Testing main module with case studies...")
try
    include("src/AnchoredDensestSubgraph.jl")
    println("   ✅ Main module loads successfully with case studies")
    
    # Test if case study functions are available through main module
    if @isdefined(ExportSimpleRs)
        println("   ✅ Case study function ExportSimpleRs is available")
    else
        println("   ❌ Case study function ExportSimpleRs not found")
    end
catch e
    println("   ❌ Main module failed to load: $e")
end

# Test 4: Test that original CS files are gone
println("\n4. Testing original CS files are moved...")
original_files = [
    "src/CS_Simple.jl",
    "src/CS_DBLP.jl",
    "src/CS_Amazon.jl",
    "src/CS_generic.jl",
    "src/CS2_plot.jl"
]

all_original_gone = true
for file in original_files
    if !isfile(file)
        println("   ✅ $file successfully moved")
    else
        println("   ❌ $file still exists (should be moved)")
        all_original_gone = false
    end
end

# Summary
println("\n=== Case Studies Test Summary ===")
println("Case study files created: $(all_case_study_files_exist ? "✅ PASS" : "❌ FAIL")")
println("Original CS files moved: $(all_original_gone ? "✅ PASS" : "❌ FAIL")")

if all_case_study_files_exist && all_original_gone
    println("\n🎉 All case studies tests PASSED!")
else
    println("\n⚠️  Some case studies tests FAILED!")
end

