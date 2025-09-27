#!/usr/bin/env julia
# Test competitors module after reorganization

println("=== Testing Competitors Module ===")

# Test 1: Check competitors directory structure
println("1. Testing competitors directory structure...")
competitor_files = [
    "src/competitors/Competitors.jl",
    "src/competitors/home_coded/mrw.jl",
    "src/competitors/home_coded/greedyl.jl",
    "src/competitors/home_coded/flowseed.jl",
    "src/competitors/external/flowseed.jl"
]

all_competitor_files_exist = true
for file in competitor_files
    if isfile(file)
        println("   ✅ $file exists")
    else
        println("   ❌ $file missing")
        all_competitor_files_exist = false
    end
end

# Test 2: Test competitors module loading
println("\n2. Testing competitors module loading...")
try
    include("src/competitors/Competitors.jl")
    println("   ✅ Competitors module loads successfully")
    
    # Test if key functions are available
    if @isdefined(Competitors)
        println("   ✅ Competitors module is defined")
    end
catch e
    println("   ❌ Competitors module failed to load: $e")
end

# Test 3: Test main module with competitors
println("\n3. Testing main module with competitors...")
try
    include("src/AnchoredDensestSubgraph.jl")
    println("   ✅ Main module loads successfully with competitors")
    
    # Test if competitor functions are available through main module
    if @isdefined(MRW)
        println("   ✅ Competitor function MRW is available")
    else
        println("   ❌ Competitor function MRW not found")
    end
    
    if @isdefined(GreedyL)
        println("   ✅ Competitor function GreedyL is available")
    else
        println("   ❌ Competitor function GreedyL not found")
    end
    
    if @isdefined(FlowSeed)
        println("   ✅ Competitor function FlowSeed is available")
    else
        println("   ❌ Competitor function FlowSeed not found")
    end
catch e
    println("   ❌ Main module failed to load: $e")
end

# Test 4: Test that original CP/CX files are gone
println("\n4. Testing original CP/CX files are moved...")
original_files = [
    "src/CP_MRW.jl",
    "src/CP_GreedyL.jl",
    "src/CP_FlowSeed.jl",
    "src/CX_FlowSeed.jl"
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
println("\n=== Competitors Test Summary ===")
println("Competitor files created: $(all_competitor_files_exist ? "✅ PASS" : "❌ FAIL")")
println("Original CP/CX files moved: $(all_original_gone ? "✅ PASS" : "❌ FAIL")")

if all_competitor_files_exist && all_original_gone
    println("\n🎉 All competitors tests PASSED!")
else
    println("\n⚠️  Some competitors tests FAILED!")
end
