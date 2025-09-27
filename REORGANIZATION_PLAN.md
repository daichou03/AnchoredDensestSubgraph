# Anchored Densest Subgraph Project Reorganization Plan

## Overview

This document outlines the reorganization plan for the Anchored Densest Subgraph research project, which implements algorithms from two academic papers in graph theory. The project was originally forked from "HypergraphFlowClustering" but has been completely transformed to focus on anchored densest subgraph algorithms.

## Project Context

- **Project Name**: Anchored Densest Subgraph
- **Original Fork**: HypergraphFlowClustering (competitor's code, now completely changed)
- **Implementation Language**: Julia
- **Research Focus**: Local community search algorithms

## Algorithm Classification

The project implements algorithms from two distinct papers:

### Paper 1: ADS (Anchored Densest Subgraph)
- **Core File**: `Core_algorithm_yd.jl` (now `src/ADS/algorithms.jl`)
- **Algorithm Type**: Flow network-based algorithms
- **Key Features**: 
  - Global algorithms (evaluate using entire graph)
  - Local algorithms (start with smaller working graph, expand as needed)
  - Better scalability with local approach

### Paper 2: GADS (Generalized Anchored Densest Subgraph)
- **Core File**: `LP_algorithm.jl` (now `src/GADS/algorithms.jl`)
- **Algorithm Type**: Linear programming-based algorithms
- **Key Features**: All algorithms are "local" in the aforementioned sense

## File Prefix Legend

Understanding the existing file naming conventions:

- **CP**: Competitors, home-coded
- **CS**: Universal case study materials (used by both papers)
- **CS2**: Newer universal case study materials
- **CX**: Competitors, external code (from original forked project)
- **LP**: Related to GADS only (Paper 2)
- **Utils**: Universal utilities

## Reorganization Strategy

### Phase 1: Directory Structure and Package Setup ✅ COMPLETED
- Create new directory structure
- Set up Julia package management with `Project.toml`
- Create placeholder directories for organized modules

### Phase 2: Source Code Reorganization ✅ IN PROGRESS

#### Commit 1: Create Julia Package Structure ✅ COMPLETED
- Created `Project.toml` with proper dependencies
- Created main module `src/AnchoredDensestSubgraph.jl`
- Set up basic package structure

#### Commit 2: Move ADS Algorithms ✅ COMPLETED
- Moved `Core_algorithm_yd.jl` → `src/ADS/algorithms.jl`
- Moved `maxflow.jl` → `src/ADS/flow_network.jl`
- Moved `Query_test_yd.jl` → `src/ADS/experiments/query_tests.jl`
- Moved `Test_degeneracy_yd.jl` → `src/ADS/experiments/degeneracy_tests.jl`
- Created `src/ADS/ADS.jl` module
- Updated include paths and module dependencies

#### Commit 3: Move GADS Algorithms ✅ COMPLETED
- Moved `LP_algorithm.jl` → `src/GADS/algorithms.jl`
- Moved `LP_consts.jl` → `src/GADS/lp_consts.jl`
- Moved `LP_load_solver.jl` → `src/GADS/lp_load_solver.jl`
- Moved `LP_evaluation.jl` → `src/GADS/lp_evaluation.jl`
- Moved `LP_compare_test.jl` → `src/GADS/experiments/lp_compare_test.jl`
- Created `src/GADS/GADS.jl` module
- Updated include paths and module dependencies

#### Commit 4: Fix Module Dependencies ✅ COMPLETED
- Fixed include paths in GADS and ADS modules
- Ensured proper module dependencies between ADS and GADS
- Resolved circular dependencies and path issues

#### Commit 5: Move Case Studies ✅ COMPLETED
- Moved all `CS_*` and `CS2_*` files to `src/case_studies/`
- Created `CaseStudies.jl` module for universal case study functionality
- Updated include paths in all case study files
- Fixed data file paths for new directory structure

### Phase 3: Remaining Reorganization Tasks

#### Commit 6: Move Competitor Algorithms (PENDING)
- Move `CP_*` files to `src/competitors/home_coded/`
- Move `CX_*` files to `src/competitors/external/`
- Create competitor modules and update references

#### Commit 7: Move Utilities (PENDING)
- Move `Utils_*` files to `src/utils/`
- Create utilities module
- Update all references to utility functions

#### Commit 8: Move Experimental Code (PENDING)
- Move experiment-specific code to appropriate modules
- Separate ADS and GADS experiments
- Organize performance evaluation code

#### Commit 9: Data and Results Reorganization (PENDING)
- Move data files to `data/` directory structure
- Organize results and performance reports
- Update data file paths throughout codebase

#### Commit 10: Documentation and Testing (PENDING)
- Create comprehensive documentation
- Add unit tests
- Create usage examples and tutorials

## Current Project Structure

```
src/
├── ADS/                           # Paper 1 algorithms (flow-based)
│   ├── algorithms.jl             # Core ADS algorithms
│   ├── flow_network.jl           # Flow network implementation
│   ├── experiments/              # ADS-specific experiments
│   │   ├── query_tests.jl
│   │   └── degeneracy_tests.jl
│   └── ADS.jl                    # ADS module
├── GADS/                         # Paper 2 algorithms (LP-based)
│   ├── algorithms.jl             # Core GADS algorithms
│   ├── lp_consts.jl              # LP constants
│   ├── lp_load_solver.jl         # LP solver loading
│   ├── lp_evaluation.jl          # LP evaluation
│   ├── experiments/              # GADS-specific experiments
│   │   └── lp_compare_test.jl
│   └── GADS.jl                   # GADS module
├── case_studies/                 # Universal case study implementations
│   ├── simple.jl                 # Simple case studies
│   ├── dblp.jl                   # DBLP case studies
│   ├── amazon.jl                 # Amazon case studies
│   ├── generic.jl                # Generic case study utilities
│   ├── plot.jl                   # Plotting utilities
│   └── CaseStudies.jl            # Case studies module
├── competitors/                  # Competitor algorithms (to be organized)
│   ├── home_coded/               # CP_* files
│   └── external/                 # CX_* files
├── utils/                        # Shared utilities (to be organized)
└── AnchoredDensestSubgraph.jl    # Main module
```

## Key Design Principles

1. **Paper-Based Organization**: Algorithms are organized by research paper (ADS vs GADS) rather than technical approach
2. **Universal vs Specific**: Case studies are universal, experiments are paper-specific
3. **Modular Design**: Each major component has its own module with clear interfaces
4. **Incremental Changes**: Each commit is testable and maintains functionality
5. **Clear Separation**: Distinct separation between algorithms, utilities, and experiments

## Testing Strategy

- Each commit is tested to ensure functionality is preserved
- Module loading tests verify proper dependencies
- Algorithm functionality tests ensure core features work
- Integration tests verify modules work together

## Git Branch Strategy

- **Branch**: `cursor-refactor`
- **Approach**: Incremental commits with testing after each major change
- **Safety**: Original code preserved in main branch

## Progress Summary

- ✅ **Commits 1-5**: Core algorithm and case study reorganization completed
- 🔄 **Commits 6-10**: Remaining reorganization tasks pending
- 📊 **Overall Progress**: ~50% complete

## Next Steps

1. Complete competitor algorithm reorganization (Commit 6)
2. Organize utilities and shared code (Commit 7)
3. Separate and organize experimental code (Commit 8)
4. Reorganize data and results (Commit 9)
5. Add documentation and testing (Commit 10)

## Benefits of Reorganization

1. **Clear Structure**: Easy to understand project organization
2. **Maintainability**: Modular design makes code easier to maintain
3. **Extensibility**: Clear separation allows for easy addition of new algorithms
4. **Documentation**: Better organization facilitates documentation
5. **Collaboration**: Clear structure makes it easier for others to contribute
6. **Testing**: Modular design enables better testing strategies

---

*Last Updated: January 12, 2025*
*Branch: cursor-refactor*
*Status: In Progress (Commits 1-5 Complete)*
