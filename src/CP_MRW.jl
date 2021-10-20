# Reference: this is a Julia implementation of MRW algorithm, :
# Bian, Yuchen, et al. "On multi-query local community detection." 2018 IEEE international conference on data mining (ICDM). IEEE, 2018.
# Single-query version only.

using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Utils.jl")


MRW_ALPHA = 0.1
MRW_BETA = 0.6
MRW_K = 5


# P: transition graph
# alpha, beta: decay parameters as per paper
# K: sliding window length
# tol: tolerance (of score difference only)
# verbose: print t etc.
function MRW(P, q, alpha=MRW_ALPHA, beta=MRW_BETA, K=MRW_K, tol=1e-6, verbose=false)
    # Initialize
    N = size(P, 1)
    t = 1
    x_0 = zeros(N)
    x_0[q] = 1.0
    x_old = copy(x_0)
    e = Any[]
    push!(e, fill(1/N, N))
    v_old = fill(1/N, N)
    converge = false
    # Theorem 2's cutoff
    cutoff_t = log(tol) / log(beta)
    while !converge && (t <= cutoff_t)
        # Update x
        x_new = alpha * (P * x_old) + (1 - alpha) * v_old
        # Add newest e
        keyWeight = 0.0
        keyIndices = []
        # TODO: the array does not look like sparse at all, performance has to be N related.
        for i = 1:N
            if almostEqual(x_new[i], keyWeight)
                push!(keyIndices, i)
            elseif keyWeight < x_new[i]
                keyWeight = x_new[i]
                keyIndices = [i]
            end
        end
        push!(e, zeros(N))
        for i in keyIndices
            e[t+1][i] = 1 / length(keyIndices)
        end
        # Calculate sum of sliding window
        e_window_sum = zeros(N)
        for k in 1:K
            if t + 2 - k < 1
                e_window_sum += x_0
            else
                e_window_sum += e[t + 2 - k]
            end
        end
        # Update v
        v_new = beta ^ (t - 1) * e_window_sum / K + (1 - beta ^ (t - 1)) * v_old
        # Check convergence
        converge = true
        for i = 1:N
            if !almostEqual(x_new[i], x_old[i], tol)
                converge = false
                break
            end
        end
        # Next loop
        t += 1
        x_old = x_new
        v_old = v_new
    end
    if verbose
        println("Number of iterations: " + (t - 1))
    end
    return x_old
end

function MRW_topK(P, q, topK, alpha=MRW_ALPHA, beta=MRW_BETA, K=MRW_K, tol=1e-6)
    x = MRW(P, q, alpha, beta, K, tol)
    return collect(partialsortperm(x, 1:topK, rev=true))
end

# partialsortperm(x, 1:k, rev=true) to get first k results.
