using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra

function GetDegree(B::SparseMatrixCSC, V::Int64)
    sum(B[V,:])
end

function GetAdjacency(B::SparseMatrixCSC, V::Int64)
    N = size(B,1)
    map(z->z[1], filter(a->a[2]>0, collect(zip(1:N,B[V,:]))))
end

function GetVolume(B::SparseMatrixCSC, S::Vector{Float64})
    sum(map(v->GetDegree(B,v), S))
end

function GetAllDegrees(B::SparseMatrixCSC)
    N = size(B,1)
    collect(zip(1:N, map(v->GetDegree(B,v), 1:N)))
end
