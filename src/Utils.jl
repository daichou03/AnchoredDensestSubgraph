# https://stackoverflow.com/a/27930367/10844976
function print_rgb(r, g, b, t)
    print("\e[1m\e[38;2;$r;$g;$b;249m",t)
end

# Assuming B is A's subset, find the array of indices of B's elements in A.
# NOTE: This function assumes A, B are SORTED!
# Thus performance is O(|B|).
function orderedSubsetIndices(A, B)
    indices = fill(-1, length(B))
    indA = 1
    indB = 1
    while indA <= length(A) && indB <= length(B)
        if A[indA] == B[indB]
            indices[indB] = indA
            indB += 1
        end
        indA += 1
    end
    return indices
end

# https://stackoverflow.com/questions/25678112/insert-item-into-a-sorted-list-with-julia-with-and-without-duplicates
insert_and_dedup!(v::Vector, x) = (splice!(v, searchsorted(v,x), [x]); v)

function emptyStringArray(len::Int64)
    ret = Array{String, 1}(undef, len)
    for i in 1:len
        ret[i] = ""
    end
    return ret
end
