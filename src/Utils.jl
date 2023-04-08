using Dates

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

# A robust function that concatenate directory strings.
# The output ends with "/" regardless, and is ready to attach to a file name.
function folderString(x::String, y::String...)
    a = x
    if string(last(a)) != "/"
        a = string(a, "/")
    end
    if length(y) == 0
        return a
    end

    b = y[1]
    if length(b) > 0 && string(b[1]) == "/"
        b = b[2:end]
    end
    return folderString(string(a, b), y[2:end]...)
end

# Timer
TIMER_DATETIME = now()

function TimerReset()
    global TIMER_DATETIME
    TIMER_DATETIME = now()
end

function TimerLapValue()
    global TIMER_DATETIME
    now_0 = now()
    ret = (now_0-TIMER_DATETIME).value
    TIMER_DATETIME = now_0
    return ret / 1000
end


# 20230220: The tombstone of attempting to write a run-with-timeout function - doesn't work.
# returns (if_func_finished_before_timeout, result_of_func_if_finished)
# func needs to have no parameter.
# EMPTY_LOOP_EPSILON = 1e-4 # 20230215: It turned out that you just need to "do something" in the while loop, otherwise it behaves as if being recognized as "unconditional blocking" by the compiler.
# WAIT_PRECISION = 1e-4 # Whenever the current elapsed time * WAIT_PRECISION > current wait interval, doubles current wait interval.

# Rewrite after works:
# if (time() - t0) * WAIT_PRECISION > wait_interval
#     wait_interval *= 2
# end

# function run_with_timeout(func, timeout)


# Only up to second
function TimeAsName()
    return split(replace(string(now()), ":"=>"-"), ".")[1]
end

# Math
ALMOST_EQUAL_TOL = 1e-6
function almostEqual(a, b, tol=ALMOST_EQUAL_TOL)
    return abs(a-b) < tol
end


function f1score(p::Float64, r::Float64)
    if p == 0 && r == 0
        return 0.0
    end
    return 2 * p * r / (p + r)
end

function f1score(set1::Union{Set, Vector}, set2::Union{Set, Vector})
    if length(set1) == 0 || length(set2) == 0
        return 0.0
    end
    tp = length(intersect(set1, set2))
    fp = length(set1) - tp
    fn = length(set2) - tp
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    return f1score(precision, recall)
end
