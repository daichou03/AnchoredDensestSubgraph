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

# Concatenate folder strings
#function folderString()

function folderString(x::String, y::String...)
    a = x
    if string(last(a)) != "/"
        a = string(a, "/")
    end
    if length(y) == 0
        return a
    end

    b = y[1]
    if string(b[1]) == "/"
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


# returns (if_func_finished_before_timeout, result_of_func_if_finished)
# func needs to have no parameter.
EMPTY_LOOP_EPSILON = 1e-16 # 20230215: It turned out that you just need to "do something" in the while loop, otherwise it behaves as if being recognized as "unconditional blocking" by the compiler.

function run_with_timeout(func, timeout)
    task_func = @async begin
        func()
    end
    
    task_timeout = @async begin
        sleep(timeout)
    end
    
    # Wait for either task to complete.

    completed_task = nothing
    task_completed = false
    while true
        for task in [task_func, task_timeout]
            if istaskdone(task)
                completed_task = task
                task_completed = true
                break
            end
        end
        if task_completed
            break
        end
        sleep(EMPTY_LOOP_EPSILON)
        # print(string(istaskdone(task_func), " ", istaskdone(task_timeout), " ", task_completed)) # This is also weird: task_completed is false by the first time it should have been true.
    end
    
    # Determine which task completed and take appropriate action
    func_complete = completed_task == task_func
    if func_complete
        return (true, fetch(task_func))
    else
        return (false, nothing)
    end
end

# Only up to second
function TimeAsName()
    return split(replace(string(now()), ":"=>"-"), ".")[1]
end

# Math
ALMOST_EQUAL_TOL = 1e-6
function almostEqual(a, b, tol=ALMOST_EQUAL_TOL)
    return abs(a-b) < tol
end


