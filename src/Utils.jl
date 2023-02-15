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


# Source: GPT-chat: "Write a function in Julia that have one thread run a function called func(), another thread waits for up to 300 seconds. If func() finishes in 300 seconds, return its result; otherwise, terminate func() and returns nothing."
# This function creates two tasks: task_func, which runs func(), and task_timeout, which sleeps for 300 seconds. The main thread waits for either of these tasks to complete using @taskwait.
# If task_func completes before task_timeout, its result is returned. Otherwise, task_func is terminated using schedule(task_func, TaskStatus.failed) and nothing is returned.
# Note that this implementation assumes that func() can be interrupted by terminating its corresponding task. If func() relies on external resources (such as network connections or file locks) that can't be released cleanly, it may not be safe to terminate it in this way.
function run_with_timeout(func, timeout)
    result = nothing
    task_func = @task func # func is passed as a callable, i.e., func(a).
    task_timeout = @task sleep(timeout)
    task_wait = @taskwait [task_func, task_timeout]

    if istaskdone(task_func)
        result = fetch(task_func)
    else
        schedule(task_func, TaskStatus.failed)
    end

    return result
end

function run_with_timeout(func, timeout)
    result = nothing
    task_func = @task func
    task_timeout = @task sleep(timeout)

    # Wait for either task to complete or timeout
    wait([task_func, task_timeout], timeout)

    task_done = istaskdone(task_func)
    if task_done
        result = task_done ? fetch(task_func) : nothing
    else
        schedule(task_func, TaskStatus.failed)
    end

    return (task_done, result)
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


