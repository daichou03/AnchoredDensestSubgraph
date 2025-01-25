# This file explores properties of a result from a single case with the same R but different \w_12 and \w_24 (x and y).
using Plots
using Glob
include("Utils.jl")
include("Utils_io.jl")
include("LP_consts.jl")
include("Utils_graph.jl")
include("LP_evaluation.jl")


X_LOG_0_SMOOTH = 0.00001
Y_LOG_0_SMOOTH = 0.00001  # y-axis as log10. If y=0, to which value it is smoothed to for plotting.


## Extract x, y, z values
function extractLPResultFileData(files::Vector{String}, solverID::Int, n_range::UnitRange{Int}, extractZFunction::Function)
    # Initialize x, y, z arrays
    data = []

    for file in files
        # Extract x and y values from the filename
        matching = match(Regex(string("-", SOLVER_NAMES[solverID], "-([0-9.]+)-([0-9.]+)-")), file)
        if matching !== nothing
            x = parse(Float64, matching.captures[1])
            y = parse(Float64, matching.captures[2])

            # Read the entire file once
            rows = readlines(file)

            # Compute z values for all rows in n_range
            for n in n_range
                if n <= length(rows)
                    z = extractZFunction(rows, n)  # Pass pre-read rows to the function
                    push!(data, (x, y, z))
                end
            end
        end
    end

    # Sort the data by (x, y) in ascending order
    sorted_data = sort(data, by = x -> (x[1], x[2]))

    # Extract x_vals, y_vals, and z_vals from sorted data
    x_vals = [x[1] for x in sorted_data]
    y_vals = [x[2] for x in sorted_data]
    z_vals = [x[3] for x in sorted_data]

    return x_vals, y_vals, z_vals
end


#########################
# Extract single result #
#########################

function lpResultLengthTo2dCluster(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int})
    files = GetParameterizedLPResultFileNames(dataName, solverID, suffixName, RESULT_TYPE_SETS)

    # Function to extract z-value (length of integers in the n-th row)
    function rowLengthZ(rows::Vector{String}, n::Int)
        return length(split(rows[n], ","))
    end

    return extractLPResultFileData(files, solverID, n_range, rowLengthZ)
end


function lpResultDistinctTo2dCluster(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int})
    files = GetParameterizedLPResultFileNames(dataName, solverID, suffixName, RESULT_TYPE_SETS)
    unique_rows = Dict{String, Int}()

    # Function to compute z-value for a single n
    function distinctRowZ(rows::Vector{String}, n::Int)
        row_data = split(rows[n], ",")
        row_key = join(sort(row_data))  # Canonicalize row data
        if !haskey(unique_rows, row_key)
            unique_rows[row_key] = length(unique_rows) + 1
        end
        return unique_rows[row_key]
    end

    return extractLPResultFileData(files, solverID, n_range, distinctRowZ)
end


# Binary, 1 if result same as when x = 0, y = 0
function lpResultMatchTo2dCluster(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int})
    files = GetParameterizedLPResultFileNames(dataName, solverID, suffixName, RESULT_TYPE_SETS)

    # Extract the reference set S_ref (when x = 0, y = 0)
    ref_file_index = findfirst(f -> occursin("-LPLAS-0.0-0.0-", f), files)
    if ref_file_index === nothing
        error("No file found for x = 0 and y = 0")
    end
    ref_file = files[ref_file_index]
    rows_ref = readlines(ref_file)

    # Function to compute z-value for a single n
    function matchZ(rows::Vector{String}, n::Int)
        S_ref = parse.(Int, split(rows_ref[n], ","))
        S = parse.(Int, split(rows[n], ","))
        return S == S_ref ? 1 : 0
    end

    return extractLPResultFileData(files, solverID, n_range, matchZ)
end


function lpResultStatsTo2dCluster(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int}, columnName::String)
    files = GetParameterizedLPResultFileNames(dataName, solverID, suffixName, RESULT_TYPE_STATS)

    # Function to compute z-value for a single n
    function columnValueZ(file::String, n::Int)
        df = CSV.read(file, DataFrame)
        return df[n, columnName]  # Get the value from the specified column
    end

    return extractLPResultFileData(files, solverID, n_range, columnValueZ)
end


function lpResultDensityTo2dCluster(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int}, B::SparseMatrixCSC)
    files = GetParameterizedLPResultFileNames(dataName, solverID, suffixName, RESULT_TYPE_SETS)

    # Function to compute z-value for a single n
    function densityZ(rows::Vector{String}, n::Int)
        S = parse.(Int, split(rows[n], ","))
        submatrix = B[S, S]
        num_edges = nnz(submatrix)
        num_nodes = length(S)
        return 2 * num_edges / num_nodes
    end

    return extractLPResultFileData(files, solverID, n_range, densityZ)
end


function lpResultConductanceTo2dCluster(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int}, B::SparseMatrixCSC)
    files = GetParameterizedLPResultFileNames(dataName, solverID, suffixName, RESULT_TYPE_SETS)

    # Function to compute z-value for a single n
    function conductanceZ(rows::Vector{String}, n::Int)
        S = parse.(Int, split(rows[n], ","))
        volume_within = GetVolume(B[S, S])
        volume_total = GetVolume(B, S)
        return 1 - volume_within / volume_total
    end

    return extractLPResultFileData(files, solverID, n_range, conductanceZ)
end


# Number of nodes in S beyond 1-hop of R.
function lpResultLengthBeyond1HopTo2dCluster(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int}, B::SparseMatrixCSC)
    files = GetParameterizedLPResultFileNames(dataName, solverID, suffixName, RESULT_TYPE_SETS)
    anchor_sets = readAnchors(dataName, "Baseline")

    # Function to compute z-value for a single n
    function externalNodesZ(rows::Vector{String}, n::Int)
        R = anchor_sets[n]
        S = parse.(Int, split(rows[n], ","))
        R_neighbors = GetComponentAdjacency(B, R)
        external_nodes = setdiff(S, R_neighbors)
        return length(external_nodes)
    end

    return extractLPResultFileData(files, solverID, n_range, externalNodesZ)
end


###########################
# Extract aggregate reult #
###########################

# standardizeLPResultFunction as an interface for each different type of signature
function standardizeLPResultFunction(lpResultFunction::Function, args...; kwargs...)
    if haskey(kwargs, :B)
        # If `B` is present, call the function with `B`
        return lpResultFunction(args..., kwargs[:B])
    elseif haskey(kwargs, :columnName)
        # If `columnName` is present, call the function with `columnName`
        return lpResultFunction(args..., kwargs[:columnName])
    else
        # If no special keywords are present, call the function directly
        return lpResultFunction(args...)
    end
end


# Note that any parameters of lpResult... functions OTHER THAN dataName, suffixName and n (as n_range)
# Must be passed as kwargs, be it a kwarg in the original lpResult function or not.
# For example: lpResultAggregateTo2dCluster(lp)
function lpResultAggregateTo2dCluster(
    lpResultFunction::Function,
    dataName::String,
    solverID::Int64,
    suffixName::String,
    n_range::UnitRange{Int};
    kwargs...
)
    # Use the standardized interface to call the lpResultFunction
    x_vals, y_vals, z_vals = standardizeLPResultFunction(lpResultFunction, dataName, solverID, suffixName, n_range; kwargs...)

    # Aggregate z_vals for each (x, y) pair
    aggregated_data = Dict{Tuple{Float64, Float64}, Vector{Float64}}()
    for (x, y, z) in zip(x_vals, y_vals, z_vals)
        key = (x, y)
        if !haskey(aggregated_data, key)
            aggregated_data[key] = []
        end
        push!(aggregated_data[key], z)
    end

    # Compute average z_vals for each (x, y) pair
    x_vals, y_vals, z_vals = Float64[], Float64[], Float64[]
    for ((x, y), zs) in aggregated_data
        push!(x_vals, x)
        push!(y_vals, y)
        push!(z_vals, mean(zs))  # Compute average z-value
    end

    # Sort the data by (x, y)
    data_tuples = [(x, y, z) for (x, y, z) in zip(x_vals, y_vals, z_vals)]
    sorted_data = sort(data_tuples, by = t -> (t[1], t[2]))
    x_vals = [t[1] for t in sorted_data]
    y_vals = [t[2] for t in sorted_data]
    z_vals = [t[3] for t in sorted_data]
    return x_vals, y_vals, z_vals
end


#######################
# Plot x, y, z values #
#######################
function visualize_cluster((x_vals, y_vals, z_vals);
    log_x::Bool = false, log_y::Bool = true,
    smooth_x::Float64 = X_LOG_0_SMOOTH, smooth_y::Float64 = Y_LOG_0_SMOOTH,
    palette = :magma)
    # Smooth y-values: Replace 0 with the specified smooth_y value
    x_smoothed = [x == 0 && log_x ? smooth_x : x for x in x_vals]
    y_smoothed = [y == 0 && log_y ? smooth_y : y for y in y_vals]

    # Convert to grid for surface plotting
    x_unique = unique(x_smoothed)
    y_unique = unique(y_smoothed)
    z_matrix = [0.0 for _ in 1:length(x_unique), _ in 1:length(y_unique)]

    for i in 1:length(x_smoothed)
        x_idx = findfirst(==(x_smoothed[i]), x_unique)
        y_idx = findfirst(==(y_smoothed[i]), y_unique)
        z_matrix[x_idx, y_idx] = z_vals[i]
    end

    x_scaled = log_x ? log10.(x_smoothed) : x_smoothed
    y_scaled = log_y ? log10.(y_smoothed) : y_smoothed

    # Generate y-axis ticks that show the original y-values
    x_ticks = (x_scaled, [string(x) for x in x_vals])
    y_ticks = (y_scaled, [string(y) for y in y_vals])

    # Scatter plot without explicit markers, and annotate the z values
    plt = scatter(
        x_scaled,
        y_scaled,
        marker_z = z_vals,           # Use z-values for color coding
        title = "2D Scatter Plot of stats per x and y value",
        xlabel = "ω_{12}",
        ylabel = "ω_{24}",
        xticks = x_ticks,
        yticks = y_ticks,
        color = palette,            # Colormap
        legend = false,
        markersize = 12               # Remove marker circles
    )

    # Annotate each point with z-values, formatted to 4 significant digits
    for (x, y, z) in zip(x_scaled, y_scaled, z_vals)
        annotate!(x, y, text(z, :white, 10))
    end

    return plt
end


function visualize_contour((x_vals, y_vals, z_vals);
    log_x::Bool = false, log_y::Bool = true,
    smooth_x::Float64 = X_LOG_0_SMOOTH, smooth_y::Float64 = Y_LOG_0_SMOOTH,
    palette = :viridis, log_z::Bool = false)
    # Smooth y-values: Replace 0 with the specified smooth_y value
    x_smoothed = [x == 0 && log_x ? smooth_x : x for x in x_vals]
    y_smoothed = [y == 0 && log_y ? smooth_y : y for y in y_vals]

    # Apply log10 to z_vals if log_z is true
    z_transformed = log_z ? log10.(z_vals) : z_vals
    x_unique = unique(x_smoothed)
    y_unique = unique(y_smoothed)

    # Create a z-matrix for contour plotting
    z_matrix = [NaN for _ in 1:length(x_unique), _ in 1:length(y_unique)]
    for i in 1:length(x_smoothed)
        x_idx = findfirst(==(x_smoothed[i]), x_unique)
        y_idx = findfirst(==(y_smoothed[i]), y_unique)
        z_matrix[x_idx, y_idx] = z_transformed[i]
    end

    x_scaled = log_x ? log10.(x_smoothed) : x_smoothed
    y_scaled = log_y ? log10.(y_smoothed) : y_smoothed

    # Generate y-axis ticks that show the original y-values
    x_ticks = (x_scaled, [string(x) for x in x_vals])
    y_ticks = (y_scaled, [string(y) for y in y_vals])
    x_unique_scaled = unique(x_scaled)
    y_unique_scaled = unique(y_scaled)

    # Plot the contour
    contour(
        x_unique_scaled,
        y_unique_scaled,
        z_matrix',
        title = "Contour Plot",
        xlabel = "ω_{12}",
        ylabel = "ω_{24}",
        xticks = x_ticks,
        yticks = y_ticks,
        color = palette,  # Colormap for trends
        fill = true,       # Filled contours
        legend = :topright
    )
end

# TODO:
# Handle FNLA
# both log
# visualize_contour((lpResultAggregateTo2dCluster(lpResultStatsTo2dCluster, dataname, suffix, 1:100, columnName="ext_time")), palette=:viridis)