# This file explores properties of a result from a single case with the same R but different \w_12 and \w_24 (x and y).
using Plots
using Glob
using LaTeXStrings
include("Utils.jl")
include("Utils_io.jl")
include("LP_consts.jl")
include("Utils_graph.jl")
include("LP_evaluation.jl")


X_LOG_0_SMOOTH = 0.001
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

function lpResultLength(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int})
    files = GetParameterizedLPResultFileNames(dataName, solverID, suffixName, RESULT_TYPE_SETS)

    # Function to extract z-value (length of integers in the n-th row)
    function rowLengthZ(rows::Vector{String}, n::Int)
        return length(split(rows[n], ","))
    end

    return extractLPResultFileData(files, solverID, n_range, rowLengthZ)
end


function lpResultDistinct(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int})
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
function lpResultMatch(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int})
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


function lpResultStats(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int}, columnName::String)
    files = GetParameterizedLPResultFileNames(dataName, solverID, suffixName, RESULT_TYPE_STATS)

    # Function to compute z-value for a single n
    function columnValueZ(file::String, n::Int)
        df = CSV.read(file, DataFrame)
        return df[n, columnName]  # Get the value from the specified column
    end

    return extractLPResultFileData(files, solverID, n_range, columnValueZ)
end


function lpResultDensity(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int}, B::SparseMatrixCSC)
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


function lpResultConductance(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int}, B::SparseMatrixCSC)
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
function lpResultLengthBeyond1Hop(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int}, B::SparseMatrixCSC)
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


function lpResultF1Score(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int})
    files = GetParameterizedLPResultFileNames(dataName, solverID, suffixName, RESULT_TYPE_SETS)
    anchor_sets = readAnchors(dataName, "Baseline")

    # Function to compute F1-score for a single n
    function computeF1Score(rows::Vector{String}, n::Int)
        R = anchor_sets[n]  # Anchor set (ground truth)
        S = parse.(Int, split(rows[n], ","))  # Extract result set

        # Compute precision, recall, and F1-score
        true_positives = length(intersect(S, R))
        precision = true_positives / length(S)  # TP / (TP + FP)
        recall = true_positives / length(R)  # TP / (TP + FN)

        # Handle cases where precision or recall are zero
        if precision + recall == 0
            return 0.0  # F1-score is 0 when both are 0
        end

        f1_score = 2 * (precision * recall) / (precision + recall)
        return f1_score
    end

    return extractLPResultFileData(files, solverID, n_range, computeF1Score)
end


function lpResultJaccardSimilarity(dataName::String, solverID::Int64, suffixName::String, n_range::UnitRange{Int})
    files = GetParameterizedLPResultFileNames(dataName, solverID, suffixName, RESULT_TYPE_SETS)
    anchor_sets = readAnchors(dataName, "Baseline")

    # Function to compute Jaccard similarity for a single n
    function computeJaccard(rows::Vector{String}, n::Int)
        R = anchor_sets[n]  # Anchor set (ground truth)
        S = parse.(Int, split(rows[n], ","))  # Extract result set

        # Compute Jaccard similarity: |S ∩ R| / |S ∪ R|
        intersection_size = length(intersect(S, R))
        union_size = length(union(S, R))

        # Handle empty union case
        if union_size == 0
            return 0.0
        end

        jaccard_index = intersection_size / union_size
        return jaccard_index
    end

    return extractLPResultFileData(files, solverID, n_range, computeJaccard)
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
# For example: lpResultAggregate(lp)
function lpResultAggregate(
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
# smooth_value is before log, so likely you want 0.0001 rather than -4.
function smooth_and_scale(values::Vector{Float64}, log_scale::Bool, smooth_value::Float64)
    # Replace zeros if log scaling is enabled
    smoothed_values = [v == 0 && log_scale ? smooth_value : v for v in values]
    
    # Apply log transformation if required
    return log_scale ? log10.(smoothed_values) : smoothed_values
end


function visualize_cluster((x_vals, y_vals, z_vals);
    log_x::Bool = false, log_y::Bool = true,
    smooth_x::Float64 = X_LOG_0_SMOOTH, smooth_y::Float64 = Y_LOG_0_SMOOTH,
    palette = :magma, 
    xlabel = L"ω_{12}", ylabel = L"-ω_{24}")
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
        xlabel = xlabel,
        ylabel = ylabel,
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
    palette = :viridis, log_z::Bool = false,
    xlabel = L"ω_{12}", ylabel = L"-ω_{24}",
    custom_xticks::Union{Nothing, Vector{Float64}} = nothing,
    custom_yticks::Union{Nothing, Vector{Float64}} = nothing
)
    # Apply smoothing and scaling in one step
    x_scaled = smooth_and_scale(x_vals, log_x, smooth_x)
    y_scaled = smooth_and_scale(y_vals, log_y, smooth_y)

    # Apply log10 transformation to z_vals if log_z is true
    z_transformed = log_z ? log10.(z_vals) : z_vals

    # Extract unique x and y values **after** both smoothing and scaling
    x_unique = unique(x_scaled)
    y_unique = unique(y_scaled)

    # Create a z-matrix for contour plotting
    z_matrix = [NaN for _ in 1:length(x_unique), _ in 1:length(y_unique)]
    for i in 1:length(x_scaled)
        x_idx = findfirst(==(x_scaled[i]), x_unique)
        y_idx = findfirst(==(y_scaled[i]), y_unique)
        z_matrix[x_idx, y_idx] = z_transformed[i]
    end

    # Generate custom x and y ticks if specified
    x_ticks = if custom_xticks === nothing
        (x_unique, [x == round(x) ? string(Int(x)) : string(x) for x in unique(x_vals)])
    else
        xticks_transformed = smooth_and_scale(custom_xticks, log_x, smooth_x)
        (xticks_transformed, [x == round(x) ? string(Int(x)) : string(x) for x in custom_xticks])
    end
    
    y_ticks = if custom_yticks === nothing
        (y_unique, [y == round(y) ? string(Int(y)) : string(y) for y in unique(y_vals)])
    else
        yticks_transformed = smooth_and_scale(custom_yticks, log_y, smooth_y)
        (yticks_transformed, [y == round(y) ? string(Int(y)) : string(y) for y in custom_yticks])
    end

    # Plot the contour
    contour(
        x_unique,
        y_unique,
        z_matrix',
        xlabel = xlabel,
        ylabel = ylabel,
        xticks = x_ticks,
        yticks = y_ticks,
        color = palette,  # Colormap for trends
        fill = true,       # Filled contours
        legend = :topright,
        tickfontsize = 16,
        guidefontsize = 20,
        right_margin = 15 * Plots.mm,
    )
end


###################
# Plot crosstable #
###################
function crossTableScatterPlot(
    result1::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}},
    result2::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}},
    xlabel::String = "Metric 1",
    ylabel::String = "Metric 2",
    title::String = "Cross Table Scatter Plot"
)
    # Unpack results
    x_vals1, y_vals1, z_vals1 = result1
    x_vals2, y_vals2, z_vals2 = result2

    # Map (x, y) -> z for both result sets
    z_map1 = Dict(zip(zip(x_vals1, y_vals1), z_vals1))
    z_map2 = Dict(zip(zip(x_vals2, y_vals2), z_vals2))

    # Find common (x, y) pairs
    common_keys = intersect(keys(z_map1), keys(z_map2))

    # Extract matched z-values
    z1 = [z_map1[key] for key in common_keys]
    z2 = [z_map2[key] for key in common_keys]

    # Plot scatter plot
    scatter(
        z1, z2,
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        marker = (:circle, 6, :blue),
        legend = false
    )
end


function computeParetoSkyline(points::Vector{Tuple{Float64, Float64}})
    # Sort by density descending, then conductance ascending
    sorted_points = sort(points, by = x -> (-x[1], x[2]))

    # Compute Pareto frontier
    pareto_front = []
    current_best_c = Inf  # Track the best (minimum) conductance

    for (d, c) in sorted_points
        if c < current_best_c  # Strictly better conductance
            push!(pareto_front, (d, c))
            current_best_c = c
        end
    end

    return pareto_front
end


function findParetoOptimalXY(
    density_results::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}},
    conductance_results::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
)
    # Unpack results
    x_vals_dens, y_vals_dens, z_vals_dens = density_results
    x_vals_cond, y_vals_cond, z_vals_cond = conductance_results

    # Map (x, y) -> z for both result sets
    z_map_dens = Dict(zip(zip(x_vals_dens, y_vals_dens), z_vals_dens))
    z_map_cond = Dict(zip(zip(x_vals_cond, y_vals_cond), z_vals_cond))

    # Find common (x, y) pairs
    common_keys = intersect(keys(z_map_dens), keys(z_map_cond))

    # Extract matched (x, y) pairs with their (density, conductance) values
    pareto_candidates = [(key[1], key[2], z_map_dens[key], z_map_cond[key]) for key in common_keys]

    # Compute Pareto frontier
    pareto_front = computeParetoSkyline(collect(Tuple{Float64, Float64}, (d, c) for (_, _, d, c) in pareto_candidates))

    # Find (x, y) pairs corresponding to Pareto-optimal points
    pareto_xy = [(x, y) for (x, y, d, c) in pareto_candidates if (d, c) in pareto_front]

    return pareto_xy
end


function compareCrossTableScatterPlots(
    result_pairs::Vector{Tuple{
        Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}},
        Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    }},
    labels::Vector{String},
    xlabel::String = "Density",
    ylabel::String = "Conductance",
    title::String = "Cross Table Comparison with Skyline"
)
    # Define a color palette for different result pairs
    palette = distinguishable_colors(length(result_pairs))

    # Initialize plot
    plt = scatter(
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        legend = :topright
    )

    all_pareto_points = []  # Store all Pareto-optimal points

    # Iterate over result pairs and plot them
    for (i, (result1, result2)) in enumerate(result_pairs)
        # Unpack results
        x_vals1, y_vals1, z_vals1 = result1
        x_vals2, y_vals2, z_vals2 = result2

        # Map (x, y) -> z for both result sets
        z_map1 = Dict(zip(zip(x_vals1, y_vals1), z_vals1))
        z_map2 = Dict(zip(zip(x_vals2, y_vals2), z_vals2))

        # Find common (x, y) pairs
        common_keys = intersect(keys(z_map1), keys(z_map2))

        # Extract matched z-values
        z1 = [z_map1[key] for key in common_keys]  # Density
        z2 = [z_map2[key] for key in common_keys]  # Conductance

        # Collect all points for Pareto frontier computation
        pareto_points = [(d, c) for (d, c) in zip(z1, z2)]
        append!(all_pareto_points, pareto_points)

        # Add points to the scatter plot with a unique color
        scatter!(
            z1,
            z2,
            marker = (:circle, 6, palette[i]),
            label = labels[i]
        )
    end

    # Compute and plot Pareto skyline
    if !isempty(all_pareto_points)
        pareto_front = computeParetoSkyline(collect(Tuple{Float64, Float64}, all_pareto_points))

        # Extract sorted Pareto-optimal points for plotting
        skyline_x = [p[1] for p in pareto_front]
        skyline_y = [p[2] for p in pareto_front]

        # Overlay Pareto frontier
        plot!(
            skyline_x,
            skyline_y,
            linewidth = 2,
            linecolor = :red,
            label = "Pareto Skyline"
        )
    end

    return plt
end

