# This file explores properties of a result from a single case with the same R but different \w_12 and \w_24 (x and y).
using Plots
using Glob
include("Utils.jl")
include("Utils_io.jl")
include("LP_consts.jl")
include("Utils_graph.jl")
include("LP_evaluation.jl")

Y_LOG_0_SMOOTH = 0.00001  # y-axis as log10. If y=0, to which value it is smoothed to for plotting.
REGEX_LPXY = string("-", SOLVER_NAMES[SOLVER_LP_ADSS],"-([0-9.]+)-([0-9.]+)-")


## Extract x, y, z values
function extractLPResultFileData(files::Vector{String}, n::Int, extractZFunction::Function)
    # Initialize x, y, z arrays
    data = []

    for file in files
        # Extract x and y values from the filename
        matching = match(Regex(REGEX_LPXY), file)
        if matching !== nothing
            x = parse(Float64, matching.captures[1])
            y = parse(Float64, matching.captures[2])

            # Compute z value using the provided function
            z = extractZFunction(file, n)

            # Append x, y, z to the respective arrays
            push!(data, (x, y, z))
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

function lpResultLengthTo2dCluster(dataName::String, suffixName::String, n::Int)
    files = GetParameterizedLPResultFileNames(dataName, suffixName, RESULT_TYPE_SETS)

    # Function to extract z-value (length of integers in the n-th row)
    function rowLengthZ(file::String, n::Int)
        return length(split(readlines(file)[n], ","))
    end

    return extractLPResultFileData(files, n, rowLengthZ)
end


function lpResultDistinctTo2dCluster(dataName::String, suffixName::String, n::Int)
    files = GetParameterizedLPResultFileNames(dataName, suffixName, RESULT_TYPE_SETS)
    unique_rows = Dict{String, Int}()

    function distinctRowZ(file::String, n::Int)
        rows = readlines(file)
        row_data = split(rows[n], ",")
        row_key = join(sort(row_data))  # Canonicalize row data
        if !haskey(unique_rows, row_key)
            unique_rows[row_key] = length(unique_rows) + 1
        end
        return unique_rows[row_key]
    end

    return extractLPResultFileData(files, n, distinctRowZ)
end


function lpResultStatsTo2dCluster(dataName::String, suffixName::String, n::Int, columnName::String)
    files = GetParameterizedLPResultFileNames(dataName, suffixName, RESULT_TYPE_STATS)

    # Function to extract z-value (column value for the n-th row)
    function columnValueZ(file::String, n::Int)
        return CSV.read(file, DataFrame)[n, columnName]  # natural exception if fails
    end

    return extractLPResultFileData(files, n, columnValueZ)
end


## Plot x, y, z values
function visualize_cluster((x_vals, y_vals, z_vals); smooth_y::Float64 = Y_LOG_0_SMOOTH, palette = :magma)
    # Smooth y-values: Replace 0 with the specified smooth_y value
    y_smoothed = [y == 0 ? smooth_y : y for y in y_vals]

    # Convert to grid for surface plotting
    x_unique = unique(x_vals)
    y_unique = unique(y_smoothed)
    z_matrix = [0.0 for _ in 1:length(x_unique), _ in 1:length(y_unique)]

    for i in 1:length(x_vals)
        x_idx = findfirst(==(x_vals[i]), x_unique)
        y_idx = findfirst(==(y_smoothed[i]), y_unique)
        z_matrix[x_idx, y_idx] = z_vals[i]
    end

    # Generate y-axis ticks that show the original y-values
    y_ticks = (log10.(y_smoothed), [string(y) for y in y_vals])

    # Scatter plot without explicit markers, and annotate the z values
    plt = scatter(
        x_vals,
        log10.(y_smoothed),          # Log-scaled smoothed y-values
        marker_z = z_vals,           # Use z-values for color coding
        title = "2D Scatter Plot of stats per x and y value",
        xlabel = "ω_{12}",
        ylabel = "ω_{24}",           # Show original y values
        yticks = y_ticks,            # Map log10 values to original labels
        color = palette,            # Colormap
        legend = false,
        markersize = 12               # Remove marker circles
    )

    # Annotate each point with z-values, formatted to 4 significant digits
    for (x, y, z) in zip(x_vals, log10.(y_smoothed), z_vals)
        annotate!(x, y, text(z, :white, 10))
    end

    return plt
end


function visualize_contour((x_vals, y_vals, z_vals); smooth_y::Float64 = Y_LOG_0_SMOOTH, palette = :magma, log_z::Bool = false)
    # Smooth y-values: Replace 0 with the specified smooth_y value
    y_smoothed = [y == 0 ? smooth_y : y for y in y_vals]

    # Apply log10 to z_vals if log_z is true
    z_transformed = log_z ? log10.(z_vals) : z_vals

    # Generate unique x and y coordinates
    x_unique = unique(x_vals)
    y_unique = unique(y_smoothed)

    # Create a z-matrix for contour plotting
    z_matrix = [NaN for _ in 1:length(x_unique), _ in 1:length(y_unique)]
    for i in 1:length(x_vals)
        x_idx = findfirst(==(x_vals[i]), x_unique)
        y_idx = findfirst(==(y_smoothed[i]), y_unique)
        z_matrix[x_idx, y_idx] = z_transformed[i]
    end

    # Plot the contour
    contour(
        x_unique,
        log10.(y_unique),  # Log-scaled y-values for better visualization
        z_matrix',
        title = "Contour Plot of Time Elapsed",
        xlabel = "ω_{12}",
        ylabel = "log10(ω_{24})",
        color = palette,  # Colormap for trends
        fill = true,       # Filled contours
        legend = :topright
    )
end

