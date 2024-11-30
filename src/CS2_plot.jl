using Plots
using Glob
include("Utils.jl")
include("Utils_io.jl")
include("LP_consts.jl")
include("Utils_graph.jl")
include("LP_evaluation.jl")

# Takes the length of the Parameterized LP results for each x (wAC) and y (wAD) value for specified anchor id.
# LP results for example are files in: (root)\LPCompResults\dblp-LPLAS-1.0-0.0-mytest.lpcompsets
function lpResultLengthTo2dCluster(dataName::String, suffixName::String, n::Int)
    # Initialize x, y, z arrays
    x_vals = []
    y_vals = []
    z_vals = []

    # Use Glob to find files matching the pattern
    files = GetParameterizedLPResultFileNames(dataName, suffixName, RESULT_TYPE_SETS)

    for file in files
        # Extract x and y values from the filename
        matching = match(r"-LPLAS-([0-9.]+)-([0-9.]+)-", file)
        if matching !== nothing
            x = parse(Float64, matching.captures[1])
            y = parse(Float64, matching.captures[2])
            y = y == 0 ? 0.00001 : y

            # Read file content and get n-th row's integer count
            rows = readlines(file)
            if n <= length(rows)
                integers = split(rows[n], ",")
                z = length(integers)
            else
                z = 0  # If the n-th row doesn't exist, default to 0
            end

            # Append x, y, z to the respective arrays
            push!(x_vals, x)
            push!(y_vals, y)
            push!(z_vals, z)
        end
    end
    return x_vals, y_vals, z_vals
end

function visualize_cluster((x_vals, y_vals, z_vals))
    # Convert to grid for surface plotting
    x_unique = unique(x_vals)
    y_unique = unique(y_vals)
    z_matrix = [0.0 for _ in 1:length(x_unique), _ in 1:length(y_unique)]

    for i in 1:length(x_vals)
        x_idx = findfirst(==(x_vals[i]), x_unique)
        y_idx = findfirst(==(y_vals[i]), y_unique)
        z_matrix[x_idx, y_idx] = z_vals[i]
    end

    # plot(
    #     x_unique,
    #     log10.(y_unique),
    #     z_matrix',
    #     st = :wireframe,  # :surface
    #     title = "Cluster Visualization",
    #     xlabel = "X-axis",
    #     ylabel = "Log10(Y-axis)",
    #     zlabel = "Z-axis",
    #     color = :rainbow
    # )
    # contour(
    #     x_unique,
    #     log10.(y_unique),
    #     z_matrix',
    #     title = "Cluster Contour Plot",
    #     xlabel = "X-axis",
    #     ylabel = "Log10(Y-axis)",
    #     color = :rainbow,  # Colormap for visualization
    #     clabels = true     # Show contour labels
    # )
    #y_ticks = [(val, "10^$(Int(round(log10(val), digits=0)))") for val in y_unique if val > 0]

    plt = scatter(
        x_vals,
        log10.(y_vals),
        marker_z = z_vals,  # Use z-values for color coding
        title = "2D Scatter Plot of result length",
        xlabel = "ω_{12}",
        ylabel = "log10 of ω_{24}",
        #yticks = y_ticks,
        color = :rainbow,   # Colormap
        legend = false,
        marker = (12, :circle),  # Marker size and style
    )
    for (x, y, z) in zip(x_vals, log10.(y_vals), z_vals)
        annotate!(x, y, text("$z", :white, 10))  # Add text annotation
    end
    return plt
end

