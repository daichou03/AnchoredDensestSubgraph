using Base

PERFORMANCE_REPORTS_DIR = "../PerformanceReports/"
PERFORMANCE_REPORTS_INTEGRATED_DIR = "../PerformanceReportsIntegrated/"
chosen_dataset_names = ["eucore","fbgov","epinion","livemocha"]
report_genre = ["time", "size"]
algorithm_names = ["ads", "iads", "slads"]

# Produce data for gnuplot from output files out of Query_test_yd.

# ---------------------
# Get Report File Names
# ---------------------

# 20210321: Specify the sub directory to retrieve reports.
# The following functions will assume ALL files under the target directory are eligible, with minimal filters that can just remove files like desktop.ini.

# Expected baseline report file name format example:
# eucore-1000-2-8-3-2
function GetBaselineReportFiles(ReportSubDir::String)
    files = map(x->split(x,"-"), readdir(string(PERFORMANCE_REPORTS_DIR, ReportSubDir)))
    files = filter(x->length(x) == 6, files)
    fileNames = map(y->join(y, "-"), files)
end

# Ordered file name by groups
# Expected anchor size report file name format example:
# eucore-1000-2-8-32-2-AnchorSizeTest
function GetAnchorSizeReportFileGroups(ReportSubDir::String)
    files = map(x->split(x,"-"), readdir(string(PERFORMANCE_REPORTS_DIR, ReportSubDir)))
    files = filter(x->(length(x) == 7) && x[7] == "AnchorSizeTest", files)
    datasetNames = unique!(map(x->x[1],files))
    fileGroups = map(x->filter(y->y[1] == x, files), datasetNames)
    fileGroups = map(fileGroup->sort(fileGroup, by=x->parse(Int64, x[5])), fileGroups)
    fileGroupNames = map(fileGroup->map(y->join(y, "-"), fileGroup), fileGroups)
end

# Try to find the baseline reports from another folder and copy to half edge report folder.
# Watch the number of datasets/files removed/copied to see if it was called properly.
function CopyBaselineReportsToHalfEdgeFolder(ReportSubDir::String, BaselineSubDir::String)
    halfEdgeFiles = map(x->split(x,"-"), readdir(string(PERFORMANCE_REPORTS_DIR, ReportSubDir)))
    halfEdgeFiles = filter(x->(length(x) == 7) && x[2][1:1] == "H", halfEdgeFiles)
    datasetNames = unique!(map(x->x[1],halfEdgeFiles))
    println(string(length(datasetNames), " unique datasets found."))
    baselineFiles = filter(y->y[1] in datasetNames, map(x->split(x,"-"), GetBaselineReportFiles(BaselineSubDir)))
    if length(baselineFiles) == 0
        println(string("Couldn't find matching baseline reports from folder: ", BaselineSubDir, ", abort action."))
    else
        # Delete all baseline reports in ReportSubDir first
        baselineFileNamesOld = map(x->string(PERFORMANCE_REPORTS_DIR, ReportSubDir, x), GetBaselineReportFiles(ReportSubDir))
        map(Base.Filesystem.rm, baselineFileNamesOld)
        println(string(length(baselineFileNamesOld), " old baseline reports removed."))
        # Copy
        baselineFileNames = map(y->join(y, "-"), baselineFiles)
        map(x->Base.Filesystem.cp(string(PERFORMANCE_REPORTS_DIR, BaselineSubDir, x), string(PERFORMANCE_REPORTS_DIR, ReportSubDir, x)), baselineFileNames)
        println(string(length(baselineFileNames), " baseline reports copied."))
    end
end

# Ordered file name by groups
# Assumes both half edge reports and their corresponding baseline reports are in the same folder.
# If not (as currently half edge tests DO NOT produce baseline reports in the same folder, by default this is the case), look at CopyBaselineReportsToHalfEdgeFolder.
# Expected half edge report file name format example:
# eucore-H1-1000-2-8-3-2
# Expected Baseline report file name format example:
# eucore-1000-2-8-3-2
function GetHalfEdgeReportFileGroups(ReportSubDir::String)
    files = map(x->split(x,"-"), readdir(string(PERFORMANCE_REPORTS_DIR, ReportSubDir)))
    files = filter(x->((length(x) == 7) && x[2][1:1] == "H") || (length(x) == 6), files)
    datasetNames = unique!(map(x->x[1],files))
    fileGroups = map(x->filter(y->y[1] == x, files), datasetNames)
    fileGroups = map(fileGroup->sort(fileGroup, by=x->length(x) * 100 + parse(Int64, x[2][2:2])), fileGroups)
    fileGroupNames = map(fileGroup->map(y->join(y, "-"), fileGroup), fileGroups)
end

# ------------------------
# Output Integrated Report
# ------------------------

function OutputIntegratedReport(ReportSubDir::String, ReportFiles::Array{String,1}, OutputDir::String, ReportGenreIndex::Integer)
    sep = (" ")
    dir = string(PERFORMANCE_REPORTS_INTEGRATED_DIR, OutputDir, "_", report_genre[ReportGenreIndex], "/")
    mkpath(dir)
    io_write = open(string(dir,"fig.txt"), "w")
    for report in ReportFiles
        if length(split(report, "-")) != 6 # Not a exclusive check, the error string itself is more explanatory.
            error(string("Unexpected report file name format: ", report, ", expected report file name format example: data-1000-2-8-3-2"))
        end
        report_name = split(report, "-")[1]
        tests = parse(Int64, split(report, "-")[2])
        io_read = open(string(PERFORMANCE_REPORTS_DIR, ReportSubDir, report))
        line = ""
        for i = 1:ReportGenreIndex
            line = readline(io_read)
        end
        close(io_read)
        nums = split(line, ",")
        nums = map(x->string(parse(Float64, x) / tests / (ReportGenreIndex == 2 ? 1048576 : 1)), nums) # For size report, bytes to megabytes
        write(io_write, string(report_name, sep, join(nums, sep), "\n"))
    end
    close(io_write)
end

function OutputIntegratedReportsByAlgorithm(ReportSubDir::String, ReportFileGroups::Array{Array{String,1},1}, OutputDir::String, ReportGenreIndex::Integer)
    sep = (" ")
    for alg = 1:3
        dir = string(PERFORMANCE_REPORTS_INTEGRATED_DIR, OutputDir, "_", algorithm_names[alg], "_", report_genre[ReportGenreIndex], "/")
        mkpath(dir)
        io_write = open(string(dir,"fig.txt"), "w")
        for fileGroup in ReportFileGroups
            report_name = split(fileGroup[1], "-")[1]
            line_print = ""
            for report in fileGroup
                tests_index = split(report, "-")[2][1:1] == "H" ? 3 : 2
                tests = parse(Int64, split(report, "-")[tests_index])
                io_read = open(string(PERFORMANCE_REPORTS_DIR, ReportSubDir, report))
                line = ""
                for i = 1:ReportGenreIndex
                    line = readline(io_read)
                end
                close(io_read)
                num = split(line, ",")[alg]
                num = string(parse(Float64, num) / tests / (ReportGenreIndex == 2 ? 1048576 : 1)) # For size report, bytes to megabytes
                if line_print != ""
                    line_print = string(line_print, sep)
                end
                line_print = string(line_print, num)
            end
            write(io_write, string(report_name, sep, line_print, "\n"))
        end
        close(io_write)
    end
end

# -----------
# Other stats
# -----------

function ReportSLADSPerformanceGain(ReportGenreIndex::Integer)
    sep = (" ")
    reportFiles = GetBaselineReportFiles()
    for report in reportFiles
        if length(split(report, "-")) != 6 # Not a exclusive check, the error string itself is more explanatory.
            error(string("Unexpected report file name format: ", report, ", expected report file name format example: data-1000-2-8-3-2"))
        end
        report_name = split(report, "-")[1]
        tests = parse(Int64, split(report, "-")[2])
        io_read = open(string(PERFORMANCE_REPORTS_DIR, report))
        line = ""
        for i = 1:ReportGenreIndex
            line = readline(io_read)
        end
        close(io_read)
        nums = map(x->parse(Float64, x), split(line, ","))
        println(string(report_name, " ", nums[2] / nums[3]))
    end
end

# ------------
# Do the thing
# ------------

# function DoCollectResults()
#     OutputIntegratedReport(GetBaselineReportFiles(), "baseline", 1)
#     OutputIntegratedReport(GetBaselineReportFiles(), "baseline", 2)
#     OutputIntegratedReportsByAlgorithm(GetAnchorSizeReportFileGroups(), "anchorsize", 1)
#     OutputIntegratedReportsByAlgorithm(GetAnchorSizeReportFileGroups(), "anchorsize", 2)
#     OutputIntegratedReportsByAlgorithm(GetHalfEdgeReportFileGroups(), "halfedge", 1)
#     OutputIntegratedReportsByAlgorithm(GetHalfEdgeReportFileGroups(), "halfedge", 2)
# end
