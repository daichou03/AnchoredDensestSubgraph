using Base
using StatsBase

PERFORMANCE_REPORTS_DIR = "../PerformanceReports/"
PERFORMANCE_REPORTS_INTEGRATED_DIR = "../PerformanceReportsIntegrated/"
DATA_POINTS_DIR = "../DataPoints/"
chosen_dataset_names = ["eucore","fbgov","epinion","livemocha"]
report_genre = ["time", "size"]
ALL_ALGORITHM_NAMES = ["ads", "iads", "slads"]

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
    files = filter(x->(length(x) >= 7) && x[7] == "AnchorSizeTest", files)
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
        for baseFile in baselineFileNames
            # Copy only SLADS
            io_read = open(string(PERFORMANCE_REPORTS_DIR, BaselineSubDir, baseFile))
            io_write = open(string(PERFORMANCE_REPORTS_DIR, ReportSubDir, baseFile), "w")
            for i = 1:2
                line = readline(io_read)
                write(io_write, string(split(line, ",")[3], "\n"))
            end
            close(io_read)
            close(io_write)
        end
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
    files = filter(x->((length(x) >= 7) && x[2][1:1] == "H") || (length(x) == 6), files)
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

function OutputIntegratedReportsByAlgorithm(ReportSubDir::String, ReportFileGroups::Array{Array{String,1},1}, OutputDir::String, ReportGenreIndex::Integer,
    AlgorithmNames::Vector{String}=ALL_ALGORITHM_NAMES)
    sep = (" ")
    for alg in 1:length(AlgorithmNames)
        dir = string(PERFORMANCE_REPORTS_INTEGRATED_DIR, OutputDir, "_", AlgorithmNames[alg], "_", report_genre[ReportGenreIndex], "/")
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
    reportFiles = GetBaselineReportFiles("/")
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

# -----------------------------------
# Integrate Data Points for SmallIADS
# -----------------------------------

function OutputIntegratedSmallIADSReports(DataPointSubDir::String)
    fileNames = readdir(string(DATA_POINTS_DIR, DataPointSubDir))
    files = map(x->split(x,"-"), fileNames)
    dir = string(PERFORMANCE_REPORTS_INTEGRATED_DIR, DataPointSubDir)
    mkpath(dir)
    io_write = open(string(dir,"overdensedR.txt"), "w")
    for i in 1:length(fileNames)
        dataName = files[i][1]
        tests = parse(Int64, files[i][2])
        Rsize = parse(Int64, files[i][5])
        # read from data point files
        io_read_dp = open(string(DATA_POINTS_DIR, DataPointSubDir, fileNames[i]))
        overdensed_sum = 0
        for j = 1:tests
            overdensed_sum += parse(Int64, split(readline(io_read_dp), ",")[5])
        end
        close(io_read_dp)
        # read from performance reports (assuming same sub folder and file name)
        io_read_pr = open(string(PERFORMANCE_REPORTS_DIR, DataPointSubDir, fileNames[i]))
        line = split(readline(io_read_pr), ",")
        IADS_time_ratio = parse(Float64, line[2]) / parse(Float64, line[1])
        close(io_read_pr)
        overdensed_mean = overdensed_sum / tests
        line_print = string(dataName, ",", Rsize, " ", overdensed_mean, " ", IADS_time_ratio) # Note need to convert overdensed_mean to non-overdensed prop later.
        write(io_write, string(line_print, "\n"))
    end
    close(io_write)
end

function AggregrateRSize(DataPointSubDir::String)
    fileNames = readdir(string(DATA_POINTS_DIR, DataPointSubDir))
    files = map(x->split(x,"-"), fileNames)
    for i in 1:length(fileNames)
        dataName = files[i][1]
        tests = parse(Int64, files[i][2])
        # read from data point files
        io_read_dp = open(string(DATA_POINTS_DIR, DataPointSubDir, fileNames[i]))
        size_R = zeros(tests)
        for j = 1:tests
            size_R[j] = parse(Int64, split(readline(io_read_dp), ",")[1])
        end
        close(io_read_dp)
        # read from performance reports (assuming same sub folder and file name)
        println(string(dataName, ",", StatsBase.mean(size_R), ",", StatsBase.std(size_R)))
    end
end

# ------------
# Do the thing
# ------------

function DoCollectResults()
    # --- baseline ---
    ReportSubDir="base/"
    ReportFiles = GetBaselineReportFiles(ReportSubDir)
    OutputIntegratedReport(ReportSubDir,ReportFiles,"20210331_baseline",1)
    OutputIntegratedReport(ReportSubDir,ReportFiles,"20210331_baseline",2)
    # --- anchorsize ---
    ReportSubDir = "ahs/" # Avoid long dir
    ReportFileGroups = GetAnchorSizeReportFileGroups(ReportSubDir)
    OutputIntegratedReportsByAlgorithm(ReportSubDir, ReportFileGroups, "20210331_anchorsize",1,["slads"])
    OutputIntegratedReportsByAlgorithm(ReportSubDir, ReportFileGroups, "20210331_anchorsize",2,["slads"]) 
    # --- halfedge ---
    ReportSubDir = "hl/"
    CopyBaselineReportsToHalfEdgeFolder(ReportSubDir, "base/")
    ReportFileGroups = GetHalfEdgeReportFileGroups(ReportSubDir)
    OutputIntegratedReportsByAlgorithm(ReportSubDir, ReportFileGroups, "20210331_halfedge",1,["slads"])
    OutputIntegratedReportsByAlgorithm(ReportSubDir, ReportFileGroups, "20210331_halfedge",2,["slads"])
end

# Example procedure of integrating half edge test results:
# ReportSubDir = "halfedge_20210322_SLADS_test/"
# ReportFileGroups = GetHalfEdgeReportFileGroups(ReportSubDir)
# OutputIntegratedReportsByAlgorithm(ReportSubDir, ReportFileGroups, "20210322_halfedge_test",1,["slads"])
# OutputIntegratedReportsByAlgorithm(ReportSubDir, ReportFileGroups, "20210322_halfedge_test",2,["slads"]) 

# function DoCollectResults()
#     OutputIntegratedReport(GetBaselineReportFiles(), "baseline", 1)
#     OutputIntegratedReport(GetBaselineReportFiles(), "baseline", 2)
#     OutputIntegratedReportsByAlgorithm(GetAnchorSizeReportFileGroups(), "anchorsize", 1)
#     OutputIntegratedReportsByAlgorithm(GetAnchorSizeReportFileGroups(), "anchorsize", 2)
#     OutputIntegratedReportsByAlgorithm(GetHalfEdgeReportFileGroups(), "halfedge", 1)
#     OutputIntegratedReportsByAlgorithm(GetHalfEdgeReportFileGroups(), "halfedge", 2)
# end
