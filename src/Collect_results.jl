PERFORMANCE_REPORTS_DIR = "../PerformanceReports/"
PERFORMANCE_REPORTS_INTEGRATED_DIR = "../PerformanceReportsIntegrated/"
chosen_dataset_names = ["eucore","fbgov","epinion","livemocha"]
report_genre = ["time", "size"]
algorithm_names = ["ads", "iads", "slads"]

function GetBaselineReportFiles()
    allReps = readdir(PERFORMANCE_REPORTS_DIR)
    files = filter(x->length(split(x,"-")) == 6, allReps)
end

# Ordered file name by groups
function GetAnchorSizeReportFileGroups(DatasetNames::Array{String,1}=chosen_dataset_names)
    allFiles = map(x->split(x,"-"), readdir(PERFORMANCE_REPORTS_DIR))
    files = filter(x->x[1] in DatasetNames && length(x) >= 7 && x[7] == "AnchorSizeTest", allFiles)
    fileGroups = map(x->filter(y->y[1] == x, files), DatasetNames)
    fileGroups = map(fileGroup->map(y->join(y, "-"), sort(fileGroup, by=x->parse(Int64, x[5]))), fileGroups)
end
# Ordered file name by groups
function GetHalfEdgeReportFileGroups(DatasetNames::Array{String,1}=chosen_dataset_names)
    allFiles = map(x->split(x,"-"), readdir(PERFORMANCE_REPORTS_DIR))
    files = filter(x->x[1] in DatasetNames && (length(x) == 6 || (length(x) >= 7 && x[2][1:1] == "H")), allFiles)
    fileGroups = map(x->filter(y->y[1] == x, files), DatasetNames)
    fileGroups = map(fileGroup->map(y->join(y, "-"), sort(fileGroup, by=x->length(x) * 100 + parse(Int64, x[2][2:2]))), fileGroups)
end

function OutputIntegratedReport(ReportFiles::Array{String,1}, OutputDir::String, ReportGenreIndex::Integer)
    sep = (" ")
    dir = string(PERFORMANCE_REPORTS_INTEGRATED_DIR, OutputDir, "_", report_genre[ReportGenreIndex], "/")
    mkpath(dir)
    io_write = open(string(dir,"fig.txt"), "w")
    for report in ReportFiles
        report_name = split(report, "-")[1]
        io_read = open(string(PERFORMANCE_REPORTS_DIR, report))
        line = ""
        for i = 1:ReportGenreIndex
            line = readline(io_read)
        end
        close(io_read)
        nums = split(line, ",")
        if ReportGenreIndex == 2
            nums = map(x->string(parse(Int64, x) ÷ 1048576), nums) # bytes to megabytes
        end       
        write(io_write, string(report_name, sep, join(nums, sep), "\n"))
    end
    close(io_write)
end

function OutputIntegratedReportsByAlgorithm(ReportFileGroups::Array{Array{String,1},1}, OutputDir::String, ReportGenreIndex::Integer)
    sep = (" ")
    for alg = 1:3
        dir = string(PERFORMANCE_REPORTS_INTEGRATED_DIR, OutputDir, "_", algorithm_names[alg], "_", report_genre[ReportGenreIndex], "/")
        mkpath(dir)
        io_write = open(string(dir,"fig.txt"), "w")
        for fileGroup in ReportFileGroups
            report_name = split(fileGroup[1], "-")[1]
            line_print = ""
            for report in fileGroup
                io_read = open(string(PERFORMANCE_REPORTS_DIR, report))
                line = ""
                for i = 1:ReportGenreIndex
                    line = readline(io_read)
                end
                close(io_read)
                num = split(line, ",")[alg]
                if ReportGenreIndex == 2
                    num = string(parse(Int64, num) ÷ 1048576) # bytes to megabytes
                end
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

OutputIntegratedReport(GetBaselineReportFiles(), "baseline", 1)
OutputIntegratedReport(GetBaselineReportFiles(), "baseline", 2)
OutputIntegratedReportsByAlgorithm(GetAnchorSizeReportFileGroups(), "anchorsize", 1)
OutputIntegratedReportsByAlgorithm(GetAnchorSizeReportFileGroups(), "anchorsize", 2)
OutputIntegratedReportsByAlgorithm(GetHalfEdgeReportFileGroups(), "halfedge", 1)
OutputIntegratedReportsByAlgorithm(GetHalfEdgeReportFileGroups(), "halfedge", 2)