library(phyfumr)
baseDir <- "~/projects/flipFlop/latestG22/tissueUnspecificSelection/nf"

runsFolder <- "runs"
#runsFilename <- "runs.tar.gz"
#mountDir <- "~/ratarmounts/"
#tarFile <- paste(sep="/",baseDir,runsFilename)

out_dir <- paste(sep="/",baseDir,"analyses/overS")
out_tree_dir <- out_dir
#out_trace_dir <- paste(sep="/",baseDir,"marginalizedRuns")

dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)

#system(paste("ratarmount -u",mountDir)) ##Just in case something was mounted before
#system(paste("ratarmount",tarFile,mountDir))


#all_results <- model_average(file_dir = paste(sep="/",mountDir,runsFolder),
all_results <- model_average(file_dir = paste(sep="/",baseDir,runsFolder),
                             out_dir = out_dir,
                             out_trace_dir = NULL,
                             out_tree_dir = out_tree_dir,
                             mle_file = "~/projects/flipFlop/latestG22/tissueUnspecificSelection/nf/analyses/MLE.csv",
                             n_cells_regex = ".*s([0-9]+).*",
                             basename_regex = "_s[0-9]+.*",
                             n_cores = 2,
                             treeannotator_args = c("-type","mcc"))

#system(paste("ratarmount -u",mountDir))
