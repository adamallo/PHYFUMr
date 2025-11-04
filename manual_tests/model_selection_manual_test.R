library(phyfumr)
baseDir <- "~/projects/flipFlop/latestG22/tissueUnspecificSelection/nf"

analysesFolder <- "analyses"
#analysesFilename <- "analyses.tar.gz"
#mountDir <- "~/ratarmounts/"
mountDir <- baseDir
#tarFile <- paste(sep="/",baseDir,analysesFilename)

out_dir <- paste(sep="/",baseDir,"analyses")
plot_dir <- paste(sep="/",baseDir,"plots/modelSelection")

dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)
dir.create(plot_dir,showWarnings = FALSE,recursive = TRUE)

#system(paste("ratarmount -u",mountDir)) ##Just in case something was mounted before
#system(paste("ratarmount",tarFile,mountDir))

test <- model_selection(in_dir = paste(sep="/",mountDir,analysesFolder),
                        file_regex = "^MLE.csv$",
                        out_dir = out_dir,
                        plot_dir = plot_dir,
                        n_cells_regex = ".*s([0-9]+).*",
                        basename_regex = "_s[0-9]+.*")

# test <- model_selection(in_dir = paste(sep="/",mountDir,analysesFolder),
#                         #file_regex = "^MLE.csv$",
#                         out_dir = out_dir,
#                         plot_dir = plot_dir,
#                         n_cells_regex = ".*s([0-9]+).*",
#                         basename_regex = "_s[0-9]+.*")

#system(paste("ratarmount -u",mountDir))
