library(phyfumr)
baseDir <- "~/projects/flipFlop/latestG22/tissueUnspecificSelection/nf"

runFolder <- "runs"
runFilename <- "runs.tar.gz"
mountDir <- "~/ratarmounts/"
tarFile <- paste(sep="/",baseDir,runFilename)

out_dir <- paste(sep="/",baseDir,"analyses")
plot_dir <- paste(sep="/",baseDir,"plots")

dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)
dir.create(plot_dir,showWarnings = FALSE,recursive = TRUE)

system(paste("ratarmount -u",mountDir)) ##Just in case something was mounted before
system(paste("ratarmount",tarFile,mountDir))

trees_files <- paste(sep="/","~/projects/flipFlop/latestG22/tissueUnspecificSelection/nf/testruns",c("666/UT_Endometrium_85_s9.trees","999/UT_Endometrium_85_s9.trees"))

mcmc_qc_patient(trees_files,out_dir = out_dir,plot_dir = plot_dir,backup=F)
#all_convergence <- mcmc_qc(file_dir = paste0(mountDir,runFolder),out_dir = out_dir, plot_dir = plot_dir, n_cores = 8, burnin_p = 0.1)

system(paste("ratarmount -u",mountDir))
