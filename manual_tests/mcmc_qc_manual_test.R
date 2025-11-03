library(phyfumr)
baseDir <- "~/Downloads/short_worflow_res/test"
inDir <- "~/projects/flipFlop/beFirstRun/"
trees_files <- c(paste(sep="/",inDir,"666/169flipflopS8.trees"),paste(sep="/",inDir,"999/169flipflopS8.trees"))
out_dir <- "~/test/out_dir"
plot_dir <- "~/test/plot_dir"

convergence_single_patient <- mcmc_qc_condition(trees_files,out_dir = out_dir,plot_dir = plot_dir,n_cores = 8)
all_convergence <- mcmc_qc(file_dir = inDir,out_dir = out_dir, plot_dir = plot_dir, n_cores = 4, burnin_p = 0.1)
all_convergence_gathered <- mcmc_qc_gather(out_dir = out_dir)
