library(phyfumr)
baseDir <- "~/projects/flipFlop/latestG22/tissueUnspecificSelection/nf"

runsFolder <- "runs"
#runsFilename <- "runs.tar.gz"
#mountDir <- "~/ratarmounts/"
#tarFile <- paste(sep="/",baseDir,runsFilename)

out_dir <- paste(sep="/",baseDir,"analysesFinal/overS")
out_tree_dir <- out_dir
#out_trace_dir <- paste(sep="/",baseDir,"marginalizedRuns")

dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)

#system(paste("ratarmount -u",mountDir)) ##Just in case something was mounted before
#system(paste("ratarmount",tarFile,mountDir))


#all_results <- model_average(file_dir = paste(sep="/",mountDir,runsFolder),
#Single function
all_results <- model_average(file_dir = paste(sep="/",baseDir,runsFolder),
                             out_dir = out_dir,
                             out_trace_dir = NULL,
                             out_tree_dir = out_tree_dir,
                             mle_file = "~/projects/flipFlop/latestG22/tissueUnspecificSelection/nf/analyses/MLE_filtered.csv",
                             n_cells_regex = ".*s([0-9]+).*",
                             basename_regex = "_s[0-9]+.*",
                             n_cores = 1,
                             treeannotator_args = c("-type","mcc"))

#Single patient
trees_files <- list.files(paste(sep="/",baseDir,runsFolder),pattern = "EAN.*.trees$",recursive = T,full.names = T)
single_result <- model_average_patient(trees_files = trees_files,
                                       out_dir = out_dir,
                                       out_trace_dir = NULL,
                                       out_tree_dir = out_tree_dir,
                                       mle_file = "~/projects/flipFlop/latestG22/tissueUnspecificSelection/nf/analyses/MLE.csv",
                                       n_cells_regex = ".*s([0-9]+).*",
                                       basename_regex = "_s[0-9]+.*",
                                       n_cores = 2,
                                       treeannotator_args = c("-type","mcc"))
#All Distributed
#Getting dataset information
trees_files <- list.files(path = paste(sep="/",baseDir,runsFolder), pattern = "*.trees$",
                          full.names = TRUE,
                          recursive = TRUE)

trees_files <- trees_files[file.size(trees_files)!=0] #Do not consider empty traces

run_names <- basename(trees_files)
runs_info <- data.table::data.table(trees_files=trees_files,run_names = gsub(pattern = ".trees",
                                                                             replacement = "",
                                                                             x = run_names))
runs_info[,`:=`(patient=gsub("_s[0-9]+.*","",run_names))]

for (this_patient in runs_info[,unique(patient)]){
  trees_files <- list.files(paste(sep="/",baseDir,runsFolder),pattern = paste0(this_patient,".*.trees$"),recursive = T,full.names = T)
  invisible(model_average_patient(trees_files = trees_files,
                                         out_dir = out_dir,
                                         out_trace_dir = NULL,
                                         out_tree_dir = out_tree_dir,
                                         mle_file = "~/projects/flipFlop/latestG22/tissueUnspecificSelection/nf/analyses/MLE.csv",
                                         n_cells_regex = ".*s([0-9]+).*",
                                         basename_regex = "_s[0-9]+.*",
                                         n_cores = 2,
                                         treeannotator_args = c("-type","mcc")))
}

#Gather after distribution
all_results_gathered <- model_average_gather(out_dir = out_dir,
                             n_cells_regex = ".*s([0-9]+).*",
                             basename_regex = "_s[0-9]+.*")

#Reformat to old format
library(data.table)
summary_long <- all_results$ptable
setnames(summary_long, old=c("meanP"), new=c("mean"))
summary_fully_long <- melt(summary_long, id.vars = c("patient","param"), variable.name = "stat")
summary_wide <- dcast(summary_fully_long[stat%in%c("mean","HDIUpper","HDILower")],patient~stat+param,sep = ".",value.var="value")
summary_wide[,`:=`(full_patient=patient)]
summary_wide[,`:=`(patient=sub("^([^_]*)_.*","\\1",full_patient),location=sub("^[^_][^_]*_([^_][^_]*)_.*","\\1",full_patient))]
summary_wide[,`:=`(location=factor(location,levels=c("Colon","SmallIntestine","Endometrium","multiple"),labels=c("Colon","Small Intestine","Endometrium","multiple")))]
check_write_csv(summary_wide,
                "summaryWithHDI.csv",
                out_dir)

#system(paste("ratarmount -u",mountDir))
