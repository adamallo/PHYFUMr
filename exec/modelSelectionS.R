#!/usr/bin/env Rscript

library(phyfumr)
#Config that is not expected to change
#######################################
basename_regex <- "\\.[0-9]+cells\\..*"

#DEBUG
#######################################
# n_cells_regex  <- ".*S([0-9]+).*" #DEBUG version
# basename_regex <- "S[0-9]+.*" #DEBUG version
# baseDir <- "~/Downloads/short_worflow_res/test"
# inDir <- "~/Downloads/short_worflow_res/test"
# args <- c(baseDir, paste(sep="/",inDir,"169flipflopS3/169flipflopS3.MLE.csv"),
#           paste(sep="/",inDir,"169flipflopS4/169flipflopS4.MLE.csv"),
#           paste(sep="/",inDir,"169flipflopS5/169flipflopS5.MLE.csv"),
#           paste(sep="/",inDir,"169flipflopS6/169flipflopS6.MLE.csv"),
#           paste(sep="/",inDir,"169flipflopS7/169flipflopS7.MLE.csv"),
#           paste(sep="/",inDir,"169flipflopS8/169flipflopS8.MLE.csv"))
#######################################

#Parsing command-line arguments and checking input
#######################################
args <- commandArgs(trailingOnly = TRUE) #Comment for debugging, uncomment to run with Rscript

if (length(args) <3) {
  stop("Usage: script outDir runNameSX.csv ... runNameSN.csv")
} else {
  ##ForPablo: Remove if we don't want this to be verbose
  print(paste0("Running the script in the MLE summaries (",paste(collapse =", ",args[-c(1)]),")"))
}

baseDir=args[1]
files=args[-c(1)]
baseName <- gsub(x = basename(files[1]),pattern = basename_regex,replacement = "")
#######################################

#Output config
#######################################
#Needed because in phyfumflow we are running this script by patient
plot_output_suffix <- paste(sep=".",baseName,"allS.pdf")
selected_output_suffix <- paste(sep=".",baseName,"selectedS.csv")
bestS_output_suffix <- paste(sep=".",baseName,"bestFitS.csv")
all_output_suffix <- paste(sep=".",baseName,"allS.csv")
extreme_selected_output_suffix <- paste(sep=".",baseName,"extremeSelectedS.csv")
extreme_valid_output_suffix <- paste(sep=".",baseName,"extremeValidS.csv")
#######################################

invisible(model_selection(in_dir = paste(sep="/",mountDir,analysesFolder),
                        out_dir = baseDir,
                        plot_output_suffix = plot_output_suffix,
                        selected_output_suffix = selected_output_suffix,
                        bestS_output_suffix = bestS_output_suffix,
                        all_output_suffix = all_output_suffix,
                        extreme_selected_output_suffix  = extreme_selected_output_suffix ,
                        extreme_valid_output_suffix = extreme_valid_output_suffix))
