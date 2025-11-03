# R/traces.R

# Simplified version of RWTY traces that does not unroot trees

#' Ensures that traces are in the expected format
#'
#' @param traces MCMC data. List of [load_traces()] values
#'
#' @returns List of [load_traces()] values
#' @keywords internal

ensure_list_traces <- function(traces){
  if(methods::is(traces,"phyfumr_trace"))
    traces <- list(traces)
  if(is.null(names(traces)))
    data.table::setattr(traces,"name",seq(1,length(traces)))
  return(traces)
}

#' Removes burnin states in lists of posterior samples
#'
#' @inheritParams ensure_list_traces
#' @param burnin_p Proportion of MCMC samples to discard as burnin
#'
#' @returns list of PHYFUMr traces with burnin samples removed
#' @keywords internal

remove_burnin_traces <- function(traces,burnin_p) {
  traces <- ensure_list_traces(traces)
  return(lapply(traces,function(trace){
    last_sample <- length(trace$trees)
    first_sample <- ceiling(burnin_p*last_sample + .phyfumr_env[['precision']])
    trace$trees <- trace$trees[seq.int(from=first_sample,to=last_sample)]
    trace$ptable <- trace$ptable[seq.int(from=first_sample,to=last_sample),]
    trace
  }))
}

#' Fixes read NEXUS content from unfinished PHYFUM runs
#'
#' @param nexus_content result from readLines(nexus)
#'
#' @returns fixed content
#' @keywords internal

fix_nexus_content <- function(nexus_content){

  n_lines <- length(nexus_content)

  if(!grepl("END;|ENDBLOCK;", nexus_content[n_lines], ignore.case = TRUE)){
    nexus_content <- c(nexus_content,"End;")
  }
  return(nexus_content)
}

#' Loads PHYFUM MCMC traces
#'
#' Loads trees, looks for a log file of tree likelihoods and parameter values,
#' returns an phyfumr_trace object containing both
#'
#' @details Automatically fixes NEXUS files of unfinished runs on the fly
#'
#' @param treefile A path to a tree file containing an MCMC chain of trees
#' @param logfile A path to a file containing model parameters and likelihoods.
#'   If no path is provided, it looks for a file with the same name as the tree
#'   file with .log extension
#' @inheritParams remove_burnin_traces
#'
#' @return output A phyfumr_trace object containing the multiPhylo and the table
#'   of values from the log file if available.
#' @seealso \code{\link[ape]{read.tree}}, \code{\link[ape]{read.nexus}},
#'   \code{\link[rwty]{load.trees}}
#'
#' @export

load_trace <- function(treefile,logfile=NULL,burnin_p=0.1){
  if(!file.exists(treefile)) {
    stop("ERROR: Input tree file is not available")
  }
  trace <- list(trees = NULL, ptable = NULL)

  #trees
  con <- file(treefile, open = "r")
  trees_content <- readLines(con)
  close(con)

  fixed_trees_content <- fix_nexus_content(trees_content)

  trees_con <- textConnection(fixed_trees_content)
  trace$trees <- ape::read.nexus(trees_con)

  #ptable
  if(is.null(logfile)){
    logfile <- gsub(pattern = ".trees",replacement = ".log",treefile)
    if(!file.exists(logfile)){
      warning(sprintf("log file for %s not found. Parameter will not be loaded",treefile))
      logfile <- NULL
    }
  }
  if(!is.null(logfile)) {
    trace$ptable <- data.table::fread(logfile)
    if (length(trace$trees) != nrow(trace$ptable))
      stop("ERROR: .trees and .log files are not of the same length!")
  }

  class(trace) <- "phyfumr_trace"

  if(!is.null(burnin_p) && burnin_p > 0){
    return(remove_burnin_traces(trace,burnin_p)[[1]])
  } else {
    return(trace)
  }
}

#' Loads multiple PHYFUM MCMC traces
#'
#' @param trees_files list of PHYFUM's output .trees files (one per independent
#'   run)
#' @param log_files list of PHYFUM's output .log files (one per independent
#'   run). If not used, this function will assume they have the same name and
#'   location as their .trees counterparts
#' @inheritParams load_trace
#'
#' @returns list of RWTY's traces
#' @keywords internal

load_traces <- function(trees_files, log_files=NULL, burnin_p=0.1){
  if(!all(file.exists(trees_files))) {
    stop("ERROR: Not all input tree files exist. Problem files: ",paste(trees_files[!file.exists(trees_files)],collapse=", "))
  }
  if(is.null(log_files)){
    log_files <- gsub(pattern = ".trees",replacement = ".log",trees_files)
  }
  if(!all(file.exists(log_files))) {
    stop("ERROR: Not all input log files exist. Problem files: ",paste(log_files[!file.exists(log_files)],collapse=", "))
  }
  if(length(trees_files) != length(log_files)) {
    stop("ERROR: The trees and log input files are not paired properly.")
  }

  traces <- lapply(seq_along(trees_files),FUN=function(i_file){
    load_trace(trees_files[i_file],log_files[i_file],burnin_p)
  })

  return(traces)
}

#' Merges a list of PHYFUM traces into one
#'
#' @inheritParams load_traces
#'
#' @details It assumes samples are taken at regular intervals to perform the
#'   burnin
#'
#' @returns list with trees and ptable with the traces concatenated and the burnin states (trees
#'   and parameters) removed.
#' @keywords internal

merge_traces <- function(traces,burnin_p=0){
  return_trace <- traces[[1]]
  return_trace$trees <- do.call(c,lapply(traces,function(trace){
    last_sample <- length(trace$trees)
    first_sample <- ceiling(burnin_p*last_sample + .phyfumr_env[['precision']]) #min 1
    trace$trees[seq.int(from=first_sample,to=last_sample)]
  }))
  return_trace$ptable <- do.call(rbind,lapply(traces,function(trace){
    last_sample <- nrow(trace$ptable)
    first_sample <- ceiling(burnin_p*last_sample + .phyfumr_env[['precision']]) #min 1
    trace$ptable[seq.int(from=first_sample,to=last_sample),]
  }))
  return(return_trace)
}

#' Load PHYFUM posterior samples from the same condition
#'
#' If more than one run was generated per condition, they will be merged
#'
#' @inheritParams load_traces
#'
#' @returns list of traces (list of trees and ptable)
#' @keywords internal

load_condition_traces <- function(trees_files, log_files=NULL, burnin_p = 0.1) {
  traces <- load_traces(trees_files, log_files, burnin_p)
  if(length(traces)!=1 && !methods::is(traces,"phyfumr_trace")) {
    return(merge_traces(traces,burnin_p = 0)) #Burnin already happened at load_traces
  }
  return(traces)
}
