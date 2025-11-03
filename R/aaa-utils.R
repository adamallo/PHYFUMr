# R/aaa-utils.R

#' Uses pblapply or lapply depending on their availability
#'
#' @inheritParams  base::lapply
#' @inheritDotParams pbapply::pblapply cl
#'
#' @keywords internal

slapply <- function(X, FUN, ...) {
  if(is.null(.phyfumr_env[["use_parallel"]]) || !.phyfumr_env[["use_parallel"]]){
    ellipsis_args <- list(...)
    in_args <- c(ellipsis_args,list(X = X, FUN = FUN))
    in_args[["cl"]] <- NULL
    return(do.call(lapply,in_args))
  } else {
    return(pbapply::pblapply(X = X, FUN = FUN, ...))
  }
}

#' Gets the variable names from a list
#'
#' @param ... objects to get the name of
#'
#' @returns vector of names
#' @export

names_from_objects <- function(...) {
  as.character(substitute(list(...)))[-1]
}

#' Run a bundled script
#'
#' Used as Rscript -e 'phyfumr::run("scriptName.R")' --args arg1 arg2
#'
#' @param script Name of the script (file in `exec/`)
#' @export
#'

run <- function(script) {
  path <- system.file("exec", script, package = "phyfumr")
  if (path == "") stop("Script not found: ", script, call. = FALSE)

  # capture all args after --args
  args <- commandArgs(trailingOnly = TRUE)

  # run the bundled script with those args
  res <- system2(command = R.home("bin/Rscript"),
                 args = c(path, args))
  invisible(res)
}

#' Writes a CSV file if flag is not null
#'
#' @param table table to write in CSV
#' @param flag print if !is.null. Usually the filename or a filename part
#' @param outdir output directory. It will be made recursively if it does not exist, without warnings
#' @param filename name of the csv file. If null, flag is used
#' @param canceled_warning warning to produce if the file is not written
#'
#' @returns TRUE if the file was written, FALSE otherwise
#'
#' @keywords internal

check_write_csv <- function(table,flag,outdir,filename=NULL,canceled_warning=NULL){
  if(is.null(filename))
    filename <- flag

  if(!is.null(flag)) {
    if(!dir.exists(outdir))
      dir.create(outdir,showWarnings = FALSE,recursive = TRUE)
    utils::write.csv(table,file = paste(sep="/",outdir,filename),quote = FALSE,row.names = FALSE)
    return(TRUE)
  } else {
    if(!is.null(canceled_warning))
      warning(canceled_warning)
  }
  return(FALSE)
}
