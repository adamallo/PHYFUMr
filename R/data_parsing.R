# R/data_parsing.R

#' Reads b-values from a Phyfum XML file
#'
#' @param xml_file xml file
#' @param matrix whether the return should be a data table (column by sample) or an alignment-like matrix (column by fCpG). All phyfumr functions expect matrix=T (Default)
#'
#' @returns b-value dataset in data.table or matrix form
#' @export

get_Bdata_from_XML <- function(xml_file,matrix=T){
  if(!file.exists(xml_file)) {
    warning(paste0(xml_file," does not exits"))
    return(NULL)
  }

  xmlContent <- XML::xmlToList(XML::xmlParse(xml_file))

  bData <- data.table::rbindlist(lapply(xmlContent$afalignment,function(thisAF){
    if(is.list(thisAF) && !is.null(thisAF$text)){
      data.table::data.table(b=as.numeric(unlist(strsplit(split = ",",x = gsub(pattern = "\n *",replacement = "",thisAF$text)))),sample=as.numeric(thisAF$taxon))}
  }))

  if(matrix==T) {
    bData[, locus := (seq_len(.N)),by=sample]
    bData <- t(as.matrix(data.table::dcast(bData,locus~sample,value.var = "b")[,`:=`(locus=NULL)]))
  }

  return(bData)
}

#' Caches b-values stored in a list of CSV files
#'
#' Do not use this function in parallel.
#'
#' @details CSV files may contain data from multiple patients but the CSV
#'   parsing is done by patient [get_Bdata_from_CSVs()]. To avoid re-loading the
#'   same files multiple times, this function can be used first to load all
#'   files in memory and then parse them efficiently. This function must not be
#'   used in parallel since it stores the data in a package environment. All
#'   downstream analyses may be conducted in parallel.
#'
#' @param csvs list of csv files to parse
#'
#' @returns NULL
#' @seealso [get_Bdata_from_CSVs()] [delete_cached_CSV_data()]
#' @export

load_Bdata_from_CSVs <- function(csvs){
  normalized_files <- normalizePath(csvs)

  if(any(!file.exists(normalized_files) | dir.exists(normalized_files)))
    stop(paste0("The requested files: ", paste0(collapse=",",normalized_files[!file.exists(normalized_files) | dir.exists(normalized_files)])," were not found or are not valid files"))

  for (csv in normalized_files) {
    if(is.null(.phyfumr_env$loadedCSVs[[csv]])){
      .phyfumr_env$loadedCSVs[[csv]] <- data.table::setnames(data.table::fread(csv),1,"CpG")
    } else {
      warning(paste0("Reusing previously-loaded CSV data from ",csv))
    }
  }

  return(NULL)
}

#' Caches b-values stored in all CSV files in a list of directories
#'
#' @inherit load_Bdata_from_CSVs details return seealso
#' @param csv_dirs list of directories for which to read all CSV files
#'
#' @seealso [load_Bdata_from_CSVs()]
#' @export

load_Bdata_from_CSV_dirs <- function(csv_dirs){
  normalized_dirs <- normalizePath(csv_dirs)

  if(any(!dir.exists(normalized_dirs)))
    stop(paste0("The requested directories: ", paste0(collapse=",",normalized_dirs[!dir.exists(normalized_dirs)])," were not found"))

  for (csv_dir in normalized_dirs) {
    theseCSVs <- list.files(csv_dir,pattern = ".csv$",full.names = T)
    load_Bdata_from_CSVs(theseCSVs)
    .phyfumr_env$loadedCSVs[[csv_dir]] <- Reduce(function(x,y) merge(x,y,by="CpG",all=T),.phyfumr_env$loadedCSVs[theseCSVs])
  }

  return(NULL)
}

#' Deletes CSV data cached by [load_Bdata_from_CSVs]
#'
#' @returns NULL
#' @seealso [get_Bdata_from_CSVs()] [load_Bdata_from_CSVs()]
#' @export

delete_cached_CSV_data <- function(){
  rm(list=c("loadedCSVs"),envir=.phyfumr_env)
}

#' Obtains b-values from CSV files (cached or not)
#'
#' @details CSV files may contain data from multiple patients but the CSV
#'   parsing is done by patient. To avoid re-loading the same files multiple
#'   times, [load_Bdata_from_CSVs()] or [load_Bdata_from_CSV_dirs()] can be used
#'   first to load all files in memory and then parse them efficiently here.
#'   This function can be used in parallel if the data has been previously
#'   cached. This function assumes it is executed in parallel and stops when it
#'   is not safe to do so. If the user is certain the function is not run in
#'   parallel, isParallel should be set to F. The cached data can be deleted
#'   using [delete_cached_CSV_data()]
#'
#' @param csvs CSV files to parse
#' @param patient_id patient ID to obtain the CSV column names to parse from
#'   sample_pairing
#' @param sample_pairing data.table with a patient column and a samplename
#'   column.
#' @inherit get_Bdata_from_XML params return
#' @param isParallel True by default, assuming this function is executed in
#'   parallel and warning when this is not safe to do.
#'
#' @seealso [load_Bdata_from_CSVs()] [load_Bdata_from_CSV_dirs()]
#'   [delete_cached_CSV_data]
#' @export

get_Bdata_from_CSVs <- function(csvs,patient_id,sample_pairing,matrix=T,isParallel=T){
  normalized_files <- normalizePath(csvs)
  if(any(!file.exists(normalized_files) & dir.exists(normalized_files)))
    stop(paste0("The requested files: ", paste0(collapse=",",normalized_files[!file.exists(normalized_files) | dir.exists(normalized_files)])," were not found or are not valid files"))

  is_loaded <- !sapply(.phyfumr_env$loadedCSVs[normalized_files],is.null)
  if(!all(is_loaded)){
    if(isParallel==F){
      load_Bdata_from_CSVs(normalized_files[!is_loaded])
    } else {
      stop(paste0("This function should not be run in parallel",
                  " unless the data has been pre-loaded with load_Bdata_from_CSVs.",
                  " If you are not executing this function in parallel, ",
                  "feel free to use the isParallel=F argument and execute it again"))
    }
  }

  is_loaded <- !sapply(.phyfumr_env$loadedCSVs[normalized_files],is.null)
  if(!all(is_loaded)){
    warning(paste0("CSV files not found: ", paste0(collapse=",",normalized_files[!is_loaded])))
    return(NULL)
  }
  if(length(sample_pairing[patient==patient_id,samplename])==0){
    warning(paste0("No patient information for patient_id ",patient_id))
    return(NULL)
  }

  theseData <- Reduce(function(x,y) merge(x,y,by="CpG",all=T),lapply(normalized_files,function(thisfile){
                      .phyfumr_env$loadedCSVs[[thisfile]][,.SD,.SDcols=sample_pairing[patient==patient_id,samplename]]}))
  cleanData <- theseData[rowSums(!is.na(theseData)) > 0,]
  if(matrix==T) return(t(as.matrix(cleanData)))
  return(data.table::melt(cleanData,variable.name = "sample",value.name = "b")[,`:=`(sample=as.character(sample))])
}

#' Obtains b-values from CSV files (cached or not)
#'
#' @inherit get_Bdata_from_CSVs details params return seealso
#' @param csv_dirs list of directories for which to read all CSV files
#'
#' @seealso [get_Bdata_from_CSVs()]
#' @export

get_Bdata_from_CSV_dirs <- function(csv_dirs,patient_id,sample_pairing,matrix=T,isParallel=T){
  normalized_dirs <- normalizePath(csv_dirs)
  if(any(!dir.exists(normalized_dirs)))
    stop(paste0("The requested directories: ", paste0(collapse=",",normalized_dirs[!dir.exists(normalized_dirs)])," were not found"))

  is_loaded <- !sapply(.phyfumr_env$loadedCSVs[normalized_dirs],is.null)

  if(!all(is_loaded)){
    if(isParallel==F){
      load_Bdata_from_CSV_dirs(normalized_dirs[!is_loaded])
    } else {
      stop(paste0("This function should not be run in parallel",
                  " unless the data has been pre-loaded with load_Bdata_from_CSV_dirs.",
                  " If you are not executing this function in parallel, ",
                  "feel free to use the isParallel=F argument and execute it again"))
    }
  }

  is_loaded <- !sapply(.phyfumr_env$loadedCSVs[normalized_dirs],is.null)
  if(!all(is_loaded)){
    warning(paste0("CSV files in the following directories cannot be loaded: ", paste0(collapse=",",normalized_dirs[!is_loaded])))
    return(NULL)
  }
  if(length(sample_pairing[patient==patient_id,samplename])==0){
    warning(paste0("No patient information for patient_id ",patient_id))
    return(NULL)
  }

  theseData <- Reduce(function(x,y) merge(x,y,by="CpG",all=T),lapply(csv_dirs,function(thisdir){
    .phyfumr_env$loadedCSVs[[thisdir]][,.SD,.SDcols=sample_pairing[patient==patient_id,samplename]]}))
  cleanData <- theseData[rowSums(!is.na(theseData)) > 0,]
  if(matrix==T) return(t(as.matrix(cleanData)))
  return(data.table::melt(cleanData,variable.name = "sample",value.name = "b")[,`:=`(sample=as.character(sample))])
}

#' Parses a phyfum .log file and summarizes a list of parameters
#'
#' @param log_file name of the log file to parse
#' @param burnin proportion of samples to be removed as burnin
#' @param sumfun summary function. One of mean or median.
#' @param params list of parameters to parse
#'
#' @returns data.table with a a column per param with its point estimate
#' @export

estimate_params <- function(log_file,
                      burnin=0.1,
                      sumfun=c("mean","median"),
                      params=c("errorModel.deltaOffset",
                               "errorModel.etaOffset",
                               "errorModel.kappaScale",
                               "flipflop.mu",
                               "flipflop.gamma",
                               "flipflop.lambda",
                               "alignment.stemCells")) {

  sumfun <- match.arg(sumfun)

  if(!file.exists(log_file))
    stop(paste0(log_file," does not exits"))

  b_data <- data.table::fread(log_file)
  b_data[, valid := .I/.N>burnin]

  return(as.list(b_data[valid==T,lapply(.SD,sumfun),.SDcols=params]))
}
