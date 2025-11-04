# R/model_averaging.R

#TODO replace 1:length(XXX) to seq_along (1D objects) or seq_length(nrow/ncol) for multidimmensional objects across the package

#' Resample posterior using a sample map
#'
#' @param samples list of posterior sample tables
#' @param map data.table with two columns. \emph{x_i} indicates the posterior
#'   sample index and \emph{s_i} indicates the index within that sample (i.e.,
#'   `samples[[x_i]][s_i]`) for each of the final samples (rows)
#'
#' @returns posterior sample (array)
#' @keywords internal

apply_samplemap_table <- function(samples,map){
  data.table::rbindlist(mapply(function(x_i,s_i) samples[[x_i]][s_i,], map$x_i, map$s_i, SIMPLIFY = FALSE, USE.NAMES = FALSE))
}

#' Resample posterior using a sample map
#'
#' @param trees list of multiphylo objects with posterior samples of trees
#' @param map data.table with two columns. \emph{x_i} indicates the posterior
#'   sample index and \emph{s_i} indicates the index within that sample (i.e.,
#'   `samples[[x_i]][s_i]`) for each of the final samples (rows)
#'
#' @returns posterior tree sample (multiphylo)
#' @keywords internal

apply_samplemap_trees <- function(trees,map){
  tip_label <- unique(lapply(trees,attr,"TipLabel"))
  if(length(tip_label)!= 1)
    stop("Incompatible trees")

  out_trees <- mapply(`[`, trees[map$x_i], map$s_i, SIMPLIFY = TRUE, USE.NAMES = FALSE)
  class(out_trees) <- "multiPhylo"
  return(out_trees)
}

#' Calculate posterior probabilities given all log joint probabilities
#'
#' @param log_joint_probabilities array of the log joint probabilities (log marginal likelihoods +
#'   log prior) of each of the posterior samples (x)
#'
#' @return probabilities
#' @keywords internal

calculate_model_averaging_posteriors <- function(log_joint_probabilities) {
  log_marginal <- matrixStats::logSumExp(log_joint_probabilities,na.rm=T)
  log_probs <- log_joint_probabilities-log_marginal
  exp(log_probs)
}

#' Resample posteriors according to probabilities
#'
#' @param x list of posterior sample tables. They do not have to have the same
#'   length
#' @param probs array of sampling probabilities (one per posterior sample x)
#' @param n_marginalized_posterior_samples length of the resampled posterior. If
#'   NULL, it is maximized sampling without resampling. If more samples than
#'   existing are required, sampling is done with replacement, otherwise
#'   without.
#'
#' @returns a list with the posterior sample and a 2D array with a row per final
#'   posterior sample and column per input posterior sample, with elements
#'   indicating the ith position of the input posterior sample to be used as
#'   output, or 0 otherwise. Only one element >0 per row is expected
#' @keywords internal

resample_posterior <- function(x,
                               probs,
                               n_marginalized_posterior_samples = NULL){
  n_x <- length(x)

  if(any(!sapply(x,is.data.frame)))
    stop("ERRO: x must be a list of data frames")

  if(length(x) != length(probs))
    stop("ERROR: probs must be of the same length as traces")

  lengths <- sapply(x,nrow,simplify = TRUE)
  max_length <- max(lengths)

  if(is.null(n_marginalized_posterior_samples)){ #if no specified length, we maximize it
    which_x <- t(stats::rmultinom(sum(lengths),1,probs)) #resampling
    cum_samples <- apply(which_x,2,cumsum) #number of samples needed from each posterior sample (col) at each final possible length (row)
    n_marginalized_posterior_samples <- min(sapply(seq_len(ncol(cum_samples)),function(i){suppressWarnings(min(which(cum_samples[,i] == lengths[i])))})) #length of the maximum length possible given the sampling
    len_x <- cum_samples[n_marginalized_posterior_samples,] #final number of samples per posterior sample
    replace <- F
  } else {
    len_x <- as.numeric(stats::rmultinom(1,n_marginalized_posterior_samples,probs))
    if(all(len_x<=max_length)){ #without replacement
      replace <- F
    } else { #with replacement
      replace <- T
    }
  }

  structured_map <- NULL
  x_i <- rep(seq_len(n_x), times = len_x) # construct row indices indicating the original posterior sample (x)
  s_i <- unlist(lapply(seq_len(n_x), function(i) {
    sample(seq_len(lengths[i]), size = len_x[i], replace = replace)
  }), use.names = FALSE) # get the values
  structured_map <- data.table::data.table(x_i,s_i)
  structured_out_x <- apply_samplemap_table(x,structured_map)
  shuffled_is <- sample.int(n_marginalized_posterior_samples)

  return(list(sample=structured_out_x[shuffled_is,],
              map=structured_map[shuffled_is,]))
}

#TODO we may want to have version of this function tweaked to do bootstrapping over different lMLEs sampling them considering their uncertainty

#' Perform model-averaging across PHYFUMr traces from different models
#'
#' @details If the joint probabilities are not provided, the marginal likelihood
#'   is estimated here using the harmonic mean estimator. If the prior
#'   probabilities for each model are provided, they are used to calculate the
#'   final joint probabilities, otherwise, a uniform prior is implicitly used.
#'
#' @inheritParams ensure_list_traces
#' @param log_joint_probabilities array of the log joint probabilities (log
#'   marginal likelihoods + log prior) of each of the posterior samples (traces)
#' @param log_prior_probabilities array of log prior probabilities for each of
#'   the models to average across (traces)
#' @param posterior_probabilities array of the posterior probabilities of each
#'   (S) (i.e., resampling probabilities)
#' @inheritParams remove_burnin_traces
#' @inheritParams resample_posterior
#' @inheritParams check_model_averaging_posteriors
#'
#' @returns model-averaged chain

model_average_traces <- function(traces,
                                 log_joint_probabilities = NULL,
                                 log_prior_probabilities = NULL,
                                 posterior_probabilities = NULL,
                                 burnin_p = 0,
                                 n_marginalized_posterior_samples = 10000,
                                 max_expected_samples = 1,
                                 warn_biased_marginalization = TRUE){

  #Checking the data
  traces <- ensure_list_traces(traces)

  #Checking probability inputs
  check_probs <- function(prob,prob_name){
    if(!is.null(prob) && length(traces) != length(prob))
      stop(paste0("ERROR: ",prob_name," must be NULL or of the same length as traces"))
  }
  probs <- list(log_joint_probabilities,
                log_prior_probabilities,
                posterior_probabilities)
  probnames <- names_from_objects(log_joint_probabilities,
                                        log_prior_probabilities,
                                        posterior_probabilities)
  for (i in seq_along(probs)){
    check_probs(probs[[i]],probnames[i])
  }

  #Applying burnin
  if(!is.null(burnin_p)&&burnin_p>0){
    traces <- remove_burnin_traces(traces,burnin_p)
  }

  #Calculating posterior probabilities of S
  if(!is.null(posterior_probabilities)){ #Provided posteriors
    if(!is.null(log_prior_probabilities) || !is.null(log_joint_probabilities))
      warning("WARNING: log_prior_probabilities and/or log_joint_probabilities will be discarded since posterior probabilities were provided")
  } else {
    if (is.null(log_joint_probabilities)){
      ##estimate marginals using HME
      lMLEs <- sapply(traces,function(chain){
        if(is.null(chain[["ptable"]]) || is.null(chain[["ptable"]][,"likelihood"]))
          return(NaN)
        LaplacesDemon::LML(LL=chain$ptable$likelihood,method="HME")$LML})
      if(!is.null(log_prior_probabilities)){ #Using provided priors to calculate the joints
        log_joint_probabilities <- lMLEs + log_prior_probabilities
      } else { #Assuming uniform prior
        log_joint_probabilities <- lMLEs
        warning("lMLEs are calculated here and log prior probabilities for each model were not provided, thus, we are assuming an improper uniform prior")
      }
    }

    posterior_probabilities <- calculate_model_averaging_posteriors(log_joint_probabilities)
  }

  #Checking final posterior probabilities for high probabilities in extreme values of S
  if(warn_biased_marginalization){
    max_p <- max_expected_samples/n_marginalized_posterior_samples
    if(posterior_probabilities[1]>max_p)
      warning("Lowest S values with high posterior probability, the marginalized posteriors will be biased")
    if(posterior_probabilities[length(posterior_probabilities)]>max_p)
      warning("Highest S values with high posterior probability, the marginalized posteriors will be biased")
  }

  #Resampling
  ptables_resampling <- resample_posterior(x = lapply(traces,`[[`,"ptable"),
                                           probs = posterior_probabilities,
                                           n_marginalized_posterior_samples = n_marginalized_posterior_samples)
  out_trees <- apply_samplemap_trees(lapply(traces,`[[`,"trees"),ptables_resampling$map)

  #Formatting output
  periods <- sapply(traces,`[[`,"gens.per.tree",simplify = TRUE)
  if(length(unique(periods)) == 1) #if they are unique, we only indicate them once, otherwise, one per input chain
    periods <- unique(periods)

  resampled_traces <- list(trees=out_trees, ptable=ptables_resampling$sample, gens.per.tree = periods)
  class(resampled_traces) <- "phyfumr_trace"

  return(resampled_traces)
}

#' Calculate posterior probabilities given all log joint probabilities
#'
#' @param mle_file filename of the tabular output with marginal likelihood
#'   estimations generated by [mcmc_qc] and related functions. It can also be a
#'   list of files and all will be parsed
#' @param mle_method sorted list of preferred methods to estimate the marginal
#'   likelihood for model selection. Choose between PS (Path Sampling), SS
#'   (Stepping Stone), and HME (Harmonic Mean Estimator)
#' @inheritParams model_selection
#' @param log_prior_probabilities array of log prior probabilities for each of
#'   the models to average across, named using the corresponding value of S
#'   (S_XX).
#'
#' @return data.table with posterior probabilities for columns S_XX, with two
#'   additional columns indicating patient and method
#' @export

calculate_model_averaging_posteriors_file <- function(mle_file,
                                                      mle_method = c("PS","SS","HME"),
                                                      n_cells_regex = ".*\\.([0-9]+)cells\\..*",
                                                      basename_regex = "\\.[0-9]+cells\\..*",
                                                      log_prior_probabilities = NULL) {
  if(!all(sapply(mle_file,file.exists)))
    stop("MLE file(s) not found")

  ljp_table <- data.table::rbindlist(lapply(mle_file,data.table::fread),idcol = "file")
  cond_col <- grep("^cond",colnames(ljp_table),value=T)
  ljp_table[,`:=`(patient=gsub(basename_regex,"",get(cond_col)),S=as.numeric(gsub(n_cells_regex,"\\1",get(cond_col))))]

  selected_ljp_table <- data.table::rbindlist(lapply(split(ljp_table,by="patient"),function(this_data){
    n_S <- length(this_data[,unique(S)])
    method_analysis <- this_data[,.N,by="method"][,`:=`(valid_method=N==n_S)][mle_method,.SD,on="method"][valid_method==T,`:=`(method_order=.I)]
    selected_method <- method_analysis[valid_method==TRUE,][order(method_order),method]
    this_result <- data.table::dcast(this_data[method==selected_method,][order(S),.(patient,lML,S,method=selected_method)],patient+method~sprintf("S_%02d",S),value.var = "lML")
    return(this_result)
  }),fill = T)

  s_cols <- grep("^S",colnames(selected_ljp_table),value = TRUE)

  if(!is.null(log_prior_probabilities)){
    if(!is.null(names(log_prior_probabilities))){
      if(!all(s_cols %in% names(log_prior_probabilities)))
        stop("log_prior_probabilities provided but missing for some needed S values, or naming is not right. Expected format S_XX with XX as the 2-digit S value (e.g., S_04)")
      log_prior_probabilities <- log_prior_probabilities[s_cols]
    } else if (length(s_cols) == length(names(log_prior_probabilities))){
      warning("log_prior_probabilities is not named but of the same length as the number of S values. This function will assume log_prior_probabilities are provided in ascending value of S")
    }
    selected_ljp_table <- cbind(selected_ljp_table[,c(1,2)],sweep(as.matrix(selected_ljp_table[,.SD,.SDcols = s_cols]),2,log_prior_probabilities))
  }

  pp_table <- cbind(selected_ljp_table[,c(1,2)],selected_ljp_table[,t(apply(.SD,1,calculate_model_averaging_posteriors)),.SDcols = s_cols])
  return(pp_table)
}


#' Check model averaging posteriors for insufficient S sampling
#'
#' @param pp_table table with posterior probabilities of S as given by
#'   [calculate_model_averaging_posteriors_file]
#' @param n_posterior_samples number of desired posterior samples to calculate the number
#'   of expected samples per S
#' @param max_expected_samples maximum allowed expected number of samples in an
#'   extreme value of S
#' @param warn_biased_marginalization Logical, indicating whether to warn or not
#'   if the limit of
#'   max_expected_samples in an extreme value of S is violated
#'
#' @returns table with information on patients with insufficient S sampling
#' @export

check_model_averaging_posteriors <- function(pp_table,
                                             n_posterior_samples = 10000,
                                             max_expected_samples = 1,
                                             warn_biased_marginalization = TRUE) {

  pp_table_long <- data.table::melt(pp_table,value.name="pp",measure.vars = data.table::patterns("^S"),variable.name="S")
  pp_table_long[,`:=`(S=as.numeric(gsub("S_","",S)))]
  expectation_table <- pp_table_long[!is.na(pp),][order(S),.(pp_min_S=data.table::first(pp),pp_max_S=data.table::last(pp)),by=patient][,`:=`(expected_samples_min_S=pp_min_S*max_expected_samples,expected_samples_max_S=pp_max_S*n_posterior_samples)]
  return_table <- expectation_table[expected_samples_max_S>max_expected_samples | expected_samples_min_S>max_expected_samples,
                                    .(patient,expected_samples_min_S=round(expected_samples_min_S),expected_samples_max_S=round(expected_samples_max_S))]

  if(warn_biased_marginalization==T && nrow(return_table) > 0){
    #warning(paste(sep="\n","Extreme S values with high posterior probability\n"),
    warning(paste(sep="\n","Extreme S values with high posterior probability",
            paste(collapse="\n",utils::capture.output(print(expectation_table[expected_samples_max_S>max_expected_samples | expected_samples_min_S>max_expected_samples,
                                                   .(patient,expected_samples_min_S=round(expected_samples_min_S),expected_samples_max_S=round(expected_samples_max_S))]))
            )))
  }
  return(return_table)
}

#'Uses BEAST's treeannotator to summarize the posterior distribution of trees
#'from a file
#'
#' @details This function uses treeannotator's defaults except for the branch
#'  length summary function. This means that using different versions of
#'  treeannotator may result in different estimates (e.g., former versions
#'  defaulted to MCC trees while newer versions default to hipstr)
#'
#' @param trees_file trees files to summarize
#' @param burnin_p proportion of samples to be discarded as burnin
#' @param burnin_trees number of samples to be discarded as burnin
#'  (burnin_trees has priority)
#' @param outname name of the output summary tree file. NULL by default (not
#'  save it)
#' @param sum_fun Name of the function to use to summarize node heights
#' @param treeannotator_args An array of additional parameters to use when calling
#'  treeannotator.
#'
#' @return the summary tree as an ape::phylo
#' @export

get_summary_tree <- function(trees_file,
                             burnin_p = 0.1,
                             burnin_trees = NULL,
                             outname = NULL,
                             sum_fun = c("mean","median","keep","ca"),
                             treeannotator_args = NULL) {
  #Checking args
  sum_fun = match.arg(sum_fun)

  #Calculating burnin_trees if needed
  if(is.null(burnin_trees)){
    n_trees <- length(ape::read.nexus(trees_file))
    burnin_trees <- ceiling(n_trees*burnin_p)
  }

  #Run treeannotator
  result_tree <- NULL

  if(is.null(outname)) {
    command_args <- c(treeannotator_args,
                      "-heights", sum_fun,
                      "-burninTrees", burnin_trees,
                      trees_file, '/dev/stdout')

    result_tree_nexus <- suppressWarnings(system2('treeannotator',
                                                  command_args,
                                                  stdout = TRUE,
                                                  stderr = FALSE))

    if(is.null(attr(result_tree_nexus,"status"))){ #The exit code is only kept when it is not 0
      result_tree_connection <- textConnection(result_tree_nexus)
      result_tree <- ape::read.nexus(result_tree_connection)
    }

  } else {
    dir.create(dirname(outname),showWarnings = FALSE,recursive = TRUE)  #Making output directory if necessary
    command_args <- c(treeannotator_args,
                      "-heights", sum_fun,
                      "-burninTrees", burnin_trees,
                      trees_file, outname)
    treeannotator_exit_code <- system2('treeannotator',
                                       command_args,
                                       stdout = FALSE,
                                       stderr = FALSE)
    #Read MCC tree and return it
    if (treeannotator_exit_code == 0) {
      result_tree <- ape::read.nexus(outname)
    }
  }

  if (is.null(result_tree))
    stop("Problem generating the summary tree estimate")

  return(result_tree)
}

#' Uses BEAST's treeannotator to summarize the posterior distribution of trees
#' from memory
#'
#' @details [get_summary_tree] wrapper that deals with temp files since
#'  treeannotator cannot read from stdin. If this function is used in parallel,
#'  different threads should use different temp_names to avoid collisions.
#'
#' @param trees trees files to summarize
#' @param temp_name name of the temporary tree file. If this function is used in
#'  parallel, each thread should use a different temp_name to avoid issues
#' @inheritParams get_summary_tree
#'
#' @return the summary tree as an ape::phylo
#' @export

get_summary_tree_from_trees <- function(trees,
                                        temp_name = "temptree",
                                        burnin_p = 0.1,
                                        burnin_trees = NULL,
                                        outname = NULL,
                                        sum_fun = c("mean","median","keep","ca"),
                                        treeannotator_args = NULL){
  #Write trees in tmp
  temp_tree_file <- tempfile(temp_name)
  ape::write.nexus(trees,file = temp_tree_file)

  #Calculate the burnin states since we have the trees in memory to avoid reading the file later again
  if(is.null(burnin_trees))
    burnin_trees <- ceiling(length(trees) * burnin_p)

  #Summarize the tree
  return_tree <- get_summary_tree(trees_file = temp_tree_file,
                                  burnin_trees = burnin_trees,
                                  outname = outname,
                                  sum_fun = sum_fun,
                                  treeannotator_args = treeannotator_args)

  #Delete temporary tree file
  unlink(temp_tree_file)

  return(return_tree)
}

#' Marginalize posteriors over S across the whole experiment
#'
#' This function uses Bayesian Model Averaging to marginalize over S. It can
#' generate posterior samples, continuous parameter estimates, and point tree
#' estimates.
#'
#' @details All .trees files within the objective directory (and
#'   sub-directories) will be analyzed. Files with the same name will be
#'   considered independent phyfum runs under the same conditions and will be
#'   merged after burning before re-sampling. If a mle_file is provided and
#'   certain conditions do not have a probability in the MLE file but they do
#'   have data, the data is not considered and a warning is generated.
#'
#' @param file_dir directory that contains .trees and .log files (with or
#'   without subdirectories).
#' @param out_dir directory that will contain a summary statistics table per
#'   patient. NULL to disable summary statistics.
#' @param out_tree_dir directory that will contain a summary tree per patient.
#'   NULL to disable writing the summary tree (it will still be calculated
#'   unless summarize_trees = FALSE).
#' @param out_trace_dir directory that will contain the resampled trace (.trees
#'   and .log files) marginalized over S. NULL to disable trace storage.
#' @param mle_method sorted list of preferred methods to use the estimated
#'   marginal likelihood for model selection. Choose between PS (Path Sampling),
#'   SS (Stepping Stone), and HME (Harmonic Mean Estimator). Only useful when
#'   mle_file is provided.
#' @inheritParams calculate_model_averaging_posteriors_file
#' @inheritParams model_average_traces
#' @inheritParams resample_posterior
#' @inheritParams mcmc_qc_condition
#' @param summarize_trees boolean indicating if trees must be summarized (e.g.,
#'   calculating the MCC tree) or not.
#' @param tree_branch_sum_fun summarization function for summary tree branches
#'   (choose from mean, median, keep, and ca). See [get_summary_tree]
#' @param summary_tree_suffix suffix to the patient ID to be used as filename
#'   when storing a summary tree
#' @param summary_params_suffix suffix to the patient ID to be used as filename
#'   when storing a summary of continuous parameters
#' @param all_summary_params_suffix suffix to the outdir to be used as filename when
#'   storing a summary of all continuous parameters across patients
#' @param n_cores number of cores to use for coarse-grained parallelism (per S
#'   value per patient)
#' @param precision maximum difference between two floating point numbers to be
#'   considered equal (used to check posterior probabilities)
#' @inheritDotParams get_summary_tree_from_trees treeannotator_args
#' @inheritDotParams get_summary_tree treeannotator_args
#'
#' @returns List with data.tables with all = convergence statistics, problematic
#'   = parameter with problems, MLE = mle estimation table data.table. Each has
#'   a column indicating the condition (i.e., filename of the input file without
#'   .trees)
#' @export

#TODO save a file with the pp_table
model_average <- function(file_dir,
                          out_dir,
                          out_tree_dir,
                          out_trace_dir,
                          mle_file = NULL,
                          mle_method = c("PS","SS","HME"),
                          log_prior_probabilities = NULL,
                          n_marginalized_posterior_samples = 10000,
                          max_expected_samples = 1,
                          warn_biased_marginalization = TRUE,
                          burnin_p = 0.1,
                          cred_mass = 0.95,
                          summarize_trees = TRUE,
                          tree_branch_sum_fun = c("mean", "median", "keep", "ca"),
                          n_cells_regex = ".*\\.([0-9]+)cells\\..*",
                          basename_regex = "\\.[0-9]+cells\\..*",
                          summary_tree_suffix = ".tree",
                          summary_params_suffix = ".csv",
                          all_summary_params_suffix = "all.csv",
                          #n_cells_regex = ".*s([0-9]+).*",
                          #basename_regex = "_s[0-9]+.*",
                          n_cores = 1,
                          precision = NULL,
                          ...){
  #Argument check and logging
  if(is.null(precision))
    precision <- .phyfumr_env[['precision']]

  args <- c(as.list(environment()), list(...))

  #If we have input mle data, we calculate the posterior probabilities, otherwise, it will be done during re-sampling
  if (!is.null(mle_file)) {
    pp_table <- calculate_model_averaging_posteriors_file(mle_file,
                                                          n_cells_regex = n_cells_regex,
                                                          basename_regex = basename_regex)
    invisible(check_model_averaging_posteriors(pp_table,
                                               n_posterior_samples = n_marginalized_posterior_samples,
                                               max_expected_samples = max_expected_samples,
                                               warn_biased_marginalization=T))
  }

  #Prepare for outputs
  for (this_dir in c(out_dir, out_trace_dir, out_tree_dir)){
    if(!is.null(this_dir))
      dir.create(this_dir,showWarnings = FALSE,recursive = TRUE)
  }

  #Warn if some of the desired outputs can't be generated
  if(!is.null(out_tree_dir) && summarize_trees == FALSE){
    warning("Summary tree output directory provided but flag to calculate summary trees was disabled. Re-enabled here.")
    summarize_trees <- TRUE
  }
  if(summarize_trees && .phyfumr_env[['treeannotator']] == FALSE){
    warning("Summary trees will not be generated since treeannotator is not in your PATH environment variable")
    out_tree_dir <- NULL
  }

  #Getting dataset information
  trees_files <- list.files(path = file_dir, pattern = "*.trees$",
                            full.names = TRUE,
                            recursive = TRUE)

  trees_files <- trees_files[file.size(trees_files)!=0] #Do not consider empty traces

  run_names <- basename(trees_files)
  runs_info <- data.table::data.table(trees_files=trees_files,run_names = gsub(pattern = ".trees",
                                                                               replacement = "",
                                                                               x = run_names))
  runs_info[,`:=`(patient=gsub(basename_regex,"",run_names))]

  # Main patient loop, returning a list with a stats data.table and summary tree per patient
  result_list <- lapply(runs_info[,unique(patient)],function(this_patient){
    #Parses phyfum output and merges them if multiple traces per condition had been used
    print(sprintf("Patient %s...",this_patient))
    cat("\tParsing data\n")
    this_patient_traces <- slapply(
      lapply(
        split(runs_info[patient==this_patient,],by="run_names"),
        '[[',"trees_files"),
      function(trees_files,burnin_p){
        load_condition_traces(trees_files = trees_files,burnin_p = burnin_p)
      },
      burnin_p = burnin_p, #burnin applied here
      cl = n_cores)
    cat("\tDone\n")

    #Performs the model-averaging, using provided MLEs or using HME to calculate them (yes we know HME is evil but in our simulations it works as well as PS or SS for S)
    cat("\tModel averaging...\n")
    if(is.null(mle_file)) { #Calculating MLEs here
      this_averaged_chain <- model_average_traces(this_patient_traces,
                                                  log_joint_probabilities = NULL,
                                                  log_prior_probabilities = NULL,
                                                  posterior_probabilities = NULL,
                                                  burnin_p = 0,
                                                  n_marginalized_posterior_samples = n_marginalized_posterior_samples,
                                                  max_expected_samples = max_expected_samples,
                                                  warn_biased_marginalization = warn_biased_marginalization)
    } else { #Using precalculated MLEs
      these_s_with_data <- sprintf("S_%02d",as.numeric(gsub(n_cells_regex,"\\1",names(this_patient_traces))))
      these_s_with_MLEs <- sapply(pp_table[patient==this_patient,.SD,.SDcols=grep("S",colnames(pp_table),value = TRUE)],function(x)!is.na(x))
      these_s_with_MLEs <- names(these_s_with_MLEs[these_s_with_MLEs])
      valid_traces <- these_s_with_data %in% these_s_with_MLEs

      if(length(these_s_with_MLEs) != length(these_s_with_data)){
        if(sum(valid_traces) < length(these_s_with_MLEs)){#Missing data for existing MLEs
          stop("ERROR: Missing trace data for existing MLEs")
        } else {
          warning("Existing data missing MLEs in provided MLE file. MLE file may be outdated or you are intentionally filtering results")
        }
      }

      these_posterior_probabilities <- as.numeric(pp_table[patient==this_patient,.SD,.SDcols=these_s_with_data[valid_traces]])
      if(abs(sum(these_posterior_probabilities,na.rm = TRUE) - 1) > precision)
        stop("ERROR: Posterior probabilities of S do not add up to 1. Most probably, there are missing Phyfum traces for S values contained in the MLE file")
      this_averaged_chain <- model_average_traces(this_patient_traces[valid_traces],
                                                  posterior_probabilities = these_posterior_probabilities,
                                                  burnin_p = 0,
                                                  warn_biased_marginalization = FALSE) #Warnings already generated in check_model_averaging_posteriors if desired
    }
    rm(this_patient_traces) #Making it easy with the memory (I assume the GC would do it but this is obvious and takes a lot of memory)

    #Write the averaged chain
    if(!is.null(out_trace_dir)) {

      this_out_trees_file <- paste(sep="/",out_trace_dir,paste0(this_patient,".trees"))
      this_out_log_file <- gsub(".trees$",".log$",this_out_trees_file)

      #Trees
      ape::write.nexus(this_averaged_chain$trees, file = this_out_trees_file)

      #Continuous params
      writeLines(c("# Marginalizing over S after Phyfum",paste0("# ",lubridate::now()),paste0(collapse=",",c("# arguments:",args))),this_out_log_file)
      suppressWarnings(utils::write.table(this_averaged_chain$ptable, file = this_out_log_file,quote = FALSE,row.names = FALSE,sep = "\t",append = TRUE))
    }

    #Summarize continuous parameters
    #Selecting params to analyze
    params <- names(this_averaged_chain$ptable)
    params <- params[!params %in% .phyfumr_env[["not_params"]]]
    params <- params[!params %in% detect_constants(this_averaged_chain$ptable,params)]

    #Calculating summary stats
    infun <- function(param,thedata,cred_mass){
      thisHDI <- HDInterval::hdi(thedata[,get(param)],credMass = cred_mass)
      data.table::data.table("param" = param,
                             "medianP"= stats::median(thedata[,get(param)]),
                             "meanP"= mean(thedata[,get(param)]),
                             "HDILower" = thisHDI[1],
                             "HDIUpper" = thisHDI[2])
    }

    this_param_table <- data.table::rbindlist(lapply(params,infun,thedata = this_averaged_chain$ptable,cred_mass = cred_mass))

    #Output (if needed)
    check_write_csv(table = this_param_table,
                    flag = out_dir,
                    outdir = out_dir,
                    filename = paste0(this_patient,"",summary_params_suffix))

    #Summarize the trees
    if(summarize_trees) {
      if(!is.null(out_tree_dir)){ #Whether we want to save the summary tree or not
        this_tree_outname <- paste(sep="/",out_tree_dir,paste0(this_patient,summary_tree_suffix))
      } else {
        this_tree_outname <- NULL
      }
      if(!is.null(out_trace_dir)){ #The trees file was written earlier so we just re-use it here
        this_summary_tree = get_summary_tree(trees_file = this_out_trees_file,
                                         burnin_trees = 0,
                                         outname = this_tree_outname,
                                         sum_fun = tree_branch_sum_fun,
                                         ...)
      } else { #Without storing the trees file
        this_summary_tree = get_summary_tree_from_trees(trees = this_averaged_chain$trees,
                                         burnin_trees = 0,
                                         outname = this_tree_outname,
                                         sum_fun = tree_branch_sum_fun,
                                         ...)
      }
    }

    cat("\tDone\nDone\n")
    return(list(ptable = this_param_table,
                trees = this_summary_tree))
  })

  names(result_list) <- runs_info[,unique(patient)]
  merged_table <- data.table::rbindlist(lapply(result_list,'[[',"ptable"),idcol = "patient")
  all_summary_trees <- lapply(result_list,'[[',"trees")
  class(all_summary_trees) <- "multiPhylo"

  check_write_csv(table = merged_table,
                  flag = all_summary_params_suffix,
                  outdir = out_dir,
                  filename = all_summary_params_suffix)


  return(list(ptable = merged_table,
              trees = all_summary_trees))
}
