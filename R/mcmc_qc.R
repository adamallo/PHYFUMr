# R/mcmc_qc.R

#Stats
######

#' Method-of-moments estimator of AICM
#'
#' Posterior simulation-based analogue of Akaike's information criterion (AIC)
#' through Markov chain Monte Carlo Raftery et al. 2007
#' @details Translated from BEAST's logMarginalLikelihoodAICM function
#'
#' @param loglikelihoods vector of posterior sample of loglikelihoods
#'
#' @returns AICM
#' @export

AICM <- function(loglikelihoods){
  return(2*stats::var(loglikelihoods)-2*mean(loglikelihoods))
}

#Utils
######

#' Detects constant parameters
#'
#' @details Defined as those parameters with sd = 0
#'
#' @param thedata data.table with a column per paramter and row per posterior
#'   sample
#' @param params list of columns to consider
#'
#' @returns vector of names of parameters with sd == 0

detect_constants <- function(thedata,params){
  params[t(thedata[,lapply(.SD,FUN = function(x){stats::sd(x)==0}),.SDcols = params])]}


#' Fixes NEXUS files of unfinished PHYFUM runs
#'
#' @param file .trees files
#' @param backup Logical. Whether a .trees_bkp with the original contents should
#'   be created before fixing the file
#'
#' @returns TRUE if the file has been fixed, FALSE if it didn't need fixing

fix_nexus <- function(file,backup = TRUE) {
  last_line <- ifelse(.phyfumr_env[["system_tail"]],
                      system2("tail", c("-n", "1", shQuote(normalizePath(file))), stdout = TRUE),
                      utils::tail(readLines(file, warn = FALSE), 1))

  if(!grepl("END;|ENDBLOCK;", last_line, ignore.case = TRUE)){
    if(backup==TRUE){
      if(file.copy(file, paste0(file,"_bkp"), overwrite = TRUE)==FALSE)
        stop("ERROR: cannot write the backup Nexus file")# overwrite=FALSE by default
    }
    cat("End;",file=file,append=T)
    return(TRUE)
  }
  return(FALSE)
}

#RWTY integration
#################

#Parsing input traces
#######################################

#' Ensures that chains are in the expected format
#'
#' @param chains MCMC data. List of [rwty::load.trees()] values
#'
#' @returns List of [rwty::load.trees()] values
#' @keywords internal

ensure_list_chains <- function(chains){
  if(methods::is(chains,"rwty.chain"))
    chains <- list(chains)
  if(is.null(names(chains)))
    data.table::setattr(chains,"name",seq(1,length(chains)))
  return(chains)
}

#' Removes burnin states in RWTY chain lists
#'
#' @inheritParams ensure_list_chains
#' @param burnin_p Proportion of MCMC samples to discard as burnin
#'
#' @returns list of RWTY chains with burnin samples removed
#' @keywords internal

remove_burnin_chains <- function(chains,burnin_p) {
  chains <- ensure_list_chains(chains)
  return(lapply(chains,function(chain){
    last_sample <- length(chain$trees)
    first_sample <- ceiling(burnin_p*last_sample + .phyfumr_env[['precision']]) #min 1
    chain$trees <- chain$trees[seq.int(from=first_sample,to=last_sample)]
    chain$ptable <- chain$ptable[seq.int(from=first_sample,to=last_sample),]
    chain
    }))
}

#' Load phyfum posterior samples
#'
#' @param trees_files array of phyfum's output .trees files (one per independent
#'   run)
#' @param log_files array of phyfum's output .log files (one per independent
#'   run). If not used, this function will assume they have the same name and
#'   location as their .trees counterparts
#' @inheritParams remove_burnin_chains
#' @inheritDotParams fix_nexus backup
#'
#' @returns list of RWTY's chains
#' @keywords internal

load_chains <- function(trees_files, log_files=NULL, burnin_p=0, ...){
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

  message("Loading traces...")
  sink(nullfile()) #shut up!
  chains <- invisible(lapply(seq_along(trees_files),FUN=function(i_file){
    fix_nexus(trees_files[i_file],...)
    chain <- rwty::load.trees(trees_files[i_file],format="BEAST",logfile = log_files[i_file]);
    chain}))
  sink()
  message("Done\n")

  if(!is.null(burnin_p)&&burnin_p!=0){
    return(remove_burnin_chains(chains,burnin_p))
  } else {
    return(chains)
  }
}

#' Obtains the parameter trace from from multiple RWTY MCMC samples
#'
#' @inheritParams ensure_list_chains
#' @param burnin_p Proportion of MCMC samples to discard as burnin. Calculated
#'   by trace and flagged or removed
#' @param remove_burnin Logical that indicates if the burnin states should be
#'   removed (True) or flagged (False, Default)
#'
#' @returns data.table with same columns as chain$ptable + a chain factor
#' @keywords internal

extract_ptable <- function(chains,burnin_p,remove_burnin=FALSE){
  chains <- ensure_list_chains(chains)

  result_table <- data.table::rbindlist(lapply(chains,function(chain){
    this_table <- data.table::data.table(chain$ptable)
    this_table[,`:=`(burnin=ifelse(state<max(state)*burnin_p,TRUE,FALSE))]
    }),idcol = "chain")
  result_table[,`:=`(chain=as.factor(chain))]

  if(remove_burnin){
    return(result_table[burnin==FALSE,])
  } else {
    return(result_table)
  }
}

#' Updates RWTY parameter tables for multiple MCMC chains
#'
#' Data can be extracted using [extract_ptable()], modified, and then this
#' function will add it back to the chains
#'
#' @inheritParams ensure_list_chains
#' @param thedata data.table with modified chain data
#' @param exclude_cols columns in thedata that should not be added to chains.
#'   Defaults to those columns added by [extract_ptable()] except burnin
#'
#' @returns new chains
#' @keywords internal

update_ptable <- function(chains,thedata,exclude_cols=c("chain")){
  the_cols <- names(thedata)
  the_cols <- the_cols[!the_cols %in% exclude_cols]
  for(ichain in 1:length(chains)){
    chains[[ichain]]$ptable <- data.table::data.table(thedata[chain==ichain,.SD,.SDcols = the_cols])
    if(nrow(chains[[ichain]]$ptable) != length(chains[[ichain]]$trees))
      stop("ERROR: ptable and trees are not of the same length")
  }
  return(chains)
}

#' Updates RWTY chains removing burnin states in parameter and tree objects for
#' multiple MCMC chains
#'
#' Data can be extracted and burnin tagged using [extract_ptable()], modified
#' manually if desired, and added back using [update_ptable()]. This function,
#' removes burnin states marked in a column called burnin (T: removed, F: kept)
#'
#' @inheritParams ensure_list_chains
#'
#' @returns new chains
#' @keywords internal

remove_marked_burnin <- function(chains) {
  for(ichain in 1:length(chains)){
    if("burnin" %in% colnames(chains[[ichain]]$ptable)){
      is_burnin <- chains[[ichain]]$ptable$burnin
      chains[[ichain]]$ptable <- chains[[ichain]]$ptable[burnin==F,]
      chains[[ichain]]$ptable[,`:=`(burnin=NULL)]
      chains[[ichain]]$trees <- chains[[ichain]]$trees[!is_burnin]
    } else {
      stop("ERROR: ptable does not contain the burnin column")
    }
  }
  return(chains)
}

#' Merges a list of RWTY chains into one
#'
#' @inheritParams extract_ptable
#'
#' @details It assumes samples are taken at regular intervals to perform the
#'   burnin
#'
#' @returns rwty.chain with the chains concatenated and the burnin states (trees
#'   and parameters) removed.
#' @keywords internal

merge_traces <- function(chains,burnin_p=0){
  returnsChain <- chains[[1]]
  returnsChain$trees <- do.call(c,lapply(chains,function(chain){
    last_sample <- length(chain$trees)
    first_sample <- ceiling(burnin_p*last_sample + .phyfumr_env[['precision']]) #min 1
    chain$trees[seq.int(from=first_sample,to=last_sample)]
  }))
  returnsChain$ptable <- do.call(rbind,lapply(chains,function(chain){
    last_sample <- nrow(chain$ptable)
    first_sample <- ceiling(burnin_p*last_sample + .phyfumr_env[['precision']]) #min 1
    chain$ptable[seq.int(from=first_sample,to=last_sample),]
  }))
  return(returnsChain)
}

#' Generates a trace of tree distances from a focal tree
#'
#' To calculate ESS using the pseudo ESS approach. Typically using
#' [tree_convergence_stats()]
#'
#' @param tree_list list of trees
#' @param treedist Tree distance to be used to calculate the distance to the
#'   focal tree ("PD" or "RF")
#' @param tree_i pre-specified element of tree_list that should be used as
#'   reference tree. If NULL (default) it is sampled at random from tree_dist
#' @param focal_tree focal tree to use instead of sampling it from the chain
#'
#' @returns array of tree distances
#' @keywords internal
#' @seealso [generate_tree_traces()] [tree_convergence_stats()]

generate_tree_trace_reftree <- function (tree_list,
                                         treedist = c("PD","RF"),
                                         tree_i = NULL,
                                         focal_tree = NULL) {
  treedist <- match.arg(treedist)
  if(!is.null(focal_tree)){
    distances <- rwty:::tree.distances(tree_list, treedist = treedist, focal.tree = focal_tree)
  } else {
    if (is.null(tree_i)){
      tree_i <- sample(1:length(tree_list), 1)
    }
    distances <- rwty:::tree.distances(tree_list, tree_i, treedist = treedist)
  }
  return(distances)
}

#' Generates a tree-topoloy trace
#'
#' To calculate ESS using the pseudo ESS approach. Typically using
#' [tree_convergence_stats()].
#'
#' @details Tree_i corresponds to the tree number after burnin if burnin_p>0
#'
#' @inheritParams extract_ptable
#' @inheritParams generate_tree_trace_reftree
#' @inheritDotParams generate_tree_trace_reftree tree_i focal_tree
#'
#' @returns data.table (long format) with a row per MCMC sample and columns for
#'   treedistance and chain
#' @keywords internal
#' @seealso [generate_tree_trace_reftree()] [tree_convergence_stats()]

generate_tree_traces <- function(chains,burnin_p = 0, treedist = "PD", remove_burnin = TRUE,...) {
  chains <- ensure_list_chains(chains)

  ellipsis_args <- list(...)

  result_table <- data.table::rbindlist(lapply(chains,function(chain,burnin_p,remove_burnin){
    last_sample <- length(chain$trees)
    first_sample <- ceiling(burnin_p*last_sample + .phyfumr_env[['precision']])

    if(burnin_p>0 && isTRUE(remove_burnin)){
        if(!is.null(ellipsis_args[["tree_i"]]))
          warning("When using a fixed tree_i in generate_tree_traces with burnin_p > 0 and remove_burnin=T, tree_i corresponds to the tree number after burnin has been removed")

      this_data <- data.table::data.table(treedistance=generate_tree_trace_reftree(chain$trees[seq.int(from=first_sample,to=last_sample)],treedist = treedist,...)$topological.distance, treedist = treedist)
      this_data[,`:=`(burnin=FALSE)]
    } else {
      this_data <- data.table::data.table(treedistance=generate_tree_trace_reftree(chain$trees,treedist = treedist,...)$topological.distance, treedist = treedist)
      this_data[,`:=`(burnin=TRUE)]
      this_data[first_sample:last_sample,`:=`(burnin=FALSE)]
      }

    return(this_data)},burnin_p = burnin_p, remove_burnin = remove_burnin),
    idcol = "chain")

  if(isTRUE(remove_burnin)){
    return(result_table[burnin==FALSE,])
  } else {
    return(result_table)
  }

}

#Convergence
############

#' Generates a table with convergence statistics
#'
#' @param stan_table MCMC data for one parameter, with as many columns as
#'   independent chains of same length
#' @param cred_mass mass of the selected high posterior density interval (HDI)
#'
#' @returns data.table with a row per parameter and colums for Rhat, ESS-bulk,
#'   and ESS-tail, median, mean, and HDI interval (HDILower and HDIUpper)
#' @seealso [rstan::Rhat()] [rstan::ess_bulk()] [rstan::ess_tail()]
#'   [HDInterval::hdi()]
#' @export

convergence_stats <- function(stan_table,cred_mass = 0.95) {
  stan_matrix <- as.matrix(stan_table)
  thisHDI <- HDInterval::hdi(stan_matrix,credMass = cred_mass)
  data.table::data.table("Rhat"=rstan::Rhat(stan_matrix),
             "essB"=rstan::ess_bulk(stan_matrix),
             "essT"=rstan::ess_tail(stan_matrix),
             "medianP"=stats::median(stan_matrix),
             "meanP"=mean(stan_matrix),
             "HDILower"=thisHDI[1],
             "HDIUpper"=thisHDI[2])
}

#' Generates a table with tree convergence statistics using the pseudo ESS
#' approach
#'
#' @details The pseudo-ESS approach generates a trace by calculating a tree
#'   distance to a focal tree (usually chosen at random). Then, that trace can
#'   be used to assess mixing, and anything else typically done with continuous
#'   parameters. To improve its accuracy, a number of replicates are generated
#'   and mean convergence statistics reported.
#'
#' @inheritParams ensure_list_chains
#' @inheritParams extract_ptable
#' @inheritParams generate_tree_trace_reftree
#' @inheritParams convergence_stats
#' @param n_focal_trees Number of focal trees to calculate tree distances
#' @param by_chain Calculates Rhat, and the two ESS parameters both using all
#'   chains and by chain (otherwise, using all of them only)
#' @param cl cl settings for [pbapply::pblapply()] if available
#'
#' @inherit convergence_stats return
#' @seealso [convergence_stats()]
#' @export

tree_convergence_stats <- function(chains,
                                   burnin_p = 0,
                                   n_focal_trees = 20,
                                   treedist = c("PD","RF"),
                                   cred_mass = 0.95,
                                   by_chain = T,
                                   cl = 1){

  chains <- ensure_list_chains(chains)
  nSamples <- sapply(chains,function(chain) length(chain$trees))
  if(any(nSamples * (1-burnin_p) < n_focal_trees)){
    stop("Not enough samples left. Revise your burnin_p parameter")
  }
  treedist <- match.arg(treedist)

  merged_chain <- merge_traces(chains,burnin_p = burnin_p)
  focal_trees <- sample(merged_chain$trees,size = n_focal_trees)

  # infun <- function(tree_i,chain,treedist,...){
  #   convergence_stats(generate_tree_traces(chains,treedist = treedist, tree_i = tree_i,...)[,.(treedistance)],cred_mass = cred_mass)
  # }

  infun <- function(i,focal_trees,chains,treedist,...){
    convergence_stats(data.table::dcast(generate_tree_traces(chains = chains,treedist = treedist, focal_tree = focal_trees[i], remove_burnin = T,...)[,idx := seq_len(.N), by=chain],idx~chain,value.var = "treedistance")[,-1][,c(list(.SD),list(complete=stats::complete.cases(.SD)))][complete==T,][,`:=`(complete=NULL)],cred_mass = cred_mass)
  }

  #WARNING: infun uses mclapply if rwrt.processors > 1. A thread explosion can easily happen here
  replicated_convergence_stats <- data.table::rbindlist(slapply(seq_along(focal_trees), infun,
                                                                focal_trees = focal_trees,
                                                                chains = chains,
                                                                burnin_p = burnin_p,
                                                                treedist = treedist,
                                                                cl = cl),idcol = "focal_tree_i")

  result_table <- replicated_convergence_stats[,lapply(.SD,mean,na.rm=T),.SDcols=names(replicated_convergence_stats)[!names(replicated_convergence_stats)%in%c("focal_tree_i")]]
  result_table[,`:=`(chain = "ALL")]

  if(isTRUE(by_chain) && length(chains)>1) {
    result_table_by_chain <- data.table::rbindlist(suppressWarnings(lapply(chains,FUN = function(chain){
      focal_trees <- sample(chain$trees[seq.int(from=ceiling(length(chain$trees)*(burnin_p) + .phyfumr_env[['precision']]),to=length(chain$trees))],size = n_focal_trees)

      replicated_convergence_stats <- data.table::rbindlist(slapply(seq_along(focal_trees), infun,
                                                                    focal_trees = focal_trees,
                                                                    chains = ensure_list_chains(chain),
                                                                    burnin_p = burnin_p,
                                                                    treedist = treedist,
                                                                    cl = cl),idcol = "focal_tree_i")

      this_chain_stats <- replicated_convergence_stats[,lapply(.SD,mean,na.rm=T),.SDcols=names(replicated_convergence_stats)[!names(replicated_convergence_stats)%in%c("focal_tree_i")]]
      return(this_chain_stats)
    })),idcol="chain")
    result_table <- cbind(result_table,result_table_by_chain[,.(medianPCV=stats::sd(medianP)/mean(medianP),meanPCV=stats::sd(meanP)/mean(meanP))])
    result_table <- rbind(result_table,result_table_by_chain,fill=T)
  }

  result_table[,`:=`(param = paste(sep="_",treedist,"treeTopology"))]
  data.table::setkey(result_table,param)

  return(result_table)
}

#' Generates a table with convergence statistics for all continuous parameters
#'
#' @inheritParams tree_convergence_stats
#' @param thedata MCMC posterior samples in data.table format after
#'   [extract_ptable()]
#' @param params List of parameter names to analyze. NULL if all are wanted
#'   (default). In that case, known phyfum output columns that are not
#'   parameters (chain, state), and constant parameters as in
#'   [detect_constants()] are removed
#'
#' @return data.table with a row per parameter and columns with convergence
#'   statistics from [convergence_stats()]
#' @seealso [extract_ptable()]

param_convergence_stats <- function(thedata,
                                 params = NULL,
                                 cred_mass = 0.95,
                                 by_chain = TRUE,
                                 cl=1){

  if(is.null(params)) {
     params <- names(thedata)
     params <- params[!params %in% .phyfumr_env[["not_params"]]]
     params <- params[!params %in% detect_constants(thedata,params)]
  }

  infun <- function(this_param,thedata) {
    convergence_stats(data.table::dcast(thedata[order(chain,state),][,idx := seq_len(.N), by=chain],idx~chain,value.var = this_param)[,-1][,c(list(.SD),list(complete=stats::complete.cases(.SD)))][complete==T,][,`:=`(complete=NULL)], cred_mass = cred_mass)[,`:=`(param = this_param)]
  }

  result_table <- data.table::rbindlist(slapply(params,infun,thedata = thedata,cl=cl))
  data.table::setkey(result_table, param)
  result_table[,`:=`(chain = "ALL")]

  if(isTRUE(by_chain) && length(thedata[,unique(chain)]) > 1){
    result_table_by_chain <- data.table::rbindlist(lapply(thedata[,unique(chain)],FUN = function(this_chain){
      this_chain_result_table <- data.table::rbindlist(slapply(params,infun,thedata = thedata[chain==this_chain,],cl=cl))
      this_chain_result_table[,`:=`(chain = this_chain)]
      return(this_chain_result_table)
    }))
    data.table::setkey(result_table_by_chain,param)
    result_table <- result_table[result_table_by_chain[,.(medianPCV=stats::sd(medianP)/mean(medianP),meanPCV=stats::sd(meanP)/mean(meanP)),by=param]]
    result_table <- data.table::rbindlist(list(result_table,result_table_by_chain),fill = TRUE)
  }
  return(result_table)
}

#Plotting: Continuous parameters
################################

#' Summary function for point estimates
#'
#' Used for plotting
#'
#' @param x data
#' @param type character indicating if we want to use the median, mean, or mode
#' @param ... arguments for the corresponding summary function
#'
#' @returns point estimate
#' @seealso [plot_continuous_parameters()]

sumfun <- function(x,type=c("mean","median","mode"),...){
  type=match.arg(type)
  if(type=="mean"){
    return(mean(x,...))
  } else if (type=="median"){
    return(stats::median(x,...))
  } else if (type=="mode") {
    return(mode(x,...))
  }
}

#' Plots posterior sample and estimated values for a parameter
#'
#' @param thedata data.table, with burnin removed if needed
#' @param params array of parameter names
#' @param ndigits number of digits to round numbers printed in the plot
#'
#' @returns list of ggplot2 plots
#' @keywords internal

plot_continuous_parameters <- function(thedata,params,ndigits=3){
  lapply(params,function(param){
    ggplot2::ggplot(data=thedata,ggplot2::aes(x=.data[[param]],group=.data$chain,fill=.data$chain,color=.data$chain)) +
      ggplot2::geom_density(alpha=0.3,linewidth=1) +
      ggplot2::stat_summary(ggplot2::aes(xintercept = ggplot2::after_stat(x),y=0,linetype="Estimated",fill=NULL),fun=sumfun,geom="vline",orientation="y", show.legend = FALSE) +
      ggplot2::stat_summary(ggplot2::aes(xintercept = ggplot2::after_stat(x),y=0,linetype="Estimated",fill=NULL,group=NULL),fun=sumfun,geom="vline",orientation="y",color="black") +
      ggplot2::stat_summary(fun=sumfun,geom=ggrepel::GeomTextRepel,ggplot2::aes(x=.data[[param]],label=round(ggplot2::after_stat(x),digits=ndigits),y=0),orientation="y",hjust=-0.5,direction="y", show.legend = FALSE) +
      ggplot2::stat_summary(fun=sumfun,geom=ggrepel::GeomTextRepel,ggplot2::aes(x=.data[[param]],label=round(ggplot2::after_stat(x),digits=ndigits),y=0),orientation="y",hjust=1.5,inherit.aes=FALSE, show.legend = FALSE) +
      ggplot2::scale_linetype(name="") +
      ggplot2::scale_color_brewer(name="Chain",type = "qual",palette = 2) +
      ggplot2::scale_fill_brewer(name="Chain", type = "qual",palette = 2) +
      ggplot2::scale_x_continuous(name="Parameter value") +
      ggplot2::scale_y_continuous(name="Density") +
      ggplot2::labs(title=param) +
      cowplot::theme_cowplot()})
}

#' Generates a pdf file with plots for a list of continuous parameters
#'
#' @inheritParams plot_continuous_parameters
#' @param out_dir output directory
#' @param out_name output name
#' @param width width in inches
#' @param height height in inches
#'
#' @returns NULL
#' @export

write_continuous_parameter_plots <- function(thedata,params,out_dir,out_name,ndigits=3,width=7,height=7){
  test <- plot_continuous_parameters(thedata = thedata,
                                      params = params,
                                      ndigits = ndigits)
  grDevices::pdf(file=paste0(out_dir,"/",out_name),width = width, height = height)
  lapply(plot_continuous_parameters(thedata = thedata,
                             params = params,
                             ndigits = ndigits),print)
  grDevices::dev.off()
  return(NULL)
}

#Integrative functions (i.e. pipelines)
#######################################

#' Complete MCMC QC per condition (patient and S parameter value)
#'
#' This function generates plots and tables to aid assessing phyfum's mixing and
#' convergence. It can use several chains per run, and automatically flags some
#' issues using Rhat and ESS.
#'
#' @details When n_cores is > 1, RNGkind("L'Ecuyer-CMRG") random number
#'   generator is used instead of the default. This allows mclapply to be
#'   reproducible.
#'
#' @inheritParams load_chains
#' @param out_dir directory for tabular outputs
#' @param plot_dir directory for graphical outputs
#' @param burnin_p proportion of input traces to be discarded as burnin
#' @param min_ess minimum ESS. Parameters with ESS under this value will be
#'   flagged as problematic
#' @param max_rhat maximum Rhat [rstan::Rhat()]. Parameters with Rhat higher
#'   than this will be flagged as problematic
#' @param cred_mass credibility mass to determine high density intervals
#' @param n_focal_trees number of focal trees to calculate a tree distance to
#'   estimate topological convergence using the pseudo-ESS method
#' @param treedist tree distance used for the topological pseudo-ESS
#' @param base_name common name to use as prefix for all outputs to share output
#'   folders between patients. By default, the name of the .trees files is used
#' @param problematic_output_suffix suffix for the filename of the tabular
#'   output with a list of parameters that may have mixing/convergence problems.
#'   If NULL, this output is not generated.
#' @param all_output_suffix suffix for the filename of the tabular output with
#'   all mixing and convergence statistics. If NULL, this output is not
#'   generated.
#' @param posterior_plots_suffix suffix for the marginal density plots file. If
#'   NULL, this output is not generated.
#' @param rwty_plots_suffix suffix for the RWTY plots. If NULL, this output is
#'   not generated.
#' @param correlation_plots_suffix suffix for the RWTY correlation plots. If
#'   NULL, this output is not generated.
#' @param mle_suffix suffix for the tabular output with marginal likelihood
#'   estimations using HME and AICm methods.
#' @param n_cores if NULL, using 1 core for rwty and let's data.table autodetect
#'   the number of cores (Default, good to run this function in series).
#'   Otherwise, it fixes the number of cores that rwty and data.table will use.
#'   They never run in parallel so both can use the maximum number of cores
#' @inheritDotParams fix_nexus backup
#' @returns list with data.tables with all = convergence statistics, problematic
#'   = parameter with problems, MLE = mle estimation table
#' @export

mcmc_qc_condition <- function(trees_files,
                            log_files=NULL,
                            out_dir,
                            plot_dir,
                            burnin_p=0.1,
                            min_ess=200,
                            max_rhat=1.1,
                            cred_mass=0.95,
                            n_focal_trees=20,
                            treedist=c("PD","RF"),
                            n_cores=NULL,
                            #maxCV=??,
                            base_name = NULL,
                            problematic_output_suffix = "problematicParams.csv",
                            all_output_suffix = "allParams.csv",
                            posterior_plots_suffix = "plotPosteriorDensities.pdf",
                            rwty_plots_suffix = "plotsRWTY.pdf",
                            correlation_plots_suffix = "plotCorrelations.jpeg",
                            mle_suffix = "MLE.csv",
                            ...
                            ){

  #arg checking and configuration based on it
  if(!is.null(n_cores) && n_cores == 1){
    data.table::setDTthreads(1) #Everything should be single thread to use coarse-grained parallelism
  } else {
    RNGkind("L'Ecuyer-CMRG") #Different random number generator, so that mclapply is reproducible
    data.table::setDTthreads(n_cores)
  }

  treedist <- match.arg(treedist)

  dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)
  dir.create(plot_dir,showWarnings = FALSE,recursive = TRUE)

  chains <- load_chains(trees_files,log_files)

  nSamples <- unique(sapply(chains,FUN=function(x){length(x$trees)}))

  if(length(nSamples) != 1){
    warning("Chains are not of the same length. This usually indicates phyfum runtime problems, but does not necessarily invalidate downstream analyses")
  }

  if(any(nSamples * (1-burnin_p) < n_focal_trees)){
    stop("Not enough samples left. Revise your burnin_p parameter")
  }
  #######################################

  #Data reorganization and augmentation
  #######################################
  these_data <- extract_ptable(chains,burnin_p = burnin_p,remove_burnin = FALSE)
  if(all(.phyfumr_env[["absolute_rates_params"]] %in% colnames(these_data)) && all(!.phyfumr_env[["relative_rates_params"]] %in% colnames(these_data))){
    these_data[,`:=`(clock.rate=flipflop.mu,flipflop.rlambda=flipflop.lambda/flipflop.mu,flipflop.MDbias=flipflop.mu/flipflop.gamma)]
  } else if (all(.phyfumr_env[["relative_rates_params"]] %in% colnames(these_data)) && all(!.phyfumr_env[["absolute_rates_params"]] %in% colnames(these_data))){
    these_data[,`:=`(flipflop.mu=clock.rate,flipflop.gamma=flipflop.mu/flipflop.MDbias,flipflop.lambda=flipflop.rlambda*clock.rate)]
  } else {
    warning("Parameterization not recognized as default phyfum with either absolute or relative rates. Parameter augmentation will not be performed.")
  }
  chains <- update_ptable(chains = chains,thedata = these_data)
  #######################################

  #Assess tree convergence
  #######################################
  message("Assessing topological convergence... ")
  topology_convergence <- tree_convergence_stats(chains = chains,
                                                 burnin_p = burnin_p,
                                                 n_focal_trees = n_focal_trees,
                                                 treedist = treedist,
                                                 cred_mass = cred_mass,
                                                 cl = n_cores,
                                                 by_chain = TRUE)

  #Make one tree-distance trace from the tree with index 1 for tree-trace plots
  tree_topology_trace <- generate_tree_traces(chains = chains,
                                              burnin_p = burnin_p,
                                              treedist = treedist,
                                              remove_burnin = FALSE,
                                              tree_i = 1)
  these_data <- cbind(these_data,tree_topology_trace[,.(treedistance)])
  message("Done\n")
  #######################################

  #Assessing convergence of continuous variables
  #######################################
  message("Assessing parameter convergence... ")
  continuous_convergence <- param_convergence_stats(thedata = these_data[burnin==FALSE,],
                          cl = n_cores)
  convergence <- rbind(continuous_convergence,topology_convergence)
  data.table::setkey(convergence,param)
  message("Done\n")
  #######################################

  #TODO add comparison between priors and posteriors and add it to the final check
  #TODO add distance to origin and entropy and add it to the final check

  #IO prep
  #######################################
  if(is.null(base_name)){
    base_name <- gsub(x = basename(trees_files[1]),pattern = ".trees",replacement = "")
  }

  #Plotting continuous parameters
  #######################################
  if(!is.null(posterior_plots_suffix)){
    message("Plotting continous parameters...")
    write_continuous_parameter_plots(thedata = these_data[burnin==FALSE,],
                                     params = names(these_data)[names(these_data) %in% continuous_convergence[,unique(param)]],
                                     out_dir = plot_dir,
                                     out_name = paste(sep="_",base_name,posterior_plots_suffix))
    message("Done\n")
  }
  #######################################

  #RWTY plots including traces
  #######################################
  if(!is.null(rwty_plots_suffix) || !is.null(correlation_plots_suffix)){
    message("Generating RWTY plots...")
    sink(nullfile())
    rwty_plots <- suppressWarnings(invisible(rwty::analyze.rwty(chains = remove_marked_burnin(chains),
                                     fill.color = 'posterior',
                                     params = colnames(chains[[1]]$ptable)[colnames(chains[[1]]$ptable) %in% continuous_convergence[,unique(param)]])))
    sink()
    message("Done\n")
  }

  #Print all but the correlations plots in pdf
  ######################################
  if(!is.null(rwty_plots_suffix)){
    message("Printing RWTY plots...")
    sink(nullfile())
    grDevices::pdf(file = paste0(plot_dir,"/",paste(sep="_",base_name,rwty_plots_suffix)), width = 10, height = 7)
    suppressWarnings(suppressMessages(print(rwty_plots[-grep(".correlations",names(rwty_plots))])))
    grDevices::dev.off()
    sink()
    message("Done\n")
  } else {
    warning("RWTY plots deactivated. If you need them, make sure to indicate a suffix using the rwty_plots_suffix argument")
  }

  #Correlations in jpeg to avoid computers crashing
  #######################################
  if(!is.null(correlation_plots_suffix)){
    message("Printing RWTY correlation plots...")
    sink(nullfile())

    for (plotname in names(rwty_plots[grep(".correlations",names(rwty_plots),value = TRUE)])){
      if (.phyfumr_env[["has_cairo"]]) {
        grDevices::jpeg(filename = paste0(plot_dir,"/",paste(sep="_",base_name,plotname,correlation_plots_suffix)),width = 1200, height = 1200,type = "cairo")
      } else if (.phyfumr_env[["sys"]] == "Darwin") {
        grDevices::jpeg(filename = paste0(plot_dir,"/",paste(sep="_",base_name,plotname,correlation_plots_suffix)),width = 1200, height = 1200,type = "quartz")
      } else {
        grDevices::jpeg(filename = paste0(plot_dir,"/",paste(sep="_",base_name,plotname,correlation_plots_suffix)),width = 1200, height = 1200)
      }
      suppressWarnings(suppressMessages(print(rwty_plots[[plotname]])))
      grDevices::dev.off()
    }
    sink()
    message("Done\n")
  } else {
    warning("Correlation plots deactivated. If you need them, make sure to indicate a suffix using the correlation_plots_suffix argument")
  }

  #Write table outputs
  #######################################

  #All convergence statistics
  check_write_csv(convergence,
                  all_output_suffix,
                  out_dir,
                  paste(sep="_",base_name,all_output_suffix),
                  "Mixing and convergence summary table deactivated. If you need it, make sure to indicate a suffix using the all_output_suffix argument")

  #Parameters that do not meet certain standards
  #For now, we consider that a param has a problem if:
  # 1. It has not achieved the minimum ess in all chains (including all combined)
  # 2. It has a Rhat > than the limit with all chains combined
  #TODO 3. The estimates from the different chains have a CV< ??? Unsure what this value should be
  problematic_params <- cbind(convergence[chain=="ALL",.(Rhat=is.na(Rhat)|Rhat>max_rhat),keyby=param],
                         convergence[,lapply(.SD,FUN = function(x){any(is.na(x)) | any(x<min_ess)}),keyby=param,.SDcols=c("essB","essT")][,-"param"])#,
                         #convergence[chain=="ALL",lapply(.SD,FUN = function(x){any(is.na(x)) | any(x>maxCV)}),keyby=param,.SDcols=c("medianPCV","meanPCV")][,-"param"])
  problematic_params <- problematic_params[,.(param,problem=any(.SD)),.SDcols=colnames(problematic_params)[!colnames(problematic_params)%in%c("param")],by=.I][problem==T,param]

  check_write_csv(convergence[param %in% problematic_params],
                  problematic_output_suffix,
                  out_dir,
                  paste(sep="_",base_name,problematic_output_suffix),
                  "Automatic parameter flagging output deactivated. If you need it, make sure to indicate a suffix using the problematic_output_suffix argument")

  #Analyze and communicate the number of problematic parameters
  #if only tree topology, warn that they should still be checked by eye and that it may be ok if the number of leaves is small. Show the number of leaves?
  #if there are more, warn that the results should be checked by hand and something in the analysis should be tweaked to improve mixing-convergence before proceeding
  if(length(problematic_params)==0){
    print("No insights of problematic parameters. Please, still review the plots to confirm the analysis worked properly before interpeting or publishing these results")
  } else if(length(grep("treeTopology",problematic_params,value = TRUE)) == length(problematic_params)){
    warning(sprintf("There were no insights of problematic continuous parameters. However, the tree topology has mixing and/or convergence issues. This may not be an issue for cases with very small trees (less than a handful of samples). For reference, the trees had %d leaves in these analyses",length(chains[[1]]$trees[[1]]$tip.label))
    )} else {
      warning("WARNING: convergence statistics show convergence problems. You must re-run the analysis after tweaking the MCMC parameters to solve these issues. These results must not be used for biological interpretation or publication")
    }

  #Get marginal likelihood estimates from the posterior sample and write them
  #######################################
  this_lml <- LaplacesDemon::LML(LL=these_data[burnin==FALSE,likelihood],method="HME")$LML
  this_AICm <- AICM(these_data[burnin==FALSE,likelihood])
  mle_table <- data.table::data.table(cond=c(base_name),method=c("HME","AICm"),lML=c(this_lml,this_AICm))

  check_write_csv(mle_table,
                  mle_suffix,
                  out_dir,
                  paste(sep="_",base_name,mle_suffix),
                  "Marginal likelihood estimation deactivated. If you need it, make sure to indicate an output file suffix using the mle_suffix argument")

  #######################################

  #RWTY's citation scheme
  #TODO do we want to do something like this here?
  # print(list(citation("phyfumr"), citation("rwty"), citation("ape"), citation("phangorn"),
  #      citation("ggplot2"), citation("coda"), citation("viridis"),
  #      citation("ggdendro"), citation("GGally"), citation("plyr"),
  #      citation("reshape2"), citation("stats")))

  return(list("all"=convergence,
              "problematic"=convergence[param %in% problematic_params],
              "mle"=mle_table))
}

#' Complete MCMC QC per experiment
#'
#' This function uses [mcmc_qc_condition()] to generate plots and tables to aid
#' assessing phyfum's mixing and convergence. It can use several chains per run,
#' and automatically flags some issues using Rhat and ESS. Depending on your
#' computer resources and the size of your experiment this function may be quite
#' slow and could be worth using [mcmc_qc_condition()] per condition in a HPC,
#' then use [mcmc_qc_gather()] to gather the results like this function does.
#'
#' @details all .trees files within the objective directory will be analyzed.
#'   Files with the same name will be considered independent phyfum runs under
#'   the same conditions.
#'
#' @param file_dir directory that contains .trees and .log files (with or
#'   without subdirectories)
#' @inheritParams mcmc_qc_condition
#' @param n_cores number of cores to use for coarse-grained parallelization (per
#'   patient)
#' @param problematic_output filename of the tabular output with a list of
#'   parameters that may have mixing/convergence problems. If NULL, this output
#'   is not generated.
#' @param all_output filename of the tabular output with all mixing and
#'   convergence statistics. If NULL, this output is not generated..
#' @param mle_output filename of the tabular output with marginal likelihood
#'   estimations using HME and AICm methods.
#' @param conditions_with_issues_output output filename for a summary table that
#'   lists conditions with topological or continuous parameter convergence
#'   issues. If NULL, this output is not generated.
#' @inheritDotParams mcmc_qc_condition burnin_p min_ess max_rhat cred_mass
#'   n_focal_trees treedist n_cores base_name problematic_output_suffix
#'   all_output_suffix posterior_plots_suffix rwty_plots_suffix
#'   correlation_plots_suffix mle_suffix
#' @returns list with data.tables with all = convergence statistics, problematic
#'   = parameter with problems, MLE = mle estimation table data.table. Each has
#'   a column indicating the condition (i.e., filename of the input file without
#'   .trees)
#' @export

mcmc_qc <- function(file_dir,
                    out_dir,
                    plot_dir,
                    n_cores,
                    problematic_output = "problematicParams.csv",
                    all_output = "allParams.csv",
                    mle_output = "MLE.csv",
                    conditions_with_issues_output = "problematicConditions.csv",
                    ...){

  trees_files <- list.files(path = file_dir, pattern = "*.trees$",
                            full.names = T,
                            recursive = T)
  run_names <- basename(trees_files)

  runs_info <- data.table::data.table(trees_files=trees_files,run_names = gsub(pattern = ".trees",
                                                                   replacement = "",
                                                                   x = run_names))

  all_results <- slapply(split(runs_info,by="run_names"),function(this_data,out_dir,plot_dir){
    this_trees_files <- this_data[,trees_files]
    mcmc_qc_condition(trees_files = this_trees_files,
                    out_dir = out_dir,
                    plot_dir = plot_dir,
                    n_cores = 1,
                    ...)
  },out_dir=out_dir,plot_dir=plot_dir,cl=n_cores)

  #binding each of the lists across conditions
  keys <- names(all_results[[1]])
  all_results_merged <- lapply(keys, function(k) {
    data.table::rbindlist(lapply(all_results, `[[`, k), idcol = "condition")
  })
  names(all_results_merged) <- keys

  #mle already had a cond field so I am removing it here
  if(!is.null(all_results_merged$mle))
    all_results_merged$mle[,`:=`(cond=NULL)]


  #writing the combined versions
  check_write_csv(all_results_merged$all,
                  all_output,
                  out_dir)
  check_write_csv(all_results_merged$problematic,
                  problematic_output,
                  out_dir)
  check_write_csv(all_results_merged$mle,
                  mle_output,
                  out_dir)

  #summary with problematic conditions
  check_write_csv(all_results_merged$problematic[,.(topology=length(grep("treeTopology",param))>1,param=length(grep("treeTopology",param,invert = TRUE))>1),by=condition],
                  conditions_with_issues_output,
                  out_dir)

  return(all_results_merged)
}

#' Gather and summarize MCMC QC per experiment
#'
#' This function collects and summarizes the results from [mcmc_qc_condition()].
#' It is an alternative to [mcmc_qc()] when [mcmc_qc_condition()] is run in
#' parallel using multiple computers (typicaly an HPC cluster).
#'
#' @details all .trees files within the objective directory will be analized.
#'   Files with the same name will be considered independent phyfum runs under
#'   the same conditions.
#'
#' @param out_dir output directory where all the [mcmc_qc_condition()] outputs
#'   are located, which will also be where the output of this function will be
#'   located (as [mcmc_qc()])
#' @inheritParams mcmc_qc_condition
#' @inherit mcmc_qc params return
#'
#' @export

mcmc_qc_gather <- function(out_dir,
                    problematic_output_suffix = "problematicParams.csv",
                    all_output_suffix = "allParams.csv",
                    mle_suffix = "MLE.csv",
                    problematic_output = "problematicParams.csv",
                    all_output = "allParams.csv",
                    mle_output = "MLE.csv",
                    conditions_with_issues_output = "problematicConditions.csv"){

  list_names <- c(all_output_suffix,problematic_output_suffix,mle_suffix)
  all_files <- lapply(list_names,function(x,path,full.names){
    these_files <- list.files(pattern=paste0('..*',x),path=path,full.names=full.names)
    names(these_files) <- gsub(pattern = paste0(".",x),replacement = "",x = basename(these_files))
    return(these_files)
    },path=out_dir,full.names=T)
  names(all_files) <- list_names

  all_results_merged <- lapply(list_names,function(type){data.table::rbindlist(lapply(all_files[[type]],data.table::fread),idcol = "condition",fill = TRUE)})
  names(all_results_merged) <- c("all","problematic","mle")

  #mle already had a cond field so I am removing it here
  if(!is.null(all_results_merged$mle))
    all_results_merged$mle[,`:=`(cond=NULL)]

  #writing the combined versions
  check_write_csv(all_results_merged$all,
                  all_output,
                  out_dir)
  check_write_csv(all_results_merged$problematic,
                  problematic_output,
                  out_dir)
  check_write_csv(all_results_merged$mle,
                  mle_output,
                  out_dir)

  #summary with problematic conditions
  check_write_csv(all_results_merged$problematic[,.(topology=length(grep("treeTopology",param))>1,param=length(grep("treeTopology",param,invert = TRUE))>1),by=condition],
                  conditions_with_issues_output,
                  out_dir)

  return(all_results_merged)
}
