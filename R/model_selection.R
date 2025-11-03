# R/model_selection.R

#' Selects candidate and best S and performs QC of these selections
#'
#' @details HME (harmonic mean estimator) is always available and, in our
#'   simulations, works well for this specific purpose despite the well-known
#'   problems it has for model selection overall. SS and PS are the preferred
#'   methods but phyfum runs with MLE (marginal likelihood estimation) activated
#'   are required and they take at least twice the time to run. This function
#'   goes down the preference list until finding a method for which it has lML
#'   data for all the S to work and uses that one.
#'
#' @param in_dir input directory where the MLE.csv files resulting from
#'   [mcmc_qc()] or [mcmc_qc_condition()] are located
#' @param out_dir output directory to write .csv files with model selection
#'   information
#' @param plot_dir plot output directory
#' @param method sorted list of preferred methods to estimate the marginal
#'   likelihood for model selection
#' @param min_bf minimum Bayes Factor to consider a model (S value) as a valid
#'   alternative to the best-fit
#' @param file_regex R regular expression to locate the files that contain the
#'   MLE estimation information. The default should work unless [mcmc_qc()]
#'   output filenames were modified. You can use "^MLE.csv$" instead of the
#'   default ".MLE.csv$ to increase the speed if [mcmc_qc()] or
#'   [mcmc_qc_gather()] were run before (instead of just [mcmc_qc_condition()]
#'   distributed)
#' @param n_cells_regex R regular expression with a group that captures the S
#'   number from the MLE.csv filename. Again, the default should work (see
#'   file_regex)
#' @param basename_regex R regular expression to substitute with "" to generate
#'   the cond name from the MLE.csv filename. Again, the default should work
#'   (see file_regex)
#' @param by_patient Logical to activate patient-specific tables (usually not
#'   needed)
#' @param warn_extremes Logical to activate warning messages indicating cases
#'   for which additional phyfum runs with other S values should be performed to
#'   improve S model selection
#' @param plot_aspect_ratio aspect ratio for the S model selection plots
#' @param all_output_suffix filename (or suffix for by_patient output) for the
#'   table with all the selection results
#' @param plot_output_suffix suffix for S model selection plots
#' @param selected_output_suffix filename (or suffix for by_patient output) for
#'   the table with the selected S values (using the BF to the best S) per
#'   patient
#' @param bestS_output_suffix filename (or suffix for by_patient output) for the
#'   table with the best S value per patient
#' @param extreme_selected_output_suffix filename (or suffix for by_patient
#'   output) for the table that contains cases for which a S value that is valid
#'   according to the BF against the best-fit model is a extreme value (i.e.,
#'   more S values should be assessed to know if more S values would be included
#'   in the final selection)
#' @param extreme_valid_output_suffix filename (or suffix for by_patient output)
#'   for the table that contains cases for which the best-fit S value was an
#'   extreme value in the current sampling (i.e., more S values should be
#'   assessed to know if a different S value will better this one)
#'
#' @returns data.table with all selection information
#' @export

model_selection <- function(in_dir,
                            out_dir,
                            plot_dir=NULL,
                            method = c("PS","SS","HME","AICm"),
                            min_bf = 10,
                            file_regex = ".MLE.csv$", #"^MLE.csv$" would use the result of mcmc_qc or mcmc_qc_gather and be a little faster, but this way we do not require that step
                            n_cells_regex = ".*\\.([0-9]+)cells\\..*",
                            basename_regex = "\\.[0-9]+cells\\..*",
                            by_patient = F,
                            warn_extremes = T,
                            plot_aspect_ratio = 4/2,
                            #n_cells_regex <- ".*S([0-9]+).*",
                            #basename_regex <- "S[0-9]+.*",
                            all_output_suffix = "allS.csv",
                            plot_output_suffix = "allS.pdf",
                            selected_output_suffix = "selectedS.csv",
                            bestS_output_suffix = "bestFitS.csv",
                            extreme_selected_output_suffix = "extremeSelectedS.csv",
                            extreme_valid_output_suffix = "extremeValidS.csv"
                            ) {

  min_lbf <- log(min_bf)

  method <- match.arg(method,several.ok = T) #sorted by preference

  if(is.null(plot_dir))
    plot_dir <- out_dir

  input_files <- list.files(path = in_dir,
                            pattern = file_regex,
                            full.names = T)
  if(length(input_files)<1)
    stop("MLE files not found")

  names(input_files) <- basename(input_files)

  #Data parsing
  all_mle_data <- data.table::rbindlist(lapply(input_files,data.table::fread),idcol = "file")
  cond_col <- grep("^cond",colnames(all_mle_data),value=T)
  all_mle_data[,`:=`(patient=gsub(basename_regex,"",get(cond_col)),S=as.numeric(gsub(n_cells_regex,"\\1",get(cond_col))))]

  #Model selection
  all_mle_data_with_selection_stats <- data.table::rbindlist(lapply(split(all_mle_data,by="patient"),function(this_data,min_lbf){
    n_S <- length(this_data[,unique(S)])
    method_analysis <- this_data[,.N,by=method][,`:=`(valid_method=N==n_S)][method,.SD,on="method"][valid_method==T,`:=`(method_order=.I)]
    new_data <- data.table::copy(this_data)
    new_data[method=="AICm",`:=`(lML=lML*-1)] #AICm's lML is the AICm not the lML and thus needs to be inverted for selection. Additionally, we can't use BFs with it.
    new_data <- new_data[order(-lML),.(file,get(cond_col),best=lML==data.table::first(lML),
                                       lML,S,selected=max(lML)-lML<=min_lbf,min_S=min(S),max_S=max(S),
                                       extreme_best=lML==data.table::first(lML)&(S==max(S)|S==min(S)),
                                       extreme_valid=max(lML)-lML<=min_lbf&(S==max(S)|S==min(S))),by=method]
    new_data[method=="AICm",`:=`(lML=lML*-1,selected=NA)]
    return(data.table::merge.data.table(new_data,method_analysis,by = "method",all.x = TRUE))
  },min_lbf=min_lbf),idcol = "patient")

  #Output
  dir.create(out_dir,recursive = T,showWarnings = F)
  dir.create(plot_dir,recursive = T,showWarnings = F)

  check_write_csv(all_mle_data_with_selection_stats[order(method_order),][best==T,.(patient,method,S,lML,extreme_best)],
                  bestS_output_suffix,
                  out_dir)

  check_write_csv(all_mle_data_with_selection_stats[order(method_order),][selected==T,.(patient,method,S,lML,extreme_valid)],
                  selected_output_suffix,
                  out_dir)

  check_write_csv(all_mle_data_with_selection_stats[order(method_order),][,.(patient,method,best,lML,S,selected,valid_method,method_order,extreme_best,extreme_valid)],
                  all_output_suffix,
                  out_dir)

  check_write_csv(all_mle_data_with_selection_stats[extreme_best==T,],
                  extreme_selected_output_suffix,
                  out_dir)

  check_write_csv(all_mle_data_with_selection_stats[extreme_valid==T & method!="AICm",],
                  extreme_valid_output_suffix,
                  out_dir)

  if(by_patient==TRUE){
    for(this_patient in all_mle_data_with_selection_stats[,unique(patient)]){
      check_write_csv(all_mle_data_with_selection_stats[order(method_order),][best==T & patient==this_patient,.(patient,method,S,lML,extreme_best)],
                      bestS_output_suffix,
                      out_dir,
                      paste(sep="_",patient,bestS_output_suffix))

      check_write_csv(all_mle_data_with_selection_stats[order(method_order),][selected==T & patient==this_patient,.(patient,method,S,lML,extreme_valid)],
                      selected_output_suffix,
                      out_dir,
                      paste(sep="_",patient,selected_output_suffix))

      check_write_csv(all_mle_data_with_selection_stats[extreme_best==T & patient==this_patient,],
                      extreme_selected_output_suffix,
                      out_dir,
                      paste(sep="_",patient,extreme_selected_output_suffix))

      check_write_csv(all_mle_data_with_selection_stats[extreme_valid==T & method!="AICm" & patient==this_patient,],
                      extreme_valid_output_suffix,
                      out_dir,
                      paste(sep="_",patient,extreme_valid_output_suffix))
    }
  }

  #Warnings
  if(warn_extremes==T){
    if(length(all_mle_data_with_selection_stats[extreme_best==T & method_order==1 & S==max_S,patient])!=0)
      warning(sprintf("IMPORTANT: The best-fit S value is a extreme of the range explored for patients: %s. The analysis should be re-run adding higher S values to find the global optimum",paste(collapse=",",all_mle_data_with_selection_stats[extreme_best==T & method_order==1 & S==max_S,patient])))

    if(length(all_mle_data_with_selection_stats[extreme_best==T & method_order==1 & S==min_S,patient])!=0)
      warning(sprintf("IMPORTANT: The best-fit S value is a extreme of the range explored for patients: %s. The analysis should be re-run adding lower S values to find the global optimum",paste(collapse=",",all_mle_data_with_selection_stats[extreme_best==T & method_order==1 & S==min_S,patient])))

    if(length(all_mle_data_with_selection_stats[extreme_valid==T & method_order==1 & S==max_S,patient])!=0)
      warning(sprintf("The largest assayed S value cannot be rejected as a valid S value for patients: %s. We recommend the analysis is re-run adding larger S values to describe better the range of S values that may be valid",paste(collapse=",",all_mle_data_with_selection_stats[extreme_valid==T & method_order==1 & S==max_S,patient])))

    if(length(all_mle_data_with_selection_stats[extreme_valid==T & method_order==1 & S==min_S,patient])!=0)
      warning(sprintf("The smallest assayed S value cannot be rejected as a valid S value for patients: %s. We recommend the analysis is re-run adding smaller S values to describe better the range of S values that may be valid",paste(collapse=",",all_mle_data_with_selection_stats[extreme_valid==T & method_order==1 & S==min_S,patient])))
  }

 #plots
  if(!is.null(plot_output_suffix)){
    for(this_patient in all_mle_data_with_selection_stats[,unique(patient)]){

      these_thresholds <- all_mle_data_with_selection_stats[patient==this_patient & method!="AICm",.(minY=max(lML)-min_lbf),by=method]

      this_patient_plot <- ggplot2::ggplot(all_mle_data_with_selection_stats[patient==this_patient & method!="AICm",],ggplot2::aes(x=S,y=lML,color=method))+
        ggplot2::geom_point(alpha=0.8) +
        ggplot2::geom_smooth(alpha=0.2,linewidth=0.3) +
        ggplot2::geom_hline(data=these_thresholds,ggplot2::aes(yintercept=minY,color=method,linetype="BF10")) +
        ggplot2::scale_y_continuous(name="Marginal Likelihood Estimation (logL)",n.breaks=10) +
        ggplot2::scale_x_continuous(name="Number of stem cells (S)") +
        ggplot2::scale_color_brewer(name="Estimation\nmethod",type = "qual",palette = 6) +
        ggplot2::scale_linetype_manual(name=NULL,values=2)+
        cowplot::theme_cowplot()

      suppressMessages(cowplot::save_plot(file=paste(sep="/",plot_dir,paste(sep="_",this_patient,plot_output_suffix)),this_patient_plot,base_height = 6,base_asp = plot_aspect_ratio))
    }
  }

  return(all_mle_data_with_selection_stats)
}
