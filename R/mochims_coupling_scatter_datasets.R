
#' mochims_coupling_scatter_datasets
#'
#' Coupling scatterplots for multiple datasets.
#'
#' @param dataset_names character vector of dataset names (required)
#' @param validation_list list of validation couplings (required)
#' @param base_dir Base directory (required)
#' @param output_dir Output directory (required)
#' @param stagenum stage number (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
mochims_coupling_scatter_datasets <- function(
  dataset_names,
  validation_list,
  base_dir,
  output_dir,
  stagenum,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  for(i in dataset_names){
    mochims_coupling_scatter(
      input_file = file.path(base_dir, paste0("002", "_mochims_structure_metrics_", i), "model_coefficients.txt"),
      validation_couplings = validation_list[[i]],
      outpath = mochims__format_dir(dir_suffix=paste0("_mochims_coupling_scatter_", i), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme)
  }

  outpath <- mochims__format_dir(dir_suffix=paste0("_mochims_coupling_scatter_datasets"), stagenum=stagenum, base_dir=output_dir)
  
  #Create output directory
  mochims__create_dir(mochims_dir = outpath)

  ### Coupling correlations - FAS vs. FAS mechanistic
  ###########################

  #FAS
  dg_dt_fas <- fread(file.path(base_dir, paste0("002", "_mochims_structure_metrics_", 'FAS'), "model_coefficients.txt"))

  #FAS - mechanistic
  dg_dt_fasm <- fread(file.path(base_dir, paste0("002", "_mochims_structure_metrics_", 'FAS-mechanistic'), "model_coefficients.txt"))

  #Coupling scatter comparison
  plot_dt <- merge(dg_dt_fas[coef_order!=0,.(id_ref, mean, ci95)], dg_dt_fasm[coef_order!=0,.(id_ref, meanm = mean, ci95m = ci95)], by = 'id_ref')
  plot_cols <- c('black', colour_scheme[["shade 0"]][1])
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(meanm, mean)) +
    ggplot2::geom_point() +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = mean-ci95/2, ymax = mean+ci95/2), alpha = 1/4) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(xmin = meanm-ci95m/2, xmax = meanm+ci95m/2), alpha = 1/4) +
    # ggplot2::geom_abline(linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(mean, meanm, use = "pairwise.complete"), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Mechanistic model coefficients") +
    ggplot2::ylab("Model coefficients") +
    ggplot2::geom_smooth(formula = 'y ~ x', linetype = 1, method = "lm", color = colour_scheme[["shade 0"]][1], se = T) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(outpath, "FAS_coupling_scatter_comparison.pdf"), d, width = 3, height = 3, useDingbats=FALSE)

  #Coupling scatter comparison - inverted
  plot_dt <- merge(dg_dt_fas[coef_order!=0,.(id_ref, mean = -mean, ci95)], dg_dt_fasm[coef_order!=0,.(id_ref, meanm = mean, ci95m = ci95)], by = 'id_ref')
  plot_cols <- c('black', colour_scheme[["shade 0"]][1])
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(meanm, mean)) +
    ggplot2::geom_point() +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = mean-ci95/2, ymax = mean+ci95/2), alpha = 1/4) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(xmin = meanm-ci95m/2, xmax = meanm+ci95m/2), alpha = 1/4) +
    # ggplot2::geom_abline(linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(mean, meanm, use = "pairwise.complete"), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Mechanistic model coefficients") +
    ggplot2::ylab("Model coefficients") +
    ggplot2::geom_smooth(formula = 'y ~ x', linetype = 1, method = "lm", color = colour_scheme[["shade 0"]][1], se = T) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(outpath, "FAS_coupling_scatter_comparison_inv.pdf"), d, width = 3, height = 3, useDingbats=FALSE)

  ### Coupling correlations - tRNA vs. tRNA linear
  ###########################

  #tRNA
  dg_dt_trna <- fread(file.path(base_dir, paste0("002", "_mochims_structure_metrics_", 'tRNA'), "model_coefficients.txt"))

  #FAS - linear
  dg_dt_trnal <- fread(file.path(base_dir, paste0("002", "_mochims_structure_metrics_", 'tRNA-linear'), "model_coefficients.txt"))

  #Coupling scatter comparison
  plot_dt <- merge(dg_dt_trna[coef_order!=0,.(id_ref, mean, ci95)], dg_dt_trnal[coef_order!=0,.(id_ref, meanm = mean, ci95m = ci95)], by = 'id_ref')
  plot_cols <- c('black', colour_scheme[["shade 0"]][1])
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(meanm, mean)) +
    ggplot2::geom_point() +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = mean-ci95/2, ymax = mean+ci95/2), alpha = 1/4) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(xmin = meanm-ci95m/2, xmax = meanm+ci95m/2), alpha = 1/4) +
    # ggplot2::geom_abline(linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(mean, meanm, use = "pairwise.complete"), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Linear model coefficients") +
    ggplot2::ylab("Model coefficients") +
    ggplot2::geom_smooth(formula = 'y ~ x', linetype = 1, method = "lm", color = colour_scheme[["shade 0"]][1], se = T) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(outpath, "tRNA_coupling_scatter_comparison.pdf"), d, width = 3, height = 3, useDingbats=FALSE)

  #Coupling scatter comparison - inverted
  plot_dt <- merge(dg_dt_trna[coef_order!=0,.(id_ref, mean = -mean, ci95)], dg_dt_trnal[coef_order!=0,.(id_ref, meanm = mean, ci95m = ci95)], by = 'id_ref')
  plot_cols <- c('black', colour_scheme[["shade 0"]][1])
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(meanm, mean)) +
    ggplot2::geom_point() +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = mean-ci95/2, ymax = mean+ci95/2), alpha = 1/4) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(xmin = meanm-ci95m/2, xmax = meanm+ci95m/2), alpha = 1/4) +
    # ggplot2::geom_abline(linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(mean, meanm, use = "pairwise.complete"), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Linear model coefficients") +
    ggplot2::ylab("Model coefficients") +
    ggplot2::geom_smooth(formula = 'y ~ x', linetype = 1, method = "lm", color = colour_scheme[["shade 0"]][1], se = T) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(outpath, "tRNA_coupling_scatter_comparison_inv.pdf"), d, width = 3, height = 3, useDingbats=FALSE)

  ### Coupling correlations - tRNA 2dim Activity1 vs. tRNA linear
  ###########################

  #tRNA 2dim
  dg_dt_trna <- fread(file.path(base_dir, paste0("002", "_mochims_structure_metrics_", 'tRNA-2dim'), "model_coefficients.txt"))[trait_name=='Activity1']

  #FAS - linear
  dg_dt_trnal <- fread(file.path(base_dir, paste0("002", "_mochims_structure_metrics_", 'tRNA-linear'), "model_coefficients.txt"))

  #Coupling scatter - backbone distance
  plot_dt <- merge(dg_dt_trna[coef_order!=0,.(id_ref, mean, ci95)], dg_dt_trnal[coef_order!=0,.(id_ref, meanm = mean, ci95m = ci95)], by = 'id_ref')
  plot_cols <- c('black', colour_scheme[["shade 0"]][1])
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(meanm, mean)) +
    ggplot2::geom_point() +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = mean-ci95/2, ymax = mean+ci95/2), alpha = 1/4) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(xmin = meanm-ci95m/2, xmax = meanm+ci95m/2), alpha = 1/4) +
    # ggplot2::geom_abline(linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(mean, meanm, use = "pairwise.complete"), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Linear model coefficients") +
    ggplot2::ylab("Model coefficients") +
    ggplot2::geom_smooth(formula = 'y ~ x', linetype = 1, method = "lm", color = colour_scheme[["shade 0"]][1], se = T) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(outpath, "tRNA-2dim_Activity1_coupling_scatter_comparison.pdf"), d, width = 3, height = 3, useDingbats=FALSE)

  ### Coupling correlations - tRNA 2dim Activity2 vs. tRNA linear
  ###########################

  #tRNA 2dim
  dg_dt_trna <- fread(file.path(base_dir, paste0("002", "_mochims_structure_metrics_", 'tRNA-2dim'), "model_coefficients.txt"))[trait_name=='Activity2']

  #FAS - linear
  dg_dt_trnal <- fread(file.path(base_dir, paste0("002", "_mochims_structure_metrics_", 'tRNA-linear'), "model_coefficients.txt"))

  #Coupling scatter - backbone distance
  plot_dt <- merge(dg_dt_trna[coef_order!=0,.(id_ref, mean, ci95)], dg_dt_trnal[coef_order!=0,.(id_ref, meanm = mean, ci95m = ci95)], by = 'id_ref')
  plot_cols <- c('black', colour_scheme[["shade 0"]][1])
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(meanm, mean)) +
    ggplot2::geom_point() +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = mean-ci95/2, ymax = mean+ci95/2), alpha = 1/4) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(xmin = meanm-ci95m/2, xmax = meanm+ci95m/2), alpha = 1/4) +
    # ggplot2::geom_abline(linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(mean, meanm, use = "pairwise.complete"), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Linear model coefficients") +
    ggplot2::ylab("Model coefficients") +
    ggplot2::geom_smooth(formula = 'y ~ x', linetype = 1, method = "lm", color = colour_scheme[["shade 0"]][1], se = T) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(outpath, "tRNA-2dim_Activity2_coupling_scatter_comparison.pdf"), d, width = 3, height = 3, useDingbats=FALSE)

}
