
#' mochims_coefficient_size_distribution
#'
#' Coupling scatterplots.
#'
#' @param input_file path to MoCHI thermo model fit results (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#'
#' @return Nothing
#' @export
#' @import data.table
mochims_coefficient_size_distribution <- function(
  input_file,
	outpath,
  colour_scheme
	){

  #Return if input_file doesn't exist
  if(!file.exists(input_file)){
    return()
  }

	#Load free energies
	dg_dt <- fread(input_file)

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

  #Display status
  message(paste("\n\n*******", paste0("running stage: mochims_coefficient_size_distribution for ", domain_name), "*******\n\n"))

  #Create output directory
  mochims__create_dir(mochims_dir = outpath)

  #Weighted mean absolute value
  dg_dt[, mean_abs := abs(mean)]

  #Coefficient size distribution
  plot_dt <- dg_dt[coef_order!=0][,.(mean_abs, coef_order)]
  plot_dt <- plot_dt[order(mean_abs, decreasing = T)]
  plot_dt[, mean_abs_rank := 1:.N]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mean_abs_rank, mean_abs)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) + 
    ggplot2::geom_vline(xintercept = 50, linetype = 2) + 
    ggplot2::xlab("Rank") +
    ggplot2::ylab("|coefficient|") +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(outpath, "coefficient_scatter_meanabsrank.pdf"), d, width = 3, height = 3, useDingbats=FALSE)

  #Coefficient size distribution
  plot_dt <- dg_dt[coef_order!=0][,.(mean_abs, coef_order)]
  plot_dt <- plot_dt[order(mean_abs, decreasing = T)]
  plot_dt[, mean_abs_rank := 1:.N]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mean_abs_rank, mean_abs)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) + 
    ggplot2::geom_vline(xintercept = 50, linetype = 2) + 
    ggplot2::xlab("Rank") +
    ggplot2::ylab("|coefficient|") +
    ggplot2::scale_x_continuous(trans='log10') +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(outpath, "coefficient_scatter_meanabsranklog10.pdf"), d, width = 3, height = 3, useDingbats=FALSE)

  #Coefficient size distribution - colour by coefficient order
  plot_dt <- dg_dt[coef_order!=0][,.(mean_abs, coef_order)]
  plot_dt <- plot_dt[order(mean_abs, decreasing = T)]
  plot_dt[, mean_abs_rank := 1:.N]
  plot_dt[, coef_order_col := as.factor(coef_order)]
  plot_cols <- c(colour_scheme[["shade 0"]][1], colour_scheme[["shade 0"]][2], colour_scheme[["shade 0"]][4], colour_scheme[["shade 0"]][3], 'grey', 'black')
  names(plot_cols) <- as.character(1:6)
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mean_abs_rank, mean_abs, color = coef_order_col)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) + 
    ggplot2::xlab("Rank") +
    ggplot2::ylab("|coefficient|") +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_manual(values = plot_cols)
  ggplot2::ggsave(file.path(outpath, "coefficient_scatter_meanabsrank_colour.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  return(dg_dt)
}
