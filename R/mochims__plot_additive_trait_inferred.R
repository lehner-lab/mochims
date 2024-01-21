
#' mochims__plot_additive_trait_inferred
#'
#' Plot inferred global epistasis.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param input_dt data.table with model free energy estimates (required)
#' @param phenotype_name phenotype name (default:'PSI')
#' @param report_outpath output path for scatterplots (required)
#' @param dataset_name dataset name and suffix for output files (required)
#' @param colour_scheme colour scheme file (default:ggplot colours)
#'
#' @return Nothing
#' @export
#' @import data.table
mochims__plot_additive_trait_inferred <- function(
  mochi_outpath,
  input_dt,
  phenotype_name = "PSI",
  report_outpath,
  dataset_name,
  colour_scheme
  ){

  #Plot data.table
  plot_dt <- copy(input_dt)

  #Return if no global epistasis
  if(!'fold_1_additive_trait0' %in% names(plot_dt)){
    return()
  }

  pred_fitness_dt <- copy(input_dt[,.(fold_1_additive_trait0, fitness=fold_1)])
  pred_fitness_dt <- pred_fitness_dt[order(fold_1_additive_trait0, decreasing = F)]

  #Binhex no facet
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(fold_1_additive_trait0, fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::xlab("Additive trait (inferred)") +
    ggplot2::ylab(paste0(phenotype_name, " fitness\n(observed)")) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,fold_1_additive_trait0][1], color = "grey", linetype = 2) +
    ggplot2::geom_line(data = pred_fitness_dt, color = colour_scheme[["shade 0"]][[1]]) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("dG_observed_", "additive_trait", "_binhex.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

}
