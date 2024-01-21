
#' mochims__plot_additive_trait_2state
#'
#' Plot additive trait for 2-state model.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param input_dt data.table with model free energy estimates (required)
#' @param RT constant (default:0.001987*(273+24))
#' @param trait_name additive trait name (default:'Folding')
#' @param phenotype_name phenotype name (default:'Abundance')
#' @param report_outpath output path for scatterplots (required)
#' @param dataset_name dataset name and suffix for output files (required)
#' @param colour_scheme colour scheme file (default:ggplot colours)
#'
#' @return Nothing
#' @export
#' @import data.table
mochims__plot_additive_trait_2state <- function(
  mochi_outpath,
  input_dt,
  RT = 0.001987*(273+24),
  trait_name = "Folding",
  phenotype_name = "Abundance",
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

  #dG
  plot_dt[, dg_pred := fold_1_additive_trait0*RT]
  plot_dt[, ddg_pred := dg_pred-plot_dt[mut_order==0,dg_pred]]

  #Model data.table (for geom_line)
  energy_range <- plot_dt[mut_order>0,range(dg_pred)]
  energy_grid <- seq(energy_range[1], energy_range[2], (energy_range[2]-energy_range[1])/500)
  #Predicted fitness
  pred_fitness_list <- mochims__predict_fitness_2state(
    mochi_outpath = mochi_outpath,
    phenotype_name = phenotype_name,
    energy = energy_grid,
    RT = RT)
  #Predicted fitness data table
  pred_fitness_dt <- data.table(
    dg_pred = rep(energy_grid, 2),
    ddg_pred = rep(energy_grid, 2)-plot_dt[mut_order==0,dg_pred][1],
    fitness = pred_fitness_list[["fitness"]],
    mut_order = rep(c(1, 2), each = length(energy_grid)))
  
  #Binhex no facet
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(dg_pred, fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::xlab(bquote(.(trait_name) ~ ""*Delta*"G (inferred)")) +
    ggplot2::ylab(paste0(phenotype_name, " fitness\n(observed)")) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,dg_pred][1], color = "grey", linetype = 2) +
    ggplot2::geom_line(data = pred_fitness_dt[mut_order==1], color = colour_scheme[["shade 0"]][[1]]) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("dG_observed_", tolower(trait_name), "_binhex.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Binhex no facet
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(dg_pred, fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::xlab(bquote(.(trait_name) ~ ""*Delta*"G (inferred)")) +
    ggplot2::ylab(paste0(phenotype_name, " fitness\n(observed)")) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,dg_pred][1], color = "grey", linetype = 2) +
    ggplot2::geom_line(data = pred_fitness_dt[mut_order==1], color = colour_scheme[["shade 0"]][[1]]) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("dG_observed_", tolower(trait_name), "_binhex_nolog.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Binhex no facet
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(dg_pred, fitness)) +
    ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::xlab(bquote(.(trait_name) ~ ""*Delta*"G (inferred)")) +
    ggplot2::ylab(paste0(phenotype_name, " fitness\n(observed)")) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,dg_pred][1], color = "grey", linetype = 2) +
    ggplot2::geom_line(data = pred_fitness_dt[mut_order==1], color = colour_scheme[["shade 0"]][[1]]) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("dG_observed_", tolower(trait_name), "_binhex_highres.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Binhex no facet
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(dg_pred, fitness)) +
    ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::xlab(bquote(.(trait_name) ~ ""*Delta*"G (inferred)")) +
    ggplot2::ylab(paste0(phenotype_name, " fitness\n(observed)")) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,dg_pred][1], color = "grey", linetype = 2) +
    ggplot2::geom_line(data = pred_fitness_dt[mut_order==1], color = colour_scheme[["shade 0"]][[1]]) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("dG_observed_", tolower(trait_name), "_binhex_highres_nolog.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Binhex no facet
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(ddg_pred, fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::xlab(bquote(.(trait_name) ~ ""*Delta*Delta*"G (inferred)")) +
    ggplot2::ylab(paste0(phenotype_name, " fitness\n(observed)")) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_line(data = pred_fitness_dt[mut_order==1], color = colour_scheme[["shade 0"]][[1]]) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("ddG_observed_", tolower(trait_name), "_binhex.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Binhex no facet
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(ddg_pred, fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::xlab(bquote(.(trait_name) ~ ""*Delta*Delta*"G (inferred)")) +
    ggplot2::ylab(paste0(phenotype_name, " fitness\n(observed)")) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_line(data = pred_fitness_dt[mut_order==1], color = colour_scheme[["shade 0"]][[1]]) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("ddG_observed_", tolower(trait_name), "_binhex_nolog.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Binhex no facet
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(ddg_pred, fitness)) +
    ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::xlab(bquote(.(trait_name) ~ ""*Delta*Delta*"G (inferred)")) +
    ggplot2::ylab(paste0(phenotype_name, " fitness\n(observed)")) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_line(data = pred_fitness_dt[mut_order==1], color = colour_scheme[["shade 0"]][[1]]) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("ddG_observed_", tolower(trait_name), "_binhex_highres.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Binhex no facet
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(ddg_pred, fitness)) +
    ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::xlab(bquote(.(trait_name) ~ ""*Delta*Delta*"G (inferred)")) +
    ggplot2::ylab(paste0(phenotype_name, " fitness\n(observed)")) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_line(data = pred_fitness_dt[mut_order==1], color = colour_scheme[["shade 0"]][[1]]) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, paste0("ddG_observed_", tolower(trait_name), "_binhex_highres_nolog.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

}
