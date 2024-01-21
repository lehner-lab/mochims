
#' mochims_position_class_violins
#'
#' Coupling scatterplots.
#'
#' @param input_file path to MoCHI thermo model fit results (required)
#' @param annotation_file path to MoCHI thermo model fit results (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#'
#' @return Nothing
#' @export
#' @import data.table
mochims_position_class_violins <- function(
  input_file,
  annotation_file,
	outpath,
  colour_scheme
	){

  #Return if input_file doesn't exist
  if(!file.exists(input_file)){
    return()
  }

	#Load free energies
	dg_dt <- fread(input_file)

  #Load annotations
  anno_dt <- fread(annotation_file)

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

  #Display status
  message(paste("\n\n*******", paste("running stage: mochims_position_class_violins for", domain_name), "*******\n\n"))

  #Create output directory
  mochims__create_dir(mochims_dir = outpath)

  #Merge annotations
  dg_dt <- merge(dg_dt, anno_dt, by = "Pos")

  #Free energy densities
  plot_dt <- copy(dg_dt)
  plot_dt[, ddg := .SD[[1]],,.SDcols = "mean_kcal/mol"]
  plot_dt[, Pos_class_plot := factor(Pos_class, levels = c("core", "salt_bridge", "far_side"))]
  plot_cols <- c(colour_scheme[["shade 0"]][1], colour_scheme[["shade 0"]][2], 'grey')
  names(plot_cols) <- c("core", "salt_bridge", "far_side")
  plot_alpha <- c(1, 0.6)
  names(plot_alpha) <- c("FOS", "JUN")
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(Pos_class_plot, ddg, fill = Pos_class_plot, alpha = protein)) +
    ggplot2::geom_violin() +
    ggplot2::xlab("Position class") +
    ggplot2::ylab("Binding free energy") +
    ggplot2::theme_classic() +
    ggplot2::scale_fill_manual(values = plot_cols) +
    ggplot2::scale_alpha_manual(values = plot_alpha)
  ggplot2::ggsave(file.path(outpath, "position_class_violins.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

}
