
#' mochims__plot_additive_trait_3state
#'
#' Plot folding additive trait.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param input_dt data.table with model free energy estimates (required)
#' @param RT constant (default:0.001987*(273+24))
#' @param trait_names additive trait name (default:c('Folding', 'Binding'))
#' @param phenotype_name phenotype name (default:'Binding')
#' @param report_outpath output path for scatterplots (required)
#' @param dataset_name dataset name and suffix for output files (required)
#' @param colour_scheme colour scheme file (default:ggplot colours)
#'
#' @return Nothing
#' @export
#' @import data.table
mochims__plot_additive_trait_3state <- function(
  mochi_outpath,
  input_dt,
  RT = 0.001987*(273+24),
  trait_names = c("Folding", "Binding"),
  phenotype_name = "Binding",
  report_outpath,
  dataset_name,
  colour_scheme
  ){

  #Plot data.table
  plot_dt <- copy(input_dt)

  #Return if no global epistasis
  if(!'fold_1_additive_trait1' %in% names(plot_dt)){
    return()
  }

  #Folding dG
  plot_dt[, dg0_pred := fold_1_additive_trait0*RT]
  plot_dt[, ddg0_pred := dg0_pred-plot_dt[mut_order==0,dg0_pred]]

  #Binding dG
  plot_dt[, dg1_pred := fold_1_additive_trait1*RT]
  plot_dt[, ddg1_pred := dg1_pred-plot_dt[mut_order==0,dg1_pred]]

  #WT dGs
  dg0_WT <- plot_dt[mut_order==0,dg0_pred][1]
  dg1_WT <- plot_dt[mut_order==0,dg1_pred][1]

  #Number of grid points
  num_grid <- 15

  ### All data
  ###########################

  #Model data.table (for geom_line)
  energy_range0 <- plot_dt[mut_order>0,range(dg0_pred)]
  energy_range1 <- plot_dt[mut_order>0,range(dg1_pred)]
  energy_grid0 <- seq(energy_range0[1], energy_range0[2], (energy_range0[2]-energy_range0[1])/num_grid)
  energy_grid1 <- seq(energy_range1[1], energy_range1[2], (energy_range1[2]-energy_range1[1])/num_grid)
  
  energy_grid_dt <- as.data.table(expand.grid(energy_grid0 = energy_grid0, energy_grid1 = energy_grid1))

  #Predicted fitness
  pred_fitness_list <- mochims__predict_fitness_3state(
    mochi_outpath = mochi_outpath,
    phenotype_name = phenotype_name,
    folding_energy = energy_grid_dt[,energy_grid0],
    binding_energy = energy_grid_dt[,energy_grid1],
    RT = RT)
  #Predicted fitness data table
  pred_fitness_dt <- data.table(
    dg0_pred = energy_grid_dt[,energy_grid0],
    dg1_pred = energy_grid_dt[,energy_grid1],
    fitness = pred_fitness_list[["fitness"]])

  # Cairo::CairoPDF(file = file.path(report_outpath, paste0("dG_observed_", tolower(phenotype_name), "_scatter.pdf")))
  # plot3D::persp3D(
  #   x = energy_grid0, 
  #   y = energy_grid1, 
  #   z = matrix(data=unlist(pred_fitness_dt[,fitness]), nrow=length(energy_grid0), ncol=length(energy_grid0)), 
  #   r=2, shade=0.4, axes=TRUE,scale=TRUE, box=TRUE, nticks=5, ticktype="detailed", colvar=F, col="white", alpha = 0, border=colour_scheme[["shade 0"]][[1]], lwd=0.2,
  #   xlab = "dG Folding",
  #   ylab = "dG Binding",
  #   zlab = "Fitness (Binding)")
  # plot3D::scatter3D(
  #   x = plot_dt[,dg0_pred], 
  #   y = plot_dt[,dg1_pred], 
  #   z = plot_dt[,fitness], 
  #   add = T, col = viridis::viridis(n = 50), alpha = 0.2, cex = 0.2)
  # dev.off()

  Cairo::CairoPDF(file = file.path(report_outpath, paste0("ddG_observed_", tolower(phenotype_name), "_scatter.pdf")))
  plot3D::persp3D(
    x = energy_grid0-dg0_WT, 
    y = energy_grid1-dg1_WT, 
    z = matrix(data=pred_fitness_dt[,fitness], nrow=length(energy_grid0), ncol=length(energy_grid0)), 
    r=2, shade=0.4, axes=TRUE,scale=TRUE, box=TRUE, nticks=5, ticktype="detailed", colvar=F, col="white", alpha = 0, border=colour_scheme[["shade 0"]][[1]], lwd=0.2,
    xlab = "dG Folding",
    ylab = "dG Binding",
    zlab = "Fitness (Binding)")
  plot3D::scatter3D(
    x = plot_dt[,dg0_pred]-dg0_WT, 
    y = plot_dt[,dg1_pred]-dg1_WT, 
    z = plot_dt[,fitness], 
    # add = T, col = viridis::viridis(n = 50), alpha = 0.2, cex = 0.2)
    add = T, col = viridis::viridis(n = 50), alpha = 1, cex = 0.2)
  dev.off()

}
