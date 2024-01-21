
#' mochims__plot_additive_trait_3state_sos
#'
#' Plot sum of sigmoids additive trait.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param input_dt data.table with model free energy estimates (required)
#' @param phenotype_name phenotype name (default:'Binding')
#' @param report_outpath output path for scatterplots (required)
#' @param dataset_name dataset name and suffix for output files (required)
#' @param colour_scheme colour scheme file (default:ggplot colours)
#'
#' @return Nothing
#' @export
#' @import data.table
mochims__plot_additive_trait_3state_sos <- function(
  mochi_outpath,
  input_dt,
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

  #Number of grid points
  num_grid <- 15

  ### All data
  ###########################

  #Model data.table (for geom_line)
  energy_range0 <- plot_dt[mut_order>0,range(fold_1_additive_trait0)]
  energy_range1 <- plot_dt[mut_order>0,range(fold_1_additive_trait1)]
  energy_grid0 <- seq(energy_range0[1], energy_range0[2], (energy_range0[2]-energy_range0[1])/num_grid)
  energy_grid1 <- seq(energy_range1[1], energy_range1[2], (energy_range1[2]-energy_range1[1])/num_grid)
  
  energy_grid_dt <- as.data.table(expand.grid(fold_1_additive_trait0 = energy_grid0, fold_1_additive_trait1 = energy_grid1))
  #Save energy grid
  write.table(energy_grid_dt, file = file.path(report_outpath, "energy_grid.txt"), sep = "\t", row.names = F, quote = F)

  #Read energy grid and predictions
  pred_fitness_dt <- fread(file.path(mochi_outpath, "task_1", "predictions", "energy_grid_fitness.txt"))

  Cairo::CairoPDF(file = file.path(report_outpath, paste0("dG_observed_", tolower(phenotype_name), "_scatter.pdf")))
  plot3D::persp3D(
    x = energy_grid0, 
    y = energy_grid1, 
    z = matrix(data=unlist(pred_fitness_dt[,fitness]), nrow=length(energy_grid0), ncol=length(energy_grid0)), 
    r=2, shade=0.4, axes=TRUE,scale=TRUE, box=TRUE, nticks=5, ticktype="detailed", colvar=F, col="white", alpha = 0, border=colour_scheme[["shade 0"]][[1]], lwd=0.2,
    xlab = "dG Folding",
    ylab = "dG Binding",
    zlab = "Fitness (Binding)")
  plot3D::scatter3D(
    x = plot_dt[fitness<max(pred_fitness_dt[,fitness]+0.1) & fitness>min(pred_fitness_dt[,fitness]-0.1),fold_1_additive_trait0], 
    y = plot_dt[fitness<max(pred_fitness_dt[,fitness]+0.1) & fitness>min(pred_fitness_dt[,fitness]-0.1),fold_1_additive_trait1], 
    z = plot_dt[fitness<max(pred_fitness_dt[,fitness]+0.1) & fitness>min(pred_fitness_dt[,fitness]-0.1),fitness], 
    # z = plot_dt[,fold_1], 
    add = T, col = viridis::viridis(n = 50), alpha = 1, cex = 0.2)
  dev.off()

}
