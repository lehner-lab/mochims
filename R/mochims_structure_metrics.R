
#' mochims_structure_metrics
#'
#' Add 3D structure metrics.
#'
#' @param input_file path to MoCHI thermo model fit results (required)
#' @param outpath output path for plots and saved objects (required)
#' @param pdb_file path to PDB file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
mochims_structure_metrics <- function(
  input_file,
  outpath,
  pdb_file,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute | !file.exists(input_file)){
    return()
  }

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

  #Display status
  message(paste("\n\n*******", paste("running stage: mochims_structure_metrics for", domain_name), "*******\n\n"))

  #Create output directory
  mochims__create_dir(mochims_dir = outpath)

  if(!is.null(pdb_file)){
    #Calculate pairwise distances - scHAmin
    dist_dt_scHA <- mochims__distance_matrix_from_PDB(
      input_file = pdb_file,
      metric = 'scHA')
    # dist_dt_calpha <- mochims__distance_matrix_from_PDB(
    #   input_file = pdb_file,
    #   metric = 'calpha')

    #Merge distances and contacts
    # dist_dt <- merge(dist_dt_scHA, dist_dt_calpha, by = c("Pos1", "Pos2"))
    dist_dt <- dist_dt_scHA
    dist_dt[, Pos_ref := paste0(Pos1, "_", Pos2)]
  }

  #Load free energies - order 0
  dg_dt <- fread(input_file)
  dg_dt[, Pos_ref := as.character(Pos_ref)]
  dg_dt_list <- list()
  dg_dt_list[['0']] <- dg_dt[coef_order==0]

  #Load free energies - order 1
  dg_dt_list[['1']] <- dg_dt[coef_order==1]

  #Load free energies - order 2
  dg_dt_list[['2']] <- dg_dt[coef_order==2]

  if(dg_dt_list[['2']][,.N]!=0){
    #Merge with free energies
    dg_dt_list[['2']][, Pos1 := as.integer(sapply(strsplit(Pos, "_"), '[', 1))]
    dg_dt_list[['2']][, Pos2 := as.integer(sapply(strsplit(Pos, "_"), '[', 2))]
    dg_dt_list[['2']][, backbone := abs(Pos1-Pos2)]
    if(!is.null(pdb_file)){
      dg_dt_list[['2']] <- merge(dg_dt_list[['2']], dist_dt[,.SD,,.SDcols = names(dist_dt)[!names(dist_dt) %in% c("Pos1", "Pos2")]], by = c("Pos_ref"), all.x = T)
    }
  }

  #Save dGs and ddGs
  write.table(rbindlist(dg_dt_list, fill = T),
    file = file.path(outpath, "model_coefficients.txt"), 
    quote = F, sep = "\t", row.names = F)
}
