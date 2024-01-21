
#' mochims_position_class_violins_datasets
#'
#' Coupling scatterplots for multiple datasets.
#'
#' @param dataset_names character vector of dataset names (required)
#' @param annotations Annotations directory (required)
#' @param base_dir Base directory (required)
#' @param output_dir Output directory (required)
#' @param stagenum stage number (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
mochims_position_class_violins_datasets <- function(
  dataset_names,
  annotations,
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
    #order 1
    mochims_position_class_violins(
      input_file = file.path(base_dir, paste0("001", "_mochims_thermo_model_results_", i), "model_coefficients.txt"),
      annotation_file = file.path(annotations, paste0(i, "_annotations.txt")),
      outpath = mochims__format_dir(dir_suffix=paste0("_mochims_position_class_violins_", i), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme)
  }

}
