
#' mochims_coefficient_size_distribution_datasets
#'
#' Coefficient size distribution plots for multiple datasets.
#'
#' @param dataset_names character vector of dataset names (required)
#' @param base_dir Base directory (required)
#' @param output_dir Output directory (required)
#' @param stagenum stage number (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
mochims_coefficient_size_distribution_datasets <- function(
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
    mochims_coefficient_size_distribution(
      input_file = file.path(base_dir, paste0("001", "_mochims_thermo_model_results_", i), "model_coefficients.txt"),
      outpath = mochims__format_dir(dir_suffix=paste0("_coefficient_size_distribution_", i), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme)
  }

}
