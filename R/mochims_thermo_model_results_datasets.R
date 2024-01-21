
#' mochims_thermo_model_results_datasets
#'
#' Evaluate thermo model results of multiple datasets.
#'
#' @param dataset_names character vector of dataset names (required)
#' @param literature_free_energies path to literature free energies (required)
#' @param base_dir Base directory (required)
#' @param output_dir Output directory (required)
#' @param stagenum stage number (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
mochims_thermo_model_results_datasets <- function(
  dataset_names,
  literature_free_energies,
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
    mochims_thermo_model_results(
      mochi_outpath = file.path(base_dir, "Data", "mochi", i),
      literature_free_energies = file.path(literature_free_energies, paste0(i, "_literature_free_energies.txt")),
      normalisation_dir = file.path(base_dir, "Data", "normalisation", i),
      outpath = mochims__format_dir(dir_suffix=paste0("_mochims_thermo_model_results_", i), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme,
      execute = execute)
  }

}
