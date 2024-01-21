
#' mochims_structure_metrics_datasets
#'
#' Add 3D structure metrics of multiple datasets.
#'
#' @param dataset_names character vector of dataset names (required)
#' @param pdb_file_list list of paths to PDB files (required)
#' @param base_dir Base directory (required)
#' @param output_dir Output directory (required)
#' @param stagenum stage number (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
mochims_structure_metrics_datasets <- function(
  dataset_names,
  pdb_file_list,
  base_dir,
  output_dir,
  stagenum,
  execute = TRUE
  ){

  for(i in dataset_names){
    mochims_structure_metrics(
      input_file = file.path(base_dir, paste0("001", paste0("_mochims_thermo_model_results_", i)), "model_coefficients.txt"),
      outpath = mochims__format_dir(dir_suffix=paste0("_mochims_structure_metrics_", i), stagenum=stagenum, base_dir=output_dir),
      pdb_file = pdb_file_list[[i]],
      execute = execute)
  }
}