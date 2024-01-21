
#' mochims__predict_fitness_3state
#'
#' Predict fitness from folding and binding energies for 3-state model.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param phenotype_name phenotype name (required)
#' @param folding_energy folding energies for 3-state model (required)
#' @param binding_energy binding energies for 3-state model (required for binding fitness prediction)
#' @param RT constant (default:0.001987*(273+24))
#'
#' @return A list with folding and binding fitness predictions
#' @export
#' @import data.table
mochims__predict_fitness_3state <- function(
  mochi_outpath,
  phenotype_name,
  folding_energy,
  binding_energy,
  RT = 0.001987*(273+24)
  ){

  #Binding fitness
  if(!is.null(folding_energy) & !is.null(binding_energy) & length(folding_energy) == length(binding_energy)){
    #Linear transformation parameters
    linears_file <- file.path(mochi_outpath, "task_1", 'weights', paste0("linears_weights_", phenotype_name, ".txt"))
    linears_params <- as.list(fread(linears_file)[1,2:3])
    #Binding fitness
    pred_list <- list()
    pred_list[["fraction_state3"]] <- mochims__fraction_bound(folding_energy, binding_energy, RT)
    pred_list[["fitness"]] <- pred_list[["fraction_state3"]] * linears_params[["kernel"]] + linears_params[["bias"]]
  }

  #Return
  return(pred_list)
}