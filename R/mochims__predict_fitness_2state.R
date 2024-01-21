
#' mochims__predict_fitness_2state
#'
#' Predict fitness from energies for 2-state model.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param phenotype_name phenotype name (required)
#' @param energy energies for 2-state model (required)
#' @param RT constant (default:0.001987*(273+24))
#'
#' @return A list with fitness predictions
#' @export
#' @import data.table
mochims__predict_fitness_2state <- function(
  mochi_outpath,
  phenotype_name,
  energy,
  RT = 0.001987*(273+24)
  ){

  #Linear transformation parameters
  linears_file <- file.path(mochi_outpath, "task_1", 'weights', paste0("linears_weights_", phenotype_name, ".txt"))
  linears_params <- as.list(fread(linears_file)[1,2:3])

  #Predicted fitness
  pred_list <- list()

  #Fitness
  pred_list[["fraction_state2"]] <- mochims__fraction_folded(energy, RT)
  pred_list[["fitness"]] <- pred_list[["fraction_state2"]] * linears_params[["kernel"]] + linears_params[["bias"]]

  #Return
  return(pred_list)
}