
#' mochims_thermo_model_results
#'
#' Evaluate thermo model results.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param temperature temperature in degrees celcuis (default:30)
#' @param literature_free_energies path to literature free energies (default:NA)
#' @param normalisation_dir path to normalisation directory (required)
#' @param position_offset residue position offset (default:0)
#' @param position_dict residue position dictionary (default:{})
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
mochims_thermo_model_results <- function(
  mochi_outpath,
  temperature = 30,
  literature_free_energies = NA,
  normalisation_dir,
  position_offset = 0,
  position_dict = list(),
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed or mochi_outpath doesn't exist
  if(!execute | !dir.exists(mochi_outpath)){
    return()
  }

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

  #Display status
  message(paste("\n\n*******", paste("running stage: mochims_thermo_model_results for", domain_name), "*******\n\n"))

  #Create output directory
  mochims__create_dir(mochims_dir = outpath)

  #Constants
  gas_constant <- 0.001987
  RT <- gas_constant*(273+temperature)

  #Load normalisation data
  normalisation_file <- ""
  if(dir.exists(normalisation_dir)){
    normalisation_file <- file.path(normalisation_dir, list.files(normalisation_dir, '.txt$')[1])
  }

  #Load model results
  model_results <- mochims__get_model_results(
    input_folder = mochi_outpath, 
    normalisation_file = normalisation_file,
    # input_dt = fitness_dt, 
    RT = RT)
  pred_dt <- model_results[['pred']]
  coef_dt <- model_results[['coef']]

  #Add model name
  coef_dt[, dataset_model := basename(mochi_outpath)]

  #Translate positions
  if(length(position_dict)!=0){
    for(i in rev(names(position_dict))){
      #Translate Pos_ref
      for(j in c(paste0('(^)', i, '($)'), paste0('(^)', i, '(_)'), paste0('(_)', i, '(_)'), paste0('(_)', i, '($)'))){
        coef_dt[, Pos_ref := gsub(j, paste0("\\1", position_dict[[i]], "\\2"), Pos_ref)]
      }
      #Translate id_ref
      coef_dt[, id_ref := gsub(paste0('([A-Za-z])', i, '([A-Za-z])'), paste0('\\1', position_dict[[i]], '\\2'), id_ref)]
    }
  }

  #Call confident ddGs
  coef_dt <- mochims__define_confident_free_energies(
    input_dt = coef_dt, 
    report_outpath = outpath, 
    highlight_colour = colour_scheme[["shade 0"]][[1]])

  #Plot model performance
  dataset_names <- names(pred_dt)[grepl("Abundance|PSI|Fitness|Brightness|Binding", names(pred_dt))]
  for(i in dataset_names){
    mochims__plot_model_performance(
      input_dt = pred_dt[get(i)==1,], 
      report_outpath = outpath, 
      suffix = i,
      highlight_colour = colour_scheme[["shade 0"]][[1]])
  }

  #Plot 2D global epistasis (folding energy vs. folding fitness) - all phenotypes
  dataset_names <- names(pred_dt)[grepl("Abundance", names(pred_dt))]
  if(sum(grepl("_trait1", names(pred_dt)))==0){
    dataset_names <- unique(c(dataset_names, names(pred_dt)[grepl("Abundance|Binding", names(pred_dt))]))
  }
  for(i in dataset_names){
    mochims__plot_additive_trait_2state(
      mochi_outpath = mochi_outpath,
      input_dt = pred_dt[get(i)==1,], 
      RT = RT,
      trait_name = gsub("Abundance", "Folding", i),
      phenotype_name = i,
      report_outpath = outpath,
      dataset_name = gsub("Abundance|Binding", "", i),
      colour_scheme = colour_scheme)
  }

  #Plot 3D global epistasis (folding+binding energies vs. binding fitness) 
  dataset_names <- names(pred_dt)[grepl("Binding", names(pred_dt))]
  for(i in dataset_names){
    mochims__plot_additive_trait_3state(
      mochi_outpath = mochi_outpath,
      input_dt = pred_dt[get(i)==1,], 
      RT = RT,
      trait_names = c(gsub("Binding", "Folding", strsplit(i, "_")[[1]][1]), i),
      phenotype_name = i,
      report_outpath = outpath,
      dataset_name = gsub("Binding", "", i),
      colour_scheme = colour_scheme)
  }

  #Plot 3D global epistasis (sum-of-sigmoids, fit loess to fitted values for grid) 
  dataset_names <- names(pred_dt)[grepl("Fitness", names(pred_dt))]
  if(basename(mochi_outpath)=="tRNA-2dim"){
    for(i in dataset_names){
      mochims__plot_additive_trait_3state_sos(
        mochi_outpath = mochi_outpath,
        input_dt = pred_dt[get(i)==1,], 
        phenotype_name = i,
        report_outpath = outpath,
        dataset_name = gsub("Fitness", "", i),
        colour_scheme = colour_scheme)
    }
  }

  #Plot inferred global epistasis 
  dataset_names <- names(pred_dt)[grepl("PSI|Fitness|Brightness", names(pred_dt))]
  for(i in dataset_names){
    mochims__plot_additive_trait_inferred(
      mochi_outpath = mochi_outpath,
      input_dt = pred_dt[get(i)==1,], 
      phenotype_name = i,
      report_outpath = outpath,
      dataset_name = gsub("PSI|Fitness|Brightness", "", i),
      colour_scheme = colour_scheme)
  }

  # #Plot folding versus binding energy coloured by fraction bound showing isochores for arbitrary double mutant
  # mochims__plot_isochore_fraction_bound(
  #   mochi_outpath = mochi_outpath,
  #   input_dt = pred_dt_conf, 
  #   RT = RT,
  #   report_outpath = outpath,
  #   colour_scheme = colour_scheme)

  #Plot correlation with validation data
  mochims__plot_validation_scatter(
    input_dt = coef_dt, 
    lit_inpath = literature_free_energies, 
    report_outpath = outpath, 
    highlight_colour = colour_scheme[["shade 0"]][[1]],
    RT = RT)

  # #Add id with reference amino acid position
  # pred_dt_conf[, id_ref := mochims__get_reference_id(pred_dt_conf[,.(id, mut_order)], position_offset)]

  # #Add residue position for singles
  # pred_dt_conf[mut_order==1, Pos := as.integer(substr(id, 2, nchar(id)-1))]
  # pred_dt_conf[mut_order==1, Pos_ref := as.integer(substr(id_ref, 2, nchar(id_ref)-1))]

  #Save dGs and ddGs
  write.table(pred_dt, 
    file = file.path(outpath, "model_predictions.txt"), 
    quote = F, sep = "\t", row.names = F)
  write.table(coef_dt, 
    file = file.path(outpath, "model_coefficients.txt"), 
    quote = F, sep = "\t", row.names = F)

  # #Add ddG significance
  # pred_dt_conf[mut_order==1 & !duplicated(id), b_ddg_pred_sig := p.adjust(mochims__pvalue(b_ddg_pred, b_ddg_pred_sd), method = "BH")<0.05]
  # pred_dt_conf[mut_order==1 & !duplicated(id), f_ddg_pred_sig := p.adjust(mochims__pvalue(f_ddg_pred, f_ddg_pred_sd), method = "BH")<0.05]

  # #Save ids of singles with significant ddGs
  # write.table(pred_dt_conf[mut_order==1 & !duplicated(id) & b_ddg_pred_sig==T,id], 
  #   file = file.path(outpath, "b_ddg_sig_singles_id.txt"), 
  #   quote = F, sep = "\n", row.names = F, col.names = "id")
  # write.table(pred_dt_conf[mut_order==1 & !duplicated(id) & f_ddg_pred_sig==T,id], 
  #   file = file.path(outpath, "f_ddg_sig_singles_id.txt"), 
  #   quote = F, sep = "\n", row.names = F, col.names = "id")
}