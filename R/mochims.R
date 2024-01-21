
#' mochims
#'
#' Main analysis script.
#'
#' @param startStage Start at a specified analysis stage (default:1)
#' @param stopStage Stop at a specified analysis stage (default:5)
#' @param base_dir Base directory for all input files (default:NB private CRG server path; change accordingly)
#' @param output_dir Output directory for all output files (default:same as base_dir)
#'
#' @return Nothing
#' @export
mochims <- function(
	startStage = 1,
	stopStage = 5,
	base_dir = "/users/project/prj004631/afaure/DMS/Results/mochims_proj",
	output_dir = ""
	){

	colour_scheme <- list(
		"shade 0" = list(
			"#F4270C",
			"#F4AD0C",
			"#1B38A6",
			"#09B636"),
		"shade 1" = list(
			"#FFB0A5",
			"#FFE4A5",
			"#9DACE3",
			"#97E9AD"),
		"shade 2" = list(
			"#FF6A56",
			"#FFCB56",
			"#4C63B7",
			"#43C766"),
		"shade 3" = list(
			"#A31300",
			"#A37200",
			"#0C226F",
			"#007A20"),
		"shade 4" = list(
			"#410800",
			"#412D00",
			"#020B2C",
			"#00300D"),
		"shade 5" = list(
			"#CCCCCC",
			"#FF991F",
			"#5CB8FF",
			"#B22222"
		))

	#First and last analysis stages
	first_stage <- startStage
	last_stage <- stopStage

	#Default output directory
	if(output_dir==""){
		output_dir <- base_dir
	}

	### Evaluate thermo model results
	###########################

	stagenum <- 1
	mochims_thermo_model_results_datasets(
		dataset_names = c(
			"FOSJUN",
			"GB1",
			"PSD95-PDZ3",
			"KRAS",
			"FAS",
			"FAS-mechanistic",
			"FAS-linear",
			"tRNA",
			"tRNA-o3",
			"tRNA-2dim",
			"tRNA-3dim",
			"tRNA-4dim",
			"tRNA-5dim",
			"tRNA-linear",
			"eqFP611-ensemble",
			"eqFP611-ensemble-top100",
			"eqFP611-ensemble-1000vars",
			"eqFP611-ensemble-1000vars-top100",
			"eqFP611-o2",
			"eqFP611-o2-1000vars"),
		literature_free_energies = file.path(base_dir, "Data", "in_vitro"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Add 3D structure metrics
	###########################

	stagenum <- 2
	mochims_structure_metrics_datasets(
		dataset_names = c(
			"FAS",
			"FAS-mechanistic",
			"FAS-linear",
			"tRNA",
			"tRNA-2dim",
			"tRNA-linear",
			"eqFP611-o2"),
		pdb_file_list = list(
			"eqFP611-o2" = file.path(base_dir, "Data", "pdb", "3m24_edit.pdb")),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Position class violin plots
	###########################

	stagenum <- 3
	mochims_position_class_violins_datasets(
		dataset_names = c(
			"FOSJUN"),
		annotations = file.path(base_dir, "Data", "annotations"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Coefficient size distribution plots
	###########################

	stagenum <- 4
	mochims_coefficient_size_distribution_datasets(
		dataset_names = c(
			"eqFP611-ensemble",
			"eqFP611-ensemble-top100",
			"eqFP611-ensemble-1000vars",
			"eqFP611-ensemble-1000vars-top100",
			"eqFP611-o2",
			"eqFP611-o2-1000vars"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Coupling scatterplots
	###########################

	stagenum <- 5
	mochims_coupling_scatter_datasets(
		dataset_names = c(
			"FAS",
			"FAS-mechanistic",
			"FAS-linear",
			"tRNA",
			"tRNA-2dim",
			"tRNA-linear",
			"eqFP611-o2"),
		validation_list = list(
			"FAS" = c("c18t_t19g", "t49c_g51c", "t24c_g26t", "c18g_t19g", "c39t_g44a", "c32t_g35t", "c39t_c41g"),
			"FAS-mechanistic" = c("c18t_t19g", "t49c_g51c", "t24c_g26t", "c18g_t19g", "c39t_g44a", "c32t_g35t", "c39t_c41g"),
			"FAS-linear" = c("c18t_t19g", "t49c_g51c", "t24c_g26t", "c18g_t19g", "c39t_g44a", "c32t_g35t", "c39t_c41g"),
			"tRNA" = c("g1a_c71t", "t2c_a70g", "t2g_a70t", "g6a_c66t", "g6t_c66a"),
			"tRNA-2dim" = c("g1a_c71t", "t2c_a70g", "t2g_a70t", "g6a_c66t", "g6t_c66a"),
			"tRNA-linear" = c("g1a_c71t", "t2c_a70g", "t2g_a70t", "g6a_c66t", "g6t_c66a"),
			"eqFP611-o2" = list(
				chrom = c("L63M_N158A", "L63M_Y197R", "L63M_A174L", "L63M_F143S"), 
				nbrhd = c("N158A_Y197R", "N158A_A174L", "F143S_N158A", "A174L_Y197R", "F143S_Y197R", "F143S_A174L"))),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

}

