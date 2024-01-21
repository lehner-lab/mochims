
#' mochims_coupling_scatter
#'
#' Coupling scatterplots.
#'
#' @param input_file path to MoCHI thermo model fit results (required)
#' @param validation_couplings validation coupling ids (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#'
#' @return Nothing
#' @export
#' @import data.table
mochims_coupling_scatter <- function(
  input_file,
  validation_couplings,
	outpath,
  colour_scheme
	){

  #Return if input_file doesn't exist
  if(!file.exists(input_file)){
    return()
  }

	#Load free energies
	dg_dt <- fread(input_file)

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

  #Display status
  message(paste("\n\n*******", paste0("running stage: mochims_coupling_scatter for ", domain_name), "*******\n\n"))

  #Create output directory
  mochims__create_dir(mochims_dir = outpath)

  #Weighted mean absolute value
  dg_dt[coef_order==2, mean_abs := abs(mean)]

  #GC compensatory mutations
  dg_dt[coef_order==2 & gsub("[0-9]", "", sapply(strsplit(id, "_"), '[', 1)) %in% c("ga", "gt", "ca", "ct"), dGC1 := TRUE]
  dg_dt[coef_order==2 & gsub("[0-9]", "", sapply(strsplit(id, "_"), '[', 2)) %in% c("ga", "gt", "ca", "ct"), dGC2 := TRUE]
  dg_dt[coef_order==2 & gsub("[0-9]", "", sapply(strsplit(id, "_"), '[', 1)) %in% c("ag", "tg", "ac", "tc"), iGC1 := TRUE]
  dg_dt[coef_order==2 & gsub("[0-9]", "", sapply(strsplit(id, "_"), '[', 2)) %in% c("ag", "tg", "ac", "tc"), iGC2 := TRUE]
  dg_dt[, gc_comp := FALSE]
  dg_dt[coef_order==2 & ((dGC1 & iGC2) | (iGC1 & dGC2)), gc_comp := TRUE]

  #Validation couplings
  dg_dt[coef_order==2, validation := id_ref %in% unlist(validation_couplings)]
  if(is.list(validation_couplings)){
    dg_dt[, validation := as.character(validation)]
    dg_dt[, validation := "remnd"]
    for(i in names(validation_couplings)){
      dg_dt[id_ref %in% validation_couplings[[i]], validation := i]
    }
  }

  #Plot colours
  plot_cols <- c('black', colour_scheme[["shade 0"]][1])
  names(plot_cols) <- c("FALSE", "TRUE")
  if(is.list(validation_couplings)){
    plot_cols <- c(colour_scheme[["shade 0"]][2], colour_scheme[["shade 0"]][4], 'black')
    names(plot_cols) <- c(names(validation_couplings), "remnd")
  }

  #Coupling scatter - backbone distance
  plot_dt <- dg_dt[coef_order==2][,.(mean_abs, ci95, backbone, validation, id_ref)]
  plot_dt[is.na(mean_abs), mean_abs := 0]
  plot_dt[is.na(ci95), ci95 := 0]
  plot_dt <- plot_dt[order(mean_abs, decreasing = T)]
  plot_dt[, mean_abs_rank := 1:.N]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mean_abs_rank, mean_abs, color = validation)) +
    ggplot2::geom_point() +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = mean_abs-ci95/2, ymax = mean_abs+ci95/2), alpha = 1/4) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) + 
    ggplot2::xlab("Rank") +
    ggplot2::ylab("|coupling term|") +
    ggrepel::geom_text_repel(ggplot2::aes(label = id_ref), show.legend = F, 
                             max.overlaps = 20, xlim = c(-10, 10), ylim = c(-10, 10), size = 1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_manual(values = plot_cols)
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_meanabsrank.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

	#Coupling scatter - backbone distance
	plot_dt <- dg_dt[coef_order==2][,.(mean_abs, ci95, backbone, validation)]
	plot_dt[is.na(mean_abs), mean_abs := 0]
	plot_dt[is.na(ci95), ci95 := 0]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(backbone, mean_abs, color = validation)) +
    ggplot2::geom_point() +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = mean_abs-ci95/2, ymax = mean_abs+ci95/2), alpha = 1/4) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_hline(yintercept = 0, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(mean_abs, backbone, use = "pairwise.complete", method = 'spearman'), 2), sep=""), validation)], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Backbone distance (residues)") +
    ggplot2::ylab("|coupling term|") +
    ggplot2::theme_classic() +
		ggplot2::scale_colour_manual(values = plot_cols)
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_backbonedist.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #3D distance metric exists
  if("scHAmin" %in% names(dg_dt)){

    #Coupling scatter - contact distance
    plot_dt <- dg_dt[coef_order==2][,.(mean_abs, ci95, scHAmin, validation)]
    plot_dt[is.na(mean_abs), mean_abs := 0]
    plot_dt[is.na(ci95), ci95 := 0]
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin, mean_abs, color = validation)) +
      ggplot2::geom_point() +
      ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = mean_abs-ci95/2, ymax = mean_abs+ci95/2), alpha = 1/4) +
      ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
      ggplot2::geom_hline(yintercept = 0, linetype = 2) + 
      ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(mean_abs, scHAmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""), validation)], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
      ggplot2::xlab("Contact distance (Angstrom)") +
      ggplot2::ylab("|coupling term|") +
      ggplot2::theme_classic() +
      ggplot2::scale_colour_manual(values = plot_cols)
    ggplot2::ggsave(file.path(outpath, "coupling_scatter_scHAmin.pdf"), d, width = 4, height = 3, useDingbats=FALSE)


  }

  #Multiple trait names
  tnames <- dg_dt[,unique(trait_name)]
  if(length(tnames)==2){

    dg_dt1 <- dg_dt[trait_name==tnames[1]]
    dg_dt2 <- dg_dt[trait_name==tnames[2]]
    plot_dt <- merge(
      dg_dt1[coef_order!=0,.SD,,.SDcols = c("id", "coef_order", "mean", "ci95", "validation", "backbone", "gc_comp")],
      dg_dt2[coef_order!=0,.SD,,.SDcols = c("id", "coef_order", "mean", "ci95", "validation", "backbone", "gc_comp")],
      by = c('id', "validation", "backbone", "coef_order", "gc_comp"), suffixes = c("_1", "_2"))

    #Contact + bbgcc
    plot_dt[coef_order==1, annotation := "order1"]
    plot_dt[coef_order==2, annotation := "remnd"]
    plot_dt[coef_order==2 & id %in% unlist(validation_couplings), annotation := "contc"]
    plot_dt[coef_order==2 & backbone<=5 & gc_comp, annotation := "bbgcc"]
    #Plot colours
    plot_cols <- c('black', 'grey', colour_scheme[["shade 0"]][1], colour_scheme[["shade 0"]][2])
    names(plot_cols) <- c("order1", "remnd", "contc", "bbgcc")
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mean_1, mean_2, color = annotation)) +
      ggplot2::geom_point(data = plot_dt[annotation == "order1" | annotation=='remnd'], size = 1) +
      ggplot2::geom_linerange(data = plot_dt[annotation == "order1" | annotation=='remnd'], ggplot2::aes(xmin = mean_1-ci95_1/2, xmax = mean_1+ci95_1/2), alpha = 1/4) +
      ggplot2::geom_linerange(data = plot_dt[annotation == "order1" | annotation=='remnd'], ggplot2::aes(ymin = mean_2-ci95_2/2, ymax = mean_2+ci95_2/2), alpha = 1/4) +
      ggplot2::geom_point(data = plot_dt[!(annotation == "order1" | annotation=='remnd')]) +
      ggplot2::geom_linerange(data = plot_dt[!(annotation == "order1" | annotation=='remnd')], ggplot2::aes(xmin = mean_1-ci95_1/2, xmax = mean_1+ci95_1/2), alpha = 1/4) +
      ggplot2::geom_linerange(data = plot_dt[!(annotation == "order1" | annotation=='remnd')], ggplot2::aes(ymin = mean_2-ci95_2/2, ymax = mean_2+ci95_2/2), alpha = 1/4) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) + 
      ggplot2::geom_hline(yintercept = 0, linetype = 2) + 
      ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(mean_1, mean_2, use = "pairwise.complete", method = 'spearman'), 2), sep=""), annotation, coef_order)], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
      ggplot2::xlab("Additive trait 1") +
      ggplot2::ylab("Additive trait 2") +
      ggplot2::theme_classic() +
      ggplot2::scale_colour_manual(values = plot_cols)
    ggplot2::ggsave(file.path(outpath, "coupling_scatter_2dim_contact_bbgcc.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

    #Contact + bb
    plot_dt[coef_order==1, annotation := "order1"]
    plot_dt[coef_order==2, annotation := "remnd"]
    plot_dt[coef_order==2 & id %in% unlist(validation_couplings), annotation := "contc"]
    plot_dt[coef_order==2 & backbone<=5, annotation := "bb"]
    #Plot colours
    plot_cols <- c('black', 'grey', colour_scheme[["shade 0"]][1], colour_scheme[["shade 0"]][2])
    names(plot_cols) <- c("order1", "remnd", "contc", "bb")
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mean_1, mean_2, color = annotation)) +
      ggplot2::geom_point(data = plot_dt[annotation == "order1" | annotation=='remnd'], size = 1) +
      ggplot2::geom_linerange(data = plot_dt[annotation == "order1" | annotation=='remnd'], ggplot2::aes(xmin = mean_1-ci95_1/2, xmax = mean_1+ci95_1/2), alpha = 1/4) +
      ggplot2::geom_linerange(data = plot_dt[annotation == "order1" | annotation=='remnd'], ggplot2::aes(ymin = mean_2-ci95_2/2, ymax = mean_2+ci95_2/2), alpha = 1/4) +
      ggplot2::geom_point(data = plot_dt[!(annotation == "order1" | annotation=='remnd')]) +
      ggplot2::geom_linerange(data = plot_dt[!(annotation == "order1" | annotation=='remnd')], ggplot2::aes(xmin = mean_1-ci95_1/2, xmax = mean_1+ci95_1/2), alpha = 1/4) +
      ggplot2::geom_linerange(data = plot_dt[!(annotation == "order1" | annotation=='remnd')], ggplot2::aes(ymin = mean_2-ci95_2/2, ymax = mean_2+ci95_2/2), alpha = 1/4) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) + 
      ggplot2::geom_hline(yintercept = 0, linetype = 2) + 
      ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(mean_1, mean_2, use = "pairwise.complete", method = 'spearman'), 2), sep=""), annotation, coef_order)], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
      ggplot2::xlab("Additive trait 1") +
      ggplot2::ylab("Additive trait 2") +
      ggplot2::theme_classic() +
      ggplot2::scale_colour_manual(values = plot_cols)
    ggplot2::ggsave(file.path(outpath, "coupling_scatter_2dim_contact_bb.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

    #Contact + gcc
    plot_dt[coef_order==1, annotation := "order1"]
    plot_dt[coef_order==2, annotation := "remnd"]
    plot_dt[coef_order==2 & id %in% unlist(validation_couplings), annotation := "contc"]
    plot_dt[coef_order==2 & gc_comp, annotation := "gcc"]
    #Plot colours
    plot_cols <- c('black', 'grey', colour_scheme[["shade 0"]][1], colour_scheme[["shade 0"]][2])
    names(plot_cols) <- c("order1", "remnd", "contc", "gcc")
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mean_1, mean_2, color = annotation)) +
      ggplot2::geom_point(data = plot_dt[annotation == "order1" | annotation=='remnd'], size = 1) +
      ggplot2::geom_linerange(data = plot_dt[annotation == "order1" | annotation=='remnd'], ggplot2::aes(xmin = mean_1-ci95_1/2, xmax = mean_1+ci95_1/2), alpha = 1/4) +
      ggplot2::geom_linerange(data = plot_dt[annotation == "order1" | annotation=='remnd'], ggplot2::aes(ymin = mean_2-ci95_2/2, ymax = mean_2+ci95_2/2), alpha = 1/4) +
      ggplot2::geom_point(data = plot_dt[!(annotation == "order1" | annotation=='remnd')]) +
      ggplot2::geom_linerange(data = plot_dt[!(annotation == "order1" | annotation=='remnd')], ggplot2::aes(xmin = mean_1-ci95_1/2, xmax = mean_1+ci95_1/2), alpha = 1/4) +
      ggplot2::geom_linerange(data = plot_dt[!(annotation == "order1" | annotation=='remnd')], ggplot2::aes(ymin = mean_2-ci95_2/2, ymax = mean_2+ci95_2/2), alpha = 1/4) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) + 
      ggplot2::geom_hline(yintercept = 0, linetype = 2) + 
      ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(mean_1, mean_2, use = "pairwise.complete", method = 'spearman'), 2), sep=""), annotation, coef_order)], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
      ggplot2::xlab("Additive trait 1") +
      ggplot2::ylab("Additive trait 2") +
      ggplot2::theme_classic() +
      ggplot2::scale_colour_manual(values = plot_cols)
    ggplot2::ggsave(file.path(outpath, "coupling_scatter_2dim_contact_gcc.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  }

  return(dg_dt)
}
