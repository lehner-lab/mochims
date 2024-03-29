
#' mochims__plot_validation_scatter
#'
#' Plot correlation with validation data.
#'
#' @param input_dt data.table with model free energy estimates (required)
#' @param lit_inpath path to literature free energies (required)
#' @param report_outpath output path for scatterplots (required)
#' @param highlight_colour colour for highlights (default:red)
#' @param RT constant (default:0.001987*(273+24))
#' @param position_offset residue position offset (default:0)
#'
#' @return Nothing
#' @export
#' @import data.table
mochims__plot_validation_scatter <- function(
  input_dt,
  lit_inpath,
  report_outpath,
  highlight_colour = "red",
  RT = 0.001987*(273+24)
  ){

  #Check validation data exists
  if(!file.exists(lit_inpath)){return()}

  lit_dt <- fread(lit_inpath)
  if("Assay/Protocol" %in% names(lit_dt)){
    #ProtaBank format
    lit_dt[, id := Description]
    names(lit_dt) <- gsub("/Protocol", "", names(lit_dt))
    #Retain single mutants only
    lit_dt <- lit_dt[!grepl(",", id)]
    #Add dataset
    lit_dt[, Dataset := gsub("dG", "", gsub("ddG", "", gsub("SD of ", "", Assay)))]
    #Split
    lit_dt_dG <- lit_dt[grepl("dG", Assay) & !grepl("ddG", Assay) & !grepl("SD of ", Assay),]
    lit_dt_ddG <- lit_dt[grepl("ddG", Assay) & !grepl("SD of ", Assay),]
    lit_dt_dG_SD <- lit_dt[grepl("dG", Assay) & !grepl("ddG", Assay) & grepl("SD of ", Assay),]
    lit_dt_ddG_SD <- lit_dt[grepl("ddG", Assay) & grepl("SD of ", Assay),]
    #Merge
    lit_dt_dG <- merge(lit_dt_dG, lit_dt_dG_SD[, .(id, Dataset, Data_sd = Data)], by = c("id", "Dataset"), all = T)
    lit_dt_ddG <- merge(lit_dt_ddG, lit_dt_ddG_SD[, .(id, Dataset, Data_sd = Data)], by = c("id", "Dataset"), all = T)
    lit_dt <- merge(
      lit_dt_dG[,.(id, Dataset, f_dg = -Data, f_dg_sd = Data_sd)],
      lit_dt_ddG[,.(id, Dataset, f_ddg = -Data, f_ddg_sd = Data_sd)], by = c("id", "Dataset"), all = T)
    lit_dt[, b_dg := NA]
    lit_dt[, b_ddg := NA]
    lit_dt[, b_dg_sd := NA]
    lit_dt[, b_ddg_sd := NA]
  }else{
    lit_dt <- fread(lit_inpath, colClasses = list("character" = 1:2, "numeric" = 3:13))
    lit_dt[!is.na(Kd), b_dg := log(Kd)*RT]
    lit_dt[!is.na(Kf) & !is.na(Ku), f_dg := -log(Kf/Ku)*RT]
    # #Logical columns to numeric
    # to.replace <- names(which(sapply(lit_dt, is.logical)))
    # for(var in to.replace){lit_dt[, (var) := as.numeric(get(var))]}
    for(i in lit_dt[,unique(Dataset)]){
      if(lit_dt[Dataset==i & id=="WT",.N]!=0){
        lit_dt[Dataset==i & !is.na(b_dg), b_ddg := b_dg-lit_dt[Dataset==i & id=="WT",b_dg]]
        lit_dt[Dataset==i & !is.na(f_dg), f_ddg := f_dg-lit_dt[Dataset==i & id=="WT",f_dg]]
      }
    }
  }

  #Remove WT
  lit_dt <- lit_dt[id!="WT"]
  #Convert to data ids
  lit_dt[, mut_pos := as.integer(substr(id, 2, nchar(id)-1))]
  lit_dt[, wtAA := substr(id, 1, 1)]
  lit_dt[, mutAA := substr(id, nchar(id), nchar(id))]
  lit_dt[, id_ref := paste0(wtAA, mut_pos, mutAA)]
  #Retain confident literature dGs and ddGs
  lit_dt <- lit_dt[f_dg_sd>=1/(1.96*2) | f_dg_sd==0, f_dg := NA]
  lit_dt <- lit_dt[f_ddg_sd>=1/(1.96*2) | f_ddg_sd==0, f_ddg := NA]
  lit_dt <- lit_dt[b_dg_sd>=1/(1.96*2) | b_dg_sd==0, b_dg := NA]
  lit_dt <- lit_dt[b_ddg_sd>=1/(1.96*2) | b_ddg_sd==0, b_ddg := NA]

  for(dname in lit_dt[,unique(Dataset)]){
    for(conf_level in c(T, F)){
      for(metric_lit in c("f_dg", "f_ddg", "b_dg", "b_ddg")){
        #Binding partner
        binding_partner <- ""
        if("Partner" %in% names(lit_dt)){
          binding_partner <- lit_dt[Dataset==dname,unique(Partner)][1]
        }
        #Merge with model results
        plot_dt <- merge(lit_dt[Dataset==dname], input_dt, by = "id_ref", allow.cartesian=TRUE)
        if(plot_dt[!is.na(get(metric_lit)),.N]==0){
          next()
        }
        metric_lit_sd <- paste0(metric_lit, "_sd")
        metric_data <- "mean_kcal/mol"
        metric_data_sd <- "std_kcal/mol"
        metric_data_conf <- "conf"
        plot_dt[, col_lit := .SD,,.SDcols = metric_lit]
        plot_dt[, col_lit_sd := .SD,,.SDcols = metric_lit_sd]
        plot_dt[, col_data := .SD,,.SDcols = metric_data]
        plot_dt[, col_data_sd := .SD,,.SDcols = metric_data_sd]
        plot_dt[, col_data_conf := .SD,,.SDcols = metric_data_conf]
        if(conf_level){plot_dt <- plot_dt[col_data_conf==T]}
        #Remove missing data
        plot_dt <- plot_dt[!is.na(col_lit) & !is.na(col_data)]
        #Subset to trait of interest
        if(metric_lit %in% c("f_dg", "f_ddg")){
          plot_dt <- plot_dt[grepl('Folding', trait_name)]
        }
        if(metric_lit %in% c("b_dg", "b_ddg")){
          if(binding_partner!=""){
            plot_dt <- plot_dt[grepl(paste0(c("Binding", binding_partner), collapse = "_"), trait_name)]
          }else{
            plot_dt <- plot_dt[grepl("Binding", trait_name)]
          }
        }        
        #Plot
        d <- ggplot2::ggplot(plot_dt,ggplot2::aes(col_lit, col_data)) +
          ggplot2::geom_hline(yintercept = 0) +
          ggplot2::geom_vline(xintercept = 0) +
          ggplot2::geom_point(alpha = 1/2, color = highlight_colour, size = 1) +
          ggplot2::geom_linerange(data = plot_dt[!is.na(col_data_sd)], ggplot2::aes(ymin = col_data-col_data_sd*1.96, ymax = col_data+col_data_sd*1.96), color = highlight_colour, alpha = 1/4) +
          ggplot2::geom_linerange(data = plot_dt[!is.na(col_lit_sd)], ggplot2::aes(xmin = col_lit-col_lit_sd*1.96, xmax = col_lit+col_lit_sd*1.96), color = highlight_colour, alpha = 1/4) +
          ggplot2::geom_smooth(formula = 'y ~ x', linetype = 1, method = "lm", color = highlight_colour, se = T) +
          ggplot2::geom_abline(linetype = 2) +
          ggplot2::annotate("text", label=paste("Pearson's r = ", plot_dt[,round(cor(col_lit, col_data), 2)], sep=""), x=-Inf, y=Inf, hjust = 0, vjust = 1) +
          ggplot2::theme_classic()
        if(grepl("^f_", metric_lit) & grepl("_ddg$", metric_lit)){
          d <- d + ggplot2::xlab(expression(Delta*Delta*"G Folding (in vitro)")) +
            ggplot2::ylab(expression(Delta*Delta*"G Folding (inferred)"))
        }else if(grepl("^b_", metric_lit) & grepl("_ddg$", metric_lit)){
          d <- d + ggplot2::xlab(expression(Delta*Delta*"G Binding (in vitro)")) +
            ggplot2::ylab(expression(Delta*Delta*"G Binding (inferred)"))
        }else if(grepl("^f_", metric_lit) & grepl("_dg$", metric_lit)){
          d <- d + ggplot2::xlab(expression(Delta*"G Folding (in vitro)")) +
            ggplot2::ylab(expression(Delta*"G Folding (inferred)"))
        }else if(grepl("^b_", metric_lit) & grepl("_dg$", metric_lit)){
          d <- d + ggplot2::xlab(expression(Delta*"G Binding (in vitro)")) +
            ggplot2::ylab(expression(Delta*"G Binding (inferred)"))
        }
        ggplot2::ggsave(file.path(report_outpath, paste0("validation_scatter_", paste0(c(dname, binding_partner), collapse = "_"), "_", metric_lit, c("", "_conf")[as.numeric(conf_level)+1], ".pdf")), d, width = 3, height = 3, useDingbats=FALSE)
        #Save plot data
        write.table(plot_dt[,.(id_ref, col_lit, col_data)], file = file.path(report_outpath, paste0("validation_scatter_", paste0(c(dname, binding_partner), collapse = "_"), "_", metric_lit, c("", "_conf")[as.numeric(conf_level)+1], ".txt")), quote = F, row.names = F)
        #Axis ranges at least 3 kcal/mol  
        if(dname=="Malagrino_2019"){
          temp_ylim_max <- max(c(plot_dt[,mean(col_data)] + 1.5, plot_dt[!is.na(col_data_sd),unlist(col_data+col_data_sd*1.96)]))
          temp_ylim_max <- max(c(temp_ylim_max, plot_dt[is.na(col_data_sd),unlist(col_data)]))
          temp_ylim_min <- min(c(plot_dt[,mean(col_data)] - 1.5, plot_dt[!is.na(col_data_sd),unlist(col_data-col_data_sd*1.96)]))
          temp_ylim_min <- min(c(temp_ylim_min, plot_dt[is.na(col_data_sd),unlist(col_data)]))
          temp_xlim_max <- max(c(plot_dt[,mean(col_lit)] + 1.5, plot_dt[!is.na(col_lit_sd),unlist(col_lit+col_lit_sd*1.96)]))
          temp_xlim_max <- max(c(temp_xlim_max, plot_dt[is.na(col_lit_sd),unlist(col_lit)]))
          temp_xlim_min <- min(c(plot_dt[,mean(col_lit)] - 1.5, plot_dt[!is.na(col_lit_sd),unlist(col_lit-col_lit_sd*1.96)]))
          temp_xlim_min <- min(c(temp_xlim_min, plot_dt[is.na(col_lit_sd),unlist(col_lit)]))
          d <- d + 
            ggplot2::scale_x_continuous(limits = c(temp_xlim_min, temp_xlim_max)) +
            ggplot2::scale_y_continuous(limits = c(temp_ylim_min, temp_ylim_max))
          ggplot2::ggsave(file.path(report_outpath, paste0("validation_scatter_", paste0(c(dname, binding_partner), collapse = "_"), "_", metric_lit, c("", "_conf")[as.numeric(conf_level)+1], "_minrange.pdf")), d, width = 3, height = 3, useDingbats=FALSE)
        }
      }
    }
  }
}
