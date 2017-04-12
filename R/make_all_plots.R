#' make_all_plots.R
#'
#' Makes all the plots for the quality control of the list(s) of genes.
#'
#' @param gene_sigs_list A list of genes representing the gene signature to be tested.
#' @param mRNA_expr_matrix A list of expression matrices
#' @param names_sigs The names of the gene signatures (e.g. Hypoxia, Invasiveness), one name per each signature in gene_sigs_list.
#' @param names_datasets The names of the different datasets contained in mRNA_expr_matrix
#' @param covariates A list containing a sub-list of 'annotations' and 'colors' which contains the annotation matrix for the given dataset and the associated colours with which to plot in the expression heatmap
#' @param thresholds A list of thresholds to be considered for each data set, default is median of the data set. A gene is considered expressed if above the threshold, non-expressed otherwise. One threshold per dataset, in the same order as the dataset list.
#' @param out_dir A path to the directory where the resulting output files are written
#' @param showResults Tells if open dialog boxes showing the computed results. Default is FALSE
#' @keywords make_all_plots
#' @export
#' @examples
# ' names = c("dataset1")
# ' data.matrix = replicate(10, rnorm(20))#random matrix - 10 genes x 20 samples
# ' mRNA_expr_matrix = list()
# ' mRNA_expr_matrix["dataset1"] = list(data.matrix)
# ' row.names(mRNA_expr_matrix$dataset1) <- as.character(1:(dim(mRNA_expr_matrix$dataset1)[1]))
# ' colnames(mRNA_expr_matrix$dataset1) <- as.character(1:(dim(mRNA_expr_matrix$dataset1)[2]))
# ' gene_sig = c('1', '4', '5')#gene ids
# ' make_all_plots(gene_sig, mRNA_expr_matrix, names)



make_all_plots <- function(gene_sigs_list, mRNA_expr_matrix, names_sigs=NULL,names_datasets=NULL , covariates=NULL, thresholds=NULL, out_dir = '~', showResults = FALSE){
  ###########Check the input
 radar_plot_values <- list();

  if(missing(gene_sigs_list)){
    stop("Need to specify a list of gene signatures. The IDs must match those in the expression matrices.")
  }
  if(missing(mRNA_expr_matrix)){
    stop("Need to specify a list of expression matrices.
         Please note that each element of the list must have a
         name matching a value in the 'names' input parameter.")
  }
  if(is.null(names_datasets)){
    if(is.null(names(mRNA_expr_matrix))){
      stop("Need to specify a list of names. Each value must match a
           name in the mRNA_expr_matrix list of expression matrices.")
    } else {
      names_datasets <- names(mRNA_expr_matrix)
    }
  }
  if(is.null(names_sigs)){
    if(is.null(names(gene_sigs_list))){
      stop("Need to specify a list of names. Each value must match a
           name in the gene_sigs_list list of gene signatures.")
    } else {
      names_sigs <- names(gene_sigs_list)
    }
  }
  for(i in 1:length(names_sigs)){
    radar_plot_values[[names_sigs[i]]] <- list();
   }
  if ((length(names_datasets)== length(mRNA_expr_matrix)) && (length(names_sigs)==length(gene_sigs_list))){
    ###########Check if the needed package exists otherwise install it
    dir.create(out_dir)
   # utils::write.table('',file=file.path(out_dir, "log.log"))
    #LOG file path
    logfile.path = file.path(out_dir, "log.log")
    #Log conn
    log.con = file(logfile.path, open = "a")
    cat(paste("LOG FILE CREATED: ",Sys.time(), sep=""), file=log.con, sep="\n")
    tryCatch(radar_plot_values <- eval_var_loc(gene_sigs_list,names_sigs, mRNA_expr_matrix,names_datasets,out_dir,file=log.con,showResults,radar_plot_values),
             error=function(err){
               cat(paste0("Error occurred in eval_var_loc: ",err), file=log.con, sep="\n")
             })
    tryCatch(radar_plot_values <- eval_expr_loc(gene_sigs_list,names_sigs,mRNA_expr_matrix,names_datasets,thresholds, out_dir,file=log.con,showResults,radar_plot_values ),
             error=function(err){
               cat(paste0("Error occurred in eval_expr_loc: ",err), file=log.con, sep="\n")
             })
    tryCatch(radar_plot_values <- eval_compactness_loc(gene_sigs_list,names_sigs,mRNA_expr_matrix,names_datasets,out_dir,file=log.con,showResults,radar_plot_values ),
             error=function(err){
               cat(paste0("Error occurred during the evaluation of compactness: ",err), file=log.con, sep="\n")
             })

    tryCatch({radar_plot_values <- compare_metrics_loc(gene_sigs_list,names_sigs,mRNA_expr_matrix,names_datasets,out_dir,file=log.con,showResults,radar_plot_values )},
             error=function(err){
            #  print(paste0("Error, likely due to inability to calculate PCA, because of missing values: ", err))
               cat(paste0("Error occurred, likely due to inability to calculate PCA, because of missing values:  ",err), file=log.con, sep="\n")
             })
    tryCatch({radar_plot_values <- eval_stan_loc(gene_sigs_list,names_sigs,mRNA_expr_matrix,names_datasets,out_dir,file=log.con,showResults,radar_plot_values )},
          error=function(err){
           cat(paste0("Error occurred in eval_stan_loc: ",err), file=log.con, sep="\n")
         })
    tryCatch(radar_plot_values <- eval_struct_loc(gene_sigs_list,names_sigs,mRNA_expr_matrix,names_datasets,covariates,out_dir,file=log.con,showResults,radar_plot_values ),
             error=function(err){
               cat(paste0("Error occurred in eval_struct_loc: ",err), file=log.con, sep="\n")
             })
    # print(radar_plot_values)
    tryCatch(make_radar_chart_loc(radar_plot_values,showResults,names_sigs, names_datasets,out_dir,file=log.con),
             error=function(err){
               cat(paste0("Error occurred in make_radar_chart_loc: ",err), file=log.con, sep="\n")
             })

    if(!showResults)
      grDevices::graphics.off()
    close(log.con)
  }else{
    print("Error: the length of names is not matching the number of elements in the expression matrices list.
          You need to have a name for every dataset used.")
  }
}
