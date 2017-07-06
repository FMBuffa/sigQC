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
#' @param origin Tells if datasets have come from different labs/experiments/machines. Is a vector of characters, with same character representing same origin. Default is assumption that all datasets come from the same source.
#' @param doNegativeControl Logical, tells the function if negative controls must be computed. TRUE by default.
#' @param numResampling Number of re-samplings (50 by default) during negative controls
#' @keywords make_all_plots
#' @export
#' @examples
# ' names = c("dataset1")
# ' data.matrix = replicate(10, rnorm(20))#random matrix - 10 genes x 20 samples
# ' mRNA_expr_matrix = list()
# ' mRNA_expr_matrix[["dataset1"]] = data.matrix
# ' row.names(mRNA_expr_matrix$dataset1) <- as.character(1:(dim(mRNA_expr_matrix$dataset1)[1]))
# ' colnames(mRNA_expr_matrix$dataset1) <- as.character(1:(dim(mRNA_expr_matrix$dataset1)[2]))
# ' #Define the signature
# ' gene_sigs_list = list()
# ' signature = "hypoxiaSig"
# ' gene_sig = c('1', '4', '5')#gene ids
# ' gene_sigs_list[[signature]] = as.matrix(gene_sig)
# ' names_sigs = c(signature)
# ' out_dir = file.path('~', "sigQC_Test")
# ' dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# ' make_all_plots(names_datasets = names, names_sigs = names_sigs, gene_sigs_list = gene_sigs_list, mRNA_expr_matrix = mRNA_expr_matrix, out_dir = out_dir)



make_all_plots <- function(gene_sigs_list, mRNA_expr_matrix,names_sigs=NULL,names_datasets=NULL , covariates=NULL, thresholds=NULL, out_dir = file.path('~', "sigQC"), showResults = FALSE,origin=NULL, doNegativeControl=TRUE, numResampling=50){

  #### print version number
  #### encoding scheme: major.minor
  #### major for large change
  #### minor for small change, whose results are expected to be similar as previous version. (two digits)
  print("-----sigQC Version 0.1.15-----")

  ###########Check the input
 radar_plot_values <- list();
  #check that there is a list of gene signatures
  if(missing(gene_sigs_list)){
    stop("Need to specify a list of gene signatures. The IDs must match those in the expression matrices.")
  }
  #check that there is a/are dataset(s)
  if(missing(mRNA_expr_matrix)){
    stop("Need to specify a list of expression matrices.
         Please note that each element of the list must have a
         name matching a value in the 'names' input parameter.")
  }

  #check that datasets have names, otherwise give them same names as in names_datasets
  if(is.null(names_datasets)){
    if(is.null(names(mRNA_expr_matrix))){
      stop("Need to specify a list of names. Each value must match a
           name in the mRNA_expr_matrix list of expression matrices.")
    } else {
      names_datasets <- names(mRNA_expr_matrix)
    }
  }
  #check that signatures have names, otherwise give them same names as in names_sigs
  if(is.null(names_sigs)){
    if(is.null(names(gene_sigs_list))){
      stop("Need to specify a list of names. Each value must match a
           name in the gene_sigs_list list of gene signatures.")
    } else {
      names_sigs <- names(gene_sigs_list)
    }
  }
  #set up the variable that will hold values for the final radarplot
  for(i in 1:length(names_sigs)){
    radar_plot_values[[names_sigs[i]]] <- list();
   }

   #need to check that the length of the origins variable is the same as the number of datasets
   if (!is.null(origin)){
    if (length(origin)!=length(names_datasets)){
      stop("Need to specify an origin for each dataset present, in the order the datasets are listed.")
    }
    if(sum(as.numeric(table(origin)) < 2) > 0){
      stop("Need to ensure that there are at least 2 datasets from each origin.")
    }
  }
  tryCatch({
    #Check if the OS is windows and allow longer path in this case
    os.name = Sys.info()["sysname"]
    os.release = Sys.info()["release"]
    useVirDrive = F
    if(doNegativeControl && nchar(out_dir)>100 && os.name=="Windows" && grepl(paste0(".*(7|8|10).*"), os.release, ignore.case=TRUE)){
      #out_dir = path.expand("~")
      # out_dir1 = file.path("/?", out_dir)
      #
      # out_dir1 = paste("\\?\\", out_dir, sep="")
      # out_dir2 = file.path("\\\\?\\C:\\Users\\Alessandro\\Documents")
      virD = "y:"
      system(paste0("subst ",virD, " ", out_dir), intern=T)
      out_dir = virD
      # outDir.file = file.path(outDir, paste0(paste(sample(letters, 251, TRUE), collapse = ''), '.txt'))
      # write(1, file = outDir.file)
      useVirDrive = T;
    }

       #check that the legnths of the names are all equal
    if ((length(names_datasets)== length(mRNA_expr_matrix)) && (length(names_sigs)==length(gene_sigs_list))){


      ###########Check if the needed package exists otherwise install it
      dir.create(out_dir)
     # utils::write.table('',file=file.path(out_dir, "log.log"))
      #LOG file path
      logfile.path = file.path(out_dir, "log.log")
      #Log conn
      log.con = file(logfile.path, open = "a") #open the logfile
      #run each of the sub functions
      cat(paste("LOG FILE CREATED: ",Sys.time(), sep=""), file=log.con, sep="\n") #start the log file

      #check that the signatures are of length >= 2, if not, we remove it
      sigs_to_remove_ind <- c()
      for (k in 1:length(names_sigs)){
        gene_sig <- gene_sigs_list[[names_sigs[k]]] #load in the gene signature
        if(length(gene_sig[,1]) < 2){
          sigs_to_remove_ind <- c(sigs_to_remove_ind,k)
  #        stop("Every signature must contain at least 2 elements.")
         }
      }
      #remove the offending signatures
      if (length(sigs_to_remove_ind) > 0){
        gene_sigs_list[names_sigs[k]] <- NULL
        for(j in 1:length(sigs_to_remove_ind)){
          cat(paste0("Need at least 2 elements per signature (signature removed): ",names_sigs[sigs_to_remove_ind[j]]), file=log.con, sep="\n")
        }
        names_sigs <- names_sigs[-sigs_to_remove_ind]#NULL
      }

      if(length(names_sigs) == 0){
        stop("No signatures contained at least 2 signature elements.") #if all signatures have been removed, then we have to stop
      }

      #first we must check that the genes of the signature are actually present in the datasets...
      datasets_to_remove_ind <- c()
      for (k in 1:length(names_sigs)){
        gene_sig <- gene_sigs_list[[names_sigs[k]]] #load in the gene signature
        for(i in 1:length(names_datasets)){
          data.matrix = mRNA_expr_matrix[[names_datasets[i]]] #load in the data matrix
          inter <- intersect(gene_sig[,1], row.names(data.matrix)) #make sure the genes are present in the dataset
          if(length(inter) < 2){
            #we need to exclude this dataset from further analysis because not enough signature genes are present
            datasets_to_remove_ind <- c(datasets_to_remove_ind,i)
          }
        }
      }
      datasets_to_remove_ind <- unique(datasets_to_remove_ind)
      if(length(datasets_to_remove_ind) > 0) {
        #now we do the removal
        mRNA_expr_matrix[names_datasets[datasets_to_remove_ind]] <- NULL
        #write it out that we've removed these datasets
        for(j in 1:length(datasets_to_remove_ind)){
          cat(paste0("Not enough signature genes present in dataset (dataset removed): ",names_datasets[datasets_to_remove_ind[j]]), file=log.con, sep="\n")
        }
        names_datasets <- names_datasets[-datasets_to_remove_ind]#NULL
      }

      if(length(names_datasets) == 0){
        stop("No datasets expressed at least 2 signature elements.") #if all datasets have been removed, then we have to stop
      }

      tryCatch(radar_plot_values <- eval_var_loc(gene_sigs_list,names_sigs, mRNA_expr_matrix,names_datasets,out_dir,file=log.con,showResults,radar_plot_values),
               error=function(err){
                 cat(paste0("Error occurred in eval_var_loc: ",err), file=log.con, sep="\n")
               })
      tryCatch(radar_plot_values <- eval_expr_loc(gene_sigs_list,names_sigs,mRNA_expr_matrix,names_datasets,thresholds, out_dir,file=log.con,showResults,radar_plot_values ),
               error=function(err){
                 cat(paste0("Error occurred in eval_expr_loc: ",err), file=log.con, sep="\n")
               })
      tryCatch(radar_plot_values <- eval_compactness_loc(gene_sigs_list,names_sigs,mRNA_expr_matrix,names_datasets,out_dir,file=log.con,showResults,radar_plot_values,logged=T,origin ),
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
      if(doNegativeControl){
        tryCatch(.sigsQcNegativeControl(genesList=gene_sigs_list, expressionMatrixList=mRNA_expr_matrix, outputDir=out_dir, logFile=logfile.path, numResampling=numResampling),
          error=function(err){
            cat(paste0("Error occurred during the execution of negative controls: ",err), file=log.con, sep="\n")
          }
        )
      }



      if(!showResults)
        grDevices::graphics.off() #make sure that everything is closed
      # close(log.con) #close connection to the file
    }else{
      print("Error: the length of names is not matching the number of elements in the expression matrices list.
            You need to have a name for every dataset and signature used.")
    }
  }, error = function(err) {
    #cat("", file=log.con, sep="\n")
    #cat(paste(Sys.time(),"Errors occurred during in sigQcNegativeControl:", err, sep=" "), file=log.con, sep="\n")
    #stop("Errors occurred during the computation of negative controls")
  }, finally = {
    #cat("---------------------------------------------------------------------------", file=log.con, sep="\n")
    if(useVirDrive){
      system(paste0("subst ",virD," /D"), intern =T)
    }
    close(log.con)
  })#END tryCatch
}
