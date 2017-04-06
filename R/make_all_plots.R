#' make_all_plots.R
#'
#' Makes the plots in each of the subfunctions.
#' @param gene_sig A list of genes representing the gene signature to be tested.
#' @param mRNA_expr_matrix A list of expression matrices
#' @param names The names of the different datasets contained in mRNA_expr_matrix
#' @param covariates A list containing a sub-list of 'annotations' and 'colors' which contains the annotation matrix for the given dataset and the associated colours with which to plot in the expression heatmap
#' @param thresholds A list of thresholds to be considered for each data set, default is median of the data set. A gene is considered expressed if above the threshold, non-expressed otherwise.
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
      stop("Neeed specify a list of names. Each value must match a
           name in the mRNA_expr_matrix list of expression matrices.")
    } else {
      names_datasets <- names(mRNA_expr_matrix)
    }
  }
  if(is.null(names_sigs)){
    if(is.null(names(gene_sigs_list))){
      stop("Neeed specify a list of names. Each value must match a
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
    source("http://bioconductor.org/biocLite.R")
    #ADDED TO RUN THE EXAMPLE WITHOUT ERRORS
    if(!require("MASS")){
      biocLite("MASS", dependencies=TRUE)
      library(MASS)
    }
    if(!require("lattice")){
      biocLite("lattice", dependencies=TRUE)
      library(lattice)
    }
    if(!require("KernSmooth")){
      biocLite("KernSmooth", dependencies=TRUE)
      library(KernSmooth)
    }
    if(!require("cluster")){
      biocLite("cluster", dependencies=TRUE)
      library(cluster)
    }
    if(!require("nnet")){
      biocLite("nnet", dependencies=TRUE)
      library(nnet)
    }
    if(!require("class")){
      biocLite("class", dependencies=TRUE)
      library(class)
    }
    #REQUIRED

    if(!require("gridGraphics")){
      biocLite("gridGraphics", dependencies=TRUE)
      library(gridGraphics)
    }
    if(!require("biclust")){
      biocLite(pkgs = "biclust", dependencies=TRUE)
      library(biclust)
    }
    if(!require("gplots")){
      biocLite("gplots", dependencies=TRUE)
      library(gplots)
    }
    # if(!require("grid")){
    #   biocLite("grid", dependencies=TRUE)
    #   library(grid)
    # }
    if(!require("ComplexHeatmap")){
      biocLite("ComplexHeatmap", dependencies=TRUE)
      library(ComplexHeatmap)
    }
    if(!require("RankProd")){
      biocLite("RankProd", dependencies=TRUE)
      library(RankProd)
    }
    if(!require("fmsb")){
      install.packages("fmsb")
      library(fmsb)
    }
    if(!require("moments")){
      install.packages("moments")
      library(moments)
    }
    dir.create(out_dir)
   # write.table('',file=file.path(out_dir, "log.log"))
    #LOG file path
    logfile.path = file.path(out_dir, "log.log")
    #Log conn
    log.con = file(logfile.path, open = "a")
    cat(paste("LOG FILE CREATED: ",Sys.time(), sep=""), file=log.con, sep="\n")
    tryCatch(radar_plot_values <- eval_var_loc(gene_sigs_list,names_sigs, mRNA_expr_matrix,names_datasets,out_dir,file=log.con,showResults,radar_plot_values),
             error=function(err){
               cat(paste0("Error occurred: ",err), file=log.con, sep="\n")
             })
    tryCatch(radar_plot_values <- eval_expr_loc(gene_sigs_list,names_sigs,mRNA_expr_matrix,names_datasets,thresholds, out_dir,file=log.con,showResults,radar_plot_values ),
             error=function(err){
               cat(paste0("Error occurred: ",err), file=log.con, sep="\n")
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
           cat(paste0("Error occurred: ",err), file=log.con, sep="\n")
         })
    tryCatch(radar_plot_values <- eval_struct_loc(gene_sigs_list,names_sigs,mRNA_expr_matrix,names_datasets,covariates,out_dir,file=log.con,showResults,radar_plot_values ),
             error=function(err){
               cat(paste0("Error occurred: ",err), file=log.con, sep="\n")
             })
    # print(radar_plot_values)
    tryCatch(make_radar_chart_loc(radar_plot_values,showResults,names_sigs, names_datasets,out_dir,file=log.con),
             error=function(err){
               cat(paste0("Error occurred: ",err), file=log.con, sep="\n")
             })

    if(!showResults)
      graphics.off()
    close(log.con)
  }else{
    print("Error: the length of names is not matching the number of elements in the expression matrices list.
          You need to have a name for every dataset used.")
  }
  # radar_plot_values
  # print(radar_plot_values)
}

draw.heatmaps <- function(hmaps,names_datasets,names_sigs){
  grid.newpage()
  num_rows <- length(names_sigs)#ceiling(sqrt(length(names)))
  num_cols <- length(names_datasets)#ceiling(length(names)/num_rows)

  pushViewport(viewport(layout = grid.layout(nrow = num_rows, ncol=num_cols)))
  count <- 1
  for (i in 1:num_rows){
    for(j in 1:num_cols){
      grid.draw(editGrob(hmaps[[count]], vp=viewport(layout.pos.row = i, layout.pos.col = j , clip=T)))
      count <- count +1
      # if (count  > length(names)){
      #   break;
      # }
    }
    # if (count  > length(names)){
    #   break;
    # }
  }
}

eval_stan_loc <- function(gene_sigs_list, names_sigs,mRNA_expr_matrix, names_datasets,out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values ){
  num_rows <- length(names_sigs)#ceiling(sqrt(length(names)))
  num_cols <- length(names_datasets)#ceiling(length(names)/num_rows)

  if (showResults){
    dev.new()
  }else{
    pdf(file.path(out_dir, 'sig_standardisation_comp.pdf'),width=10,height=10)
  }
  par(mfrow=c(num_rows,num_cols),oma=c(2,2,2,2),mar=c(4,4,4,4))
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets)){

      z_transf_mRNA <- mRNA_expr_matrix[[names_datasets[i]]][gene_sig,]
      for (j in 1:length(gene_sig)){
        z_transf_mRNA[gene_sig[j],] <- (as.numeric(z_transf_mRNA[gene_sig[j],]) - mean(as.numeric(z_transf_mRNA[gene_sig[j],]),na.rm=T)) / sd(as.numeric(z_transf_mRNA[gene_sig[j],]),na.rm=T)
      }
      z_transf_scores <- apply(z_transf_mRNA[gene_sig,],2,function(x) median(na.omit(x)))


      med_scores <- apply(mRNA_expr_matrix[[names_datasets[i]]][gene_sig,],2,function(x) median(na.omit(x)))

      combined_scores <- as.data.frame(cbind(med_scores,z_transf_scores))
      colnames(combined_scores) <- c('Median',"Z_Median")

      jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      smoothScatter(med_scores,z_transf_scores,colramp = jet.colors,xlab=NA,ylab=NA,main='Median vs Z-median')
      points(med_scores,z_transf_scores,pch='.')
      par(new=T,srt=45)
      plot(density(med_scores), axes=F, xlab=NA, ylab=NA,  col='red',main=NA)
      axis(side = 4)
      mtext(side = 4, line = 2, 'Density',cex=0.8)
      mtext(side = 2, line = 2, 'Z-median',cex=0.8)
      mtext(side = 1, line = 2, 'Median',cex=0.8)
      mtext(side=3,line=2.5,paste0(names_datasets[i],', ',names_sigs[k] ))

      rho <- cor(med_scores,z_transf_scores,method='spearman')
      rho_mean_med <- rho
      mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(med_scores))
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['standardization_comp'] <- rho
    }
  }
  cat('Standardisation compared successfully.\n', file=file)
  if(showResults){
    dev.copy(pdf,file.path(out_dir, 'sig_standardisation_comp.pdf'),width=10,height=10)
  }
  dev.off()
  radar_plot_values
}

eval_struct_loc <- function(gene_sigs_list,names_sigs, mRNA_expr_matrix,names_datasets,covariates, out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values){
  # library(gplots)
  require(biclust)
  require(ComplexHeatmap)
  # par(cex.main=0.7,cex.lab = 0.6,oma=c(2,2,3,2),mar=c(2,2,2,2))
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    hmaps <- sapply(1:length(names_datasets),function(i) {
      sig_scores <- (as.matrix(mRNA_expr_matrix[[names_datasets[i]]][gene_sig,]))
    #  print(sig_scores)
      tryCatch({
        if (length(covariates[[names_datasets[i]]]) ==0){

          #here we will set the parameters for the font size etc in the heatmaps

          dim.pdf = dim(sig_scores);
          # w = dim.pdf[2];
          h = dim.pdf[1];
          if(h<20){
            row_names.fontsize = 12
          }else{
            row_names.fontsize=5/log10(h)
          }

          ans_hmap <- Heatmap((sig_scores),show_column_dend = F,
                  show_column_names = F,
                  name=names_datasets[i],
                  heatmap_legend_param = list(title = names_datasets[i], color_bar = "continuous",legend_direction='vertical'),
                  column_title = paste0(names_datasets[i]),
                  row_names_gp =  gpar(fontsize = row_names.fontsize),
                  row_title = 'Genes')#,
        }else{
          if (is.vector(covariates[[names_datasets[i]]][['annotations']])){
            ha1 = HeatmapAnnotation(df = as.data.frame(covariates[[names_datasets[i]]][['annotations']][intersect(names(covariates[[names_datasets[i]]][['annotations']]),colnames(sig_scores))]),
                                    col=covariates[[names_datasets[i]]][['colors']],
                                    na_col="grey")#,
            #show_annotation_name = TRUE)#, col = list(type = c("a" = "red", "b" = "blue")
          }else{
            ha1 = HeatmapAnnotation(df = as.data.frame(covariates[[names_datasets[i]]][['annotations']][intersect(rownames(covariates[[names_datasets[i]]][['annotations']]),colnames(sig_scores)),]),
                                    col=covariates[[names_datasets[i]]][['colors']],
                                    na_col="grey",
                                    which="column")#,
            #show_annotation_name = TRUE)#, col = list(type = c("a" = "red", "b" = "blue"),
          }
           dim.pdf = dim(sig_scores);
          # w = dim.pdf[2];
          h = dim.pdf[1];
          if(h<20){
            row_names.fontsize = 12
          }else{
            row_names.fontsize=5/log10(h)
          }

          # row_names_gp = gpar(fontsize = row_names.fontsize)

          ans_hmap <- Heatmap((sig_scores),show_column_dend = F,
                  show_column_names = F,
                  name=names_datasets[i],
                  heatmap_legend_param = list(title = names_datasets[i], color_bar = "continuous",legend_direction='vertical'),
                  column_title = paste0(names_datasets[i]),
                  row_names_gp = gpar(fontsize = row_names.fontsize),
                  row_title = 'Genes', top_annotation=ha1)#,
    
        }
      },
      error=function(err){
        #print(paste0("There was an error, likely due to NA values in ",names[i] ," : ", err))
        cat(paste0('Error when creating expression heatmaps for ',names_datasets[i],' ', names_sigs[k],': ',err,'\n'), file=file)

        plot.new()
        title(paste0('\n \n \n',names_datasets[i],' ',names_sigs[k]))
      })
      ans_hmap

      #grab_grob()
    })
    for ( i in 1:length(names_datasets)){

      str_to_eval <- 'draw(hmaps[[1]]'
      if (length(names_datasets) > 1) {
        for (j in 2:length(names_datasets)){
          str_to_eval <- paste0(str_to_eval,' + hmaps[[',j,']]')
        }
      }

      str_to_eval <- paste0(str_to_eval,',heatmap_legend_side = \"left\",annotation_legend_side = \"left\", main_heatmap = \"',names_datasets[i],'\")')
      if (showResults){
        dev.new()
      } else{
       pdf(file.path(out_dir, paste0('sig_eval_struct_clustering_',names_datasets[i],'_',names_sigs[k],'.pdf')),width=10,height=10)#paste0(out_dir,'/sig_autocor_hmps.pdf'))  
      }

      eval(parse(text=str_to_eval))

      if(showResults){
        dev.copy(pdf,file.path(out_dir, paste0('sig_eval_struct_clustering_',names_datasets[i],'_',names_sigs[k],'.pdf')),width=10,height=10)#paste0(out_dir,'/sig_autocor_hmps.pdf'))
      }
      cat('Expression heatmaps saved successfully.\n', file=file)
      dev.off()
    }
  }

  dev.off()

  # draw.heatmaps(hmaps,names)

  if (showResults){
    dev.new()
  }
  
  par(cex.main=0.8,cex.lab = 0.8,oma=c(4,2,2,2),mar=c(4,4,4,4))
 
  hmaps <- lapply(1:(length(names_datasets)* length(names_sigs)),function(i) {

    dataset_ind <- i %% length(names_datasets)
    if (dataset_ind == 0 ){
      dataset_ind <- length(names_datasets)
    }

    sig_ind <- ceiling(i/length(names_datasets))
    gene_sig <- gene_sigs_list[[names_sigs[sig_ind]]]

    sig_scores <- as.matrix(mRNA_expr_matrix[[names_datasets[dataset_ind]]][gene_sig,])
    sig_scores[!is.finite(sig_scores)] <- NA

    #need to standardize here the matrix
    for (j in 1:length(gene_sig)){
      sig_scores[gene_sig[j],] <- (as.numeric(sig_scores[gene_sig[j],]) - mean(as.numeric(sig_scores[gene_sig[j],]),na.rm=T)) / sd(as.numeric(sig_scores[gene_sig[j],]),na.rm=T)
    }


    x <- binarize(na.omit(t(sig_scores)))#discretize(na.omit(t(sig_scores)),nof=10,quant=F)
    if (dim(sig_scores)[2] > 40) {
      num_cols_chosen <- 20
    }else if (dim(sig_scores)[2] > 20){
      num_cols_chosen <- 10

    }else if (dim(sig_scores)[2] > 10){
      num_cols_chosen <- 5
    }else{
      num_cols_chosen <- 2
    }

    # Xmotif <- biclust(x, method=BCXmotifs(), number=50, alpha=0.5, nd=floor(dim(sig_scores)/num_cols_chosen), ns=num_cols_chosen, sd=floor(dim(sig_scores)/num_cols_chosen))
    Xmotif <- biclust(x, method=BCCC(), delta=1,alpha=1.5, number=50)# alpha=0.5, nd=floor(dim(sig_scores)/num_cols_chosen), ns=num_cols_chosen, sd=floor(dim(sig_scores)/num_cols_chosen))

    if(Xmotif@Number > 1){
      heatmapBC(na.omit(t(sig_scores)),bicResult=Xmotif,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')
      #	heatmapBC(x,bicResult=Xmotif,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')

      title(paste0('\n \n \nBivariate clustering\n',names_datasets[dataset_ind],' ',names_sigs[sig_ind]))
      axis(1,at=1:length(rownames(sig_scores)), labels=rownames(sig_scores),las=2,tck=0,cex.axis=0.6)
      #mtext(rownames(sig_scores))
      # }else if (Xmotif@Number == 1){
      # 	heatmapBC(na.omit(t(sig_scores)),bicResult=Xmotif,number=1,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')
    }else{
      #print(paste0("Zero or one co-clusters found: ", names[i]))
      plot.new()
      title(paste0('\n\n\n <=1 bivariate clusters for\n',names_datasets[dataset_ind],' ',names_sigs[sig_ind]))
      cat(paste0('<= 1 bi-clusters found for: ', names_datasets[dataset_ind],' ',names_sigs[sig_ind],'\n'), file=file)

    }
    grab_grob()
  })

  draw.heatmaps(hmaps,names_datasets,names_sigs)

  dev.copy(pdf,file.path(out_dir,'sig_eval_bivariate_clustering.pdf'),width=10,height=10)
  cat('Bi-clustering completed successfully\n', file=file)

  dev.off()
  #same thing, but with binarized heatmaps underneath
  if (showResults){
    dev.new()
  }
  # pdf(file.path(out_dir,'sig_eval_bivariate_clustering_binarized_maps.pdf'),width=10,height=10)

  par(cex.main=0.8,cex.lab = 0.8,oma=c(4,2,2,2),mar=c(4,4,4,4))
  hmaps <- lapply(1:(length(names_datasets) * length(names_sigs)),function(i) {
    dataset_ind <- i %% length(names_datasets)
    if (dataset_ind == 0 ){
      dataset_ind <- length(names_datasets)
    }
    sig_ind <- ceiling(i/length(names_datasets))
    
    gene_sig <- gene_sigs_list[[names_sigs[sig_ind]]]

    sig_scores <- as.matrix(mRNA_expr_matrix[[names_datasets[dataset_ind]]][gene_sig,])
    sig_scores[!is.finite(sig_scores)] <- NA
    #standardise by z-transform
    for (j in 1:length(gene_sig)){
      sig_scores[gene_sig[j],] <- (as.numeric(sig_scores[gene_sig[j],]) - mean(as.numeric(sig_scores[gene_sig[j],]),na.rm=T)) / sd(as.numeric(sig_scores[gene_sig[j],]),na.rm=T)
    }
    x <- binarize(na.omit(t(sig_scores)))#discretize(na.omit(t(sig_scores)),nof=10,quant=F)
    if (dim(sig_scores)[2] > 40) {
      num_cols_chosen <- 20
    }else if (dim(sig_scores)[2] > 20){
      num_cols_chosen <- 10

    }else if (dim(sig_scores)[2] > 10){
      num_cols_chosen <- 5
    }else{
      num_cols_chosen <- 2
    }
    # Xmotif <- biclust(x, method=BCXmotifs(), number=50, alpha=0.5, nd=floor(dim(sig_scores)/num_cols_chosen), ns=num_cols_chosen, sd=floor(dim(sig_scores)/num_cols_chosen))
    Xmotif <- biclust(x, method=BCCC(), delta=1,alpha=1.5, number=50)
    if(Xmotif@Number > 1){
      #heatmapBC(na.omit(t(sig_scores)),bicResult=Xmotif,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')
      heatmapBC(x,bicResult=Xmotif,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')

      title(paste0('\n \n \nBivariate clustering\n',names_datasets[dataset_ind],' ',names_sigs[sig_ind]))
      axis(1,at=1:length(rownames(sig_scores)), labels=rownames(sig_scores),las=2,tck=0,cex.axis=0.6)
      #mtext(rownames(sig_scores))
      # }else if (Xmotif@Number == 1){
      # 	heatmapBC(na.omit(t(sig_scores)),bicResult=Xmotif,number=1,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')
    }else{
     # print(paste0("Zero or one co-clusters found: ", names[i]))
      plot.new()
      title(paste0('\n\n\n <=1 bivariate clusters for\n',names_datasets[dataset_ind],' ',names_sigs[sig_ind]))
    }
    grab_grob()
  })

  draw.heatmaps(hmaps,names_datasets,names_sigs)

  dev.copy(pdf,file.path(out_dir,'sig_eval_bivariate_clustering_binarized_maps.pdf'),width=10,height=10)
  dev.off()
  radar_plot_values
}

eval_var_loc <- function(gene_sigs_list,names_sigs, mRNA_expr_matrix,names_datasets, out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values){
  #calculate the number of rows and columns in the image
  num_rows <- length(names_sigs)#ceiling(sqrt(length(names)))
  num_cols <- length(names_datasets)#ceiling(length(names)/num_rows)
  # pdf(paste0(out_dir,'/sig_expr_var.pdf'))
  if (showResults){
    dev.new()
  }else{
    pdf(file.path(out_dir,'sig_mean_vs_sd.pdf'),width=10,height=10) 
  }
 
  gene_sig_mean_sd_table <- list()
  par(mfrow=c(num_rows,num_cols))
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets)){
      sd_genes <- apply(mRNA_expr_matrix[[names_datasets[i]]],1,function(x) sd(as.numeric(x),na.rm=T) )
      mean_genes <- apply(mRNA_expr_matrix[[names_datasets[i]]],1,function(x) mean(as.numeric(x),na.rm=T))
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['sd_median_ratio'] = median(na.omit(sd_genes[gene_sig]))/(median(na.omit(sd_genes)) +median(na.omit(sd_genes[gene_sig])))
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['abs_skewness_ratio'] = abs(skewness(na.omit(mean_genes[gene_sig])))/(abs(skewness(na.omit(mean_genes))) +  abs(skewness(na.omit(mean_genes[gene_sig]))))

      plot(mean_genes,sd_genes,pch=19,col='grey',
           main=paste0('Mean vs SD for all genes and signature genes\n',names_datasets[i], ' ',names_sigs[k]),
           xlab='Mean',
           ylab='Standard deviation')
      points(mean_genes[gene_sig],sd_genes[gene_sig],pch=19,col='red')
      quants_mean <- quantile(mean_genes*is.finite(mean_genes),probs=c(0.1,0.25,0.5,0.75,0.9),na.rm=T)
      abline(v=quants_mean[1],lty=3) #add line at 10%
      abline(v=quants_mean[2],lty=3) #add line at 25%
      abline(v=quants_mean[3],lty=3) #add line at 50%
      abline(v=quants_mean[4],lty=3) #add line at 75%
      abline(v=quants_mean[5],lty=3) #add line at 90%

      mtext(side = 3, line = 0, at=quants_mean[1], '10%',cex=0.4)
      mtext(side = 3, line = 0.4, at= quants_mean[2], '25%',cex=0.4)
      mtext(side = 3, line = 0, at= quants_mean[3], '50%',cex=0.4)
      mtext(side = 3, line = 0.4, at=quants_mean[4], '75%',cex=0.4)
      mtext(side = 3, line = 0, at = quants_mean[5], '90%',cex=0.4)

      quants_sd <- quantile(sd_genes*is.finite(sd_genes),probs=c(0.1,0.25,0.5,0.75,0.9),na.rm=T)
      abline(h=quants_sd[1],lty=3) #add line at 10%
      abline(h=quants_sd[2],lty=3) #add line at 25%
      abline(h=quants_sd[3],lty=3) #add line at 50%
      abline(h=quants_sd[4],lty=3) #add line at 75%
      abline(h=quants_sd[5],lty=3) #add line at 90%

      mtext(side = 4, line = 0, at=quants_sd[1], '10%',cex=0.4)
      mtext(side = 4, line = 0.4, at= quants_sd[2], '25%',cex=0.4)
      mtext(side = 4, line = 0.8, at= quants_sd[3], '50%',cex=0.4)
      mtext(side = 4, line = 0.4, at=quants_sd[4], '75%',cex=0.4)
      mtext(side = 4, line = 0, at = quants_sd[5], '90%',cex=0.4)

      gene_sig_mean_sd_table[[names_sigs[k]]][[names_datasets[i]]] <- cbind(mean_genes[gene_sig],sd_genes[gene_sig])
      colnames(gene_sig_mean_sd_table[[names_sigs[k]]][[names_datasets[i]]]) <- c("Mean","SD")

    }
  }
  if(showResults){
   dev.copy(pdf,file.path(out_dir,'sig_mean_vs_sd.pdf'),width=10,height=10)  
  }
  dev.off()
  cat('Mean vs SD graphs created successfully.\n', file=file)


  #now let's output the tables, one for each dataset considered
  dir.create(file.path(out_dir,'mean_sd_tables'))
  for(k in 1:length(names_sigs)){
    for (i in 1:length(names_datasets)){
      write.table(gene_sig_mean_sd_table[[names_sigs[k]]][[names_datasets[i]]],file=file.path(out_dir,'mean_sd_tables', paste0('mean_sd_table_',names_sigs[k],'_',names_datasets[i],'.txt')),quote=F,sep='\t')

    }
  }
   cat('Mean vs SD tables written to file successfully.\n', file=file)


  #the following is the code for the original variable boxplots
  # 	dev.new()
  # par(mfrow=c(num_rows,num_cols))
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets)){
    	coeff_of_var <- apply(mRNA_expr_matrix[[names_datasets[i]]],1,function(x) sd(as.numeric(na.omit(x)),na.rm=T) / mean(as.numeric(na.omit(x)),na.rm=T))
    	coeff_of_var_gene_sig <- coeff_of_var[gene_sig]
      quantiles_considered <- quantile(coeff_of_var,probs=c(0.9,0.75,0.5),na.rm=T)
      prop_top_10_percent <- sum(na.omit(coeff_of_var_gene_sig) >= quantiles_considered[1]) / length(na.omit(coeff_of_var_gene_sig))
      prop_top_25_percent <- sum(na.omit(coeff_of_var_gene_sig) >= quantiles_considered[2]) / length(na.omit(coeff_of_var_gene_sig))
      prop_top_50_percent <- sum(na.omit(coeff_of_var_gene_sig) >= quantiles_considered[3]) / length(na.omit(coeff_of_var_gene_sig))
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['prop_top_10_percent'] <- prop_top_10_percent
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['prop_top_25_percent'] <- prop_top_25_percent
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['prop_top_50_percent'] <- prop_top_50_percent
      
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['coeff_of_var_ratio'] <- median(na.omit(coeff_of_var_gene_sig))/(median(na.omit(coeff_of_var)) + median(na.omit(coeff_of_var_gene_sig)))

    	# boxplot(na.omit(coeff_of_var),na.omit(coeff_of_var_gene_sig),#log="y",
    		# names=c('All Genes','Gene Signature'),
    		# ylab='Coefficient of Variation',
    		# main=paste0('Variance of signature genes vs. all genes\n',names[i]))
    }
  }
  # dev.copy(pdf,paste0(out_dir,'/sig_expr_var.pdf'))

  # dev.off()
  radar_plot_values
}


eval_expr_loc <- function(gene_sigs_list,names_sigs, mRNA_expr_matrix,names_datasets, thresholds = NULL, out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values){
  #calculate the number of rows and columns in the image
  num_rows <- length(names_sigs)#ceiling(sqrt(length(names)))
  num_cols <- length(names_datasets)#ceiling(length(names)/num_rows)

  if (showResults){
    dev.new()
  }else{
    pdf(file.path(out_dir,'sig_expr_barcharts_NA_values.pdf'),width=10,height=10)

  }

  par(mfrow=c(num_rows,num_cols),cex=0.7, cex.axis=0.5)
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets)){
      #calculate the porportion of non-NA expression data in the matrix
      genes_expr <- mRNA_expr_matrix[[names_datasets[i]]][gene_sig,]
      gene_expr_vals <- (rowSums(is.na(genes_expr)) / dim(genes_expr)[2])
      gene_expr_vals <- -sort(-gene_expr_vals)
      if(max(gene_expr_vals)==0){
        bar_expr <- barplot(gene_expr_vals,xlab="Signature Gene IDs",ylab="Proportion of NA expression ",main=paste0("Signature gene expression across samples\n",names_datasets[i],' ',names_sigs[k]), axisnames=F,axis=F,ylim=c(0,1))
      }else{
        bar_expr <- barplot(gene_expr_vals,xlab="Signature Gene IDs",ylab="Proportion of NA expression ",main=paste0("Signature gene expression across samples\n",names_datasets[i],' ',names_sigs[k]), axisnames=F,axis=F)

      }
      text(bar_expr, par("usr")[3], labels = names(gene_expr_vals), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.5)
      axis(2)
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['med_prop_na'] <- median(1-gene_expr_vals)

    }
  }
  if(showResults){
    dev.copy(pdf,file.path(out_dir,'sig_expr_barcharts_NA_values.pdf'),width=10,height=10)
  }
  dev.off()

  if (showResults){
    dev.new()
  }else{
    pdf(file.path(out_dir,'sig_expr_barcharts.pdf'),width=10,height=10) 
  }
 
  par(mfrow=c(num_rows,num_cols),cex=0.7, cex.axis=0.5)
  if (length(thresholds) == 0) {
    thresholds <- rep(0,length(names_datasets))
    for ( i in 1:length(names_datasets)){
      thresholds[i] <- median(unlist(na.omit(mRNA_expr_matrix[[names_datasets[i]]])))
    }
    #names(thresholds) <- names_datasets
  }
  if(length(names(thresholds))==0){
    names(thresholds) <- names_datasets
  }
  for (k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets)){
      #calculate the porportion of nonzero expression data in the matrix
      genes_expr <- mRNA_expr_matrix[[names_datasets[i]]][gene_sig,]
      gene_expr_vals <- 1 - ((rowSums(genes_expr < thresholds[i])) / (dim(genes_expr)[2]))
      gene_expr_vals <- sort(gene_expr_vals)
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['med_prop_above_med'] <- median(gene_expr_vals)
      bar_expr <- barplot(gene_expr_vals,xlab="Signature Gene IDs",ylab="Proportion with expression above threshold",main=paste0("Signature gene expression across samples\n",names_datasets[i],' ',names_sigs[k]), axisnames=F,axis=F)
      text(bar_expr, par("usr")[3], labels = names(gene_expr_vals), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.5)
      axis(2)
    }
  }
  if(showResults){
    dev.copy(pdf,file.path(out_dir,'sig_expr_barcharts.pdf'),width=10,height=10)
  }
  dev.off()


  if (showResults){
    dev.new()
  }else{
    pdf(file.path(out_dir,'sig_expr_density_plots.pdf'),width=10,height=10)
  }

  par(mfrow=c(num_rows,num_cols))
  for (k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets)){
      #calculate the porportion of nonzero expression data in the matrix
      genes_expr <- mRNA_expr_matrix[[names_datasets[i]]][gene_sig,]
      gene_expr_vals <- 1 - (rowSums(genes_expr < thresholds[i]) / (dim(genes_expr)[2]))
      plot(density(na.omit(gene_expr_vals),adjust=0.25),main=paste0("Signature gene expression across samples\n",names_datasets[i], ' ',names_sigs[k]),ylab="Density")
      
    }
  }
  if(showResults){
    dev.copy(pdf,file.path(out_dir,'sig_expr_density_plots.pdf'),width=10,height=10)
  }
  dev.off()
  cat('Expression and density graphs created successfully.\n', file=file)

 # print(paste0("Min expression occurs for gene ID ",names(gene_expr_vals)[1], ", with ", round(gene_expr_vals[1]*100,2),"% non-zero expression."))
  radar_plot_values
}


compare_metrics_loc <- function(gene_sigs_list,names_sigs, mRNA_expr_matrix, names_datasets, out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values){
  require(gplots)
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    hmaps <- list()
    if (showResults){
      dev.new()
    }else{
      pdf(file.path(out_dir,paste0('sig_compare_metrics_',names_sigs[k],'.pdf')),width=10,height=10)
    }

    par(mfcol = c(4,length(names_datasets)),mar=c(4,4,4,4))
    for ( i in 1:length(names_datasets)){
      mRNA_expr_matrix[[names_datasets[i]]][!(is.finite(as.matrix(mRNA_expr_matrix[[names_datasets[i]]])))] <- NA

      med_scores <- apply(mRNA_expr_matrix[[names_datasets[i]]][gene_sig,],2,function(x) median(na.omit(x)))
      mean_scores <- apply(mRNA_expr_matrix[[names_datasets[i]]][gene_sig,],2,function(x) mean(na.omit(x)))
      pca1_scores <- prcomp(na.omit(t(mRNA_expr_matrix[[names_datasets[i]]][gene_sig,])),retx=T)
      props_of_variances <- pca1_scores$sdev^2/(sum(pca1_scores$sdev^2))
      pca1_scores <- pca1_scores$x[,1]

      common_score_cols <- intersect(names(med_scores),intersect(names(mean_scores),names(pca1_scores)))
      jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      smoothScatter(med_scores[common_score_cols],mean_scores[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='Median vs Mean')
      points(med_scores[common_score_cols],mean_scores[common_score_cols],pch='.')
      par(new=T)#,srt=45)
      plot(density(med_scores[common_score_cols]), axes=F, xlab=NA, ylab=NA,  col='red',main=NA)
      axis(side = 4)
      mtext(side = 4, line = 2, 'Density',cex=0.8)
      mtext(side = 2, line = 2, 'Mean',cex=0.8)
      mtext(side = 1, line = 2, 'Median',cex=0.8)
      mtext(side=3,line=2.5,paste0(names_datasets[i],' ',names_sigs[k]))

      rho <- cor(med_scores[common_score_cols],mean_scores[common_score_cols],method='spearman')
      rho_mean_med <- rho
      mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(med_scores[common_score_cols]))
      smoothScatter(mean_scores[common_score_cols],pca1_scores[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='Mean vs PCA1')
      points(mean_scores[common_score_cols],pca1_scores[common_score_cols],pch='.')
      par(new=T)
      plot(density(mean_scores[common_score_cols]), axes=F, xlab=NA, ylab=NA, col='red',main=NA)
      axis(side = 4)
      mtext(side = 4, line = 2, 'Density',cex=0.8)
      rho <- cor(mean_scores[common_score_cols],pca1_scores[common_score_cols],method='spearman')
      rho_mean_pca1 <- rho
      mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(mean_scores[common_score_cols]))
      mtext(side = 2, line = 2, 'PCA1',cex=0.8)
      mtext(side = 1, line = 2, 'Mean',cex=0.8)

      smoothScatter(pca1_scores[common_score_cols],med_scores[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='PCA1 vs Median')
      points(pca1_scores[common_score_cols],med_scores[common_score_cols],pch='.')
      par(new=T)
      plot(density(pca1_scores[common_score_cols]), axes=F, xlab=NA, ylab=NA, col='red',main=NA)
      axis(side = 4)
      mtext(side = 4, line = 2, 'Density',cex=0.8)
      rho <- cor(pca1_scores[common_score_cols],med_scores[common_score_cols],method='spearman')
      rho_pca1_med <- rho
      mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(pca1_scores[common_score_cols]))
      mtext(side = 2, line = 2, 'Median',cex=0.8)
      mtext(side = 1, line = 2, 'PCA1',cex=0.8)

      #here we put the cor plot of the mean median and pca1

      autocors_mat <- matrix(0,nrow=3,ncol=3)
      row.names(autocors_mat) <- c('Mean' , 'Median','PCA1')
      colnames(autocors_mat) <- c('Mean' , 'Median','PCA1')
      autocors_mat[1,1] <- 1
      autocors_mat[1,2] <- rho_mean_med
      autocors_mat[1,3] <- rho_mean_pca1

      autocors_mat[2,1] <- rho_mean_med
      autocors_mat[2,2] <- 1
      autocors_mat[2,3] <- rho_pca1_med

      autocors_mat[3,1] <- rho_mean_pca1
      autocors_mat[3,2] <- rho_pca1_med
      autocors_mat[3,3] <- 1

      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['rho_mean_med'] <- rho_mean_med
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['rho_pca1_med'] <- rho_pca1_med
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['rho_mean_pca1'] <- rho_mean_pca1
      

      # hmaps[[i]] <- heatmap.2(autocors_mat,
      #         	col = colorpanel(100,"blue","white","red"),#colorpanel(100,"red","yellow","green"),
      #          	trace = "none",
      #          	xlab = NA,
      #          	na.color="grey",
      #          	labRow=rownames(autocors_mat),
      #       	labCol=colnames(autocors_mat),
      #          	main = paste("Scoring metric correlation"),
      #          	dendrogram = "none",
      #          	breaks=seq(-1,1,length.out=101) ,
      #         	symbreaks = T,
      #        	Rowv = NA,Colv=NA)# ,margins = c(4,4))

      bars_plot <- props_of_variances[1:min(10,length(props_of_variances))]
      barplot(bars_plot,main="PCA vs proportion\n of variance") #ylim= c(0,1),
      mtext(side = 1, line = 2, 'PCA',cex=0.8)
    }
    if(showResults){
      dev.copy(pdf,file.path(out_dir,paste0('sig_compare_metrics_',names_sigs[k],'.pdf')),width=10,height=10)
    }
    dev.off()
  }
  cat('Metrics compared successfully.\n', file=file)
  radar_plot_values
}


grab_grob <- function(){
  require(gridGraphics)
  require(grid)
  grid.echo()
  grid.grab()
}

eval_compactness_loc <- function(gene_sigs_list,names_sigs, mRNA_expr_matrix, names_datasets, out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values){
  require(gplots)
  if (showResults){
    dev.new()
  }
  # pdf(file.path(out_dir,'sig_autocor_hmps.pdf'),width=10,height=10)

  par(cex.main=0.8,cex.lab = 0.6,oma=c(2,0,0,0),mar=c(0,0,0,0))
  hmaps <- lapply(1:(length(names_sigs) *length(names_datasets)),function(i) {
    dataset_ind <- i %% length(names_datasets)
    if (dataset_ind == 0 ){
      dataset_ind <- length(names_datasets)
    }
    sig_ind <- ceiling(i/length(names_datasets))
    gene_sig <- gene_sigs_list[[names_sigs[sig_ind]]]
    autocors <- cor(t(na.omit(mRNA_expr_matrix[[names_datasets[dataset_ind]]][intersect(gene_sig,rownames(mRNA_expr_matrix[[names_datasets[dataset_ind]]])),])),method='spearman')

    tryCatch({
      heatmap.2( na.omit(autocors),
                 col = colorpanel(100,"blue","white","red"), #redgreen(100),#colorpanel(100,"red","yellow","green"),
                 trace = "none",
                 xlab = "Gene ID",
                 ylab="Gene ID",
                 na.color="grey",
                 labRow=rownames(autocors),
                 labCol=colnames(autocors),#gene_sig,
                 main = paste0("\n\nAutocorrelation ", names_datasets[[dataset_ind]] ,' ',names_sigs[[sig_ind]]),
                 dendrogram = "col",
                 symbreaks = T,
                 Rowv = T,Colv=T ,key.xlab='Rho',key.ylab=NA,  key.title=NA,cexRow=0.5,cexCol=0.5,margins=c(4,4))
    },
    error=function(err){
#      print(paste0("There was an error, likely due to NA values in data: ", err))
      cat(paste0("There was an error, likely due to NA values in data, for dataset: ",names_datasets[dataset_ind]," ", names_sigs[sig_ind]," ", err,'\n'), file=file)

    })
    grab_grob()
  })
  draw.heatmaps(hmaps,names_datasets,names_sigs)
  dev.copy(pdf,file.path(out_dir,'sig_autocor_hmps.pdf'),width=10,height=10)
  dev.off()

  if (showResults){
    dev.new()
  }else{
    pdf(file.path(out_dir,'sig_autocor_dens.pdf'),width=10,height=10)
  }

  par(cex.main=0.8,cex.lab = 0.6,oma=c(2,2,2,2),mar=c(2,2,2,2),mfrow=c(1,1))
  max_dens <- -9999
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets) ){
      autocors <- cor(na.omit(t(mRNA_expr_matrix[[names_datasets[i]]][intersect(gene_sig,rownames(mRNA_expr_matrix[[names_datasets[i]]])),])),method='spearman')
      cur_max <- max(density(unlist(na.omit(autocors)))$y)
      if (max_dens  < cur_max){
        max_dens <- cur_max
      }
    }
  }

  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets) ){
      autocors <- cor(t(na.omit(mRNA_expr_matrix[[names_datasets[i]]][gene_sig,])),method='spearman')
      if ((i ==1) && (k==1)){
        plot(density(unlist(na.omit(autocors))),ylim=c(0,ceiling(max_dens)),col=i,main=NA,lwd=2,lty=k)
      }else{
        lines(density(unlist(na.omit(autocors))),col=i,main=NA,lwd=2,lty=k)
      }
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['autocor_median'] <- median(na.omit(autocors))

    }
  }

  mtext(side = 2, line = 2, 'Density',cex=0.8)
  mtext(side = 1, line = 2, 'Rho',cex=0.8)
  mtext(side=3,line=2,'Autocorrelation Density')

  op <- par(cex=0.6)
  legend_names <- c()
  legend_cols <- c()
  legend_lty <- c()
  for(k in 1:length(names_sigs)){
    for (i in 1:length(names_datasets) ){
      legend_names <- c(legend_names,paste0(names_datasets[i],' ',names_sigs[k]))
      legend_cols <- c(legend_cols,i)
      legend_lty <- c(legend_lty,k)
    }
  }
  legend("topright",legend_names,col=legend_cols,lty=legend_lty,lwd=rep(2,times=(length(names_datasets) * length(names_sigs))))

  if(showResults){
    dev.copy(pdf,file.path(out_dir,'sig_autocor_dens.pdf'),width=10,height=10)
  }
  dev.off()
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (showResults){
      dev.new()
    }else{
      pdf(file.path(out_dir,paste0('sig_autocor_rankProd_',names_sigs[k],'.pdf')),width=10,height=10)
    }

    par(cex.main=0.8,cex.lab = 0.6,oma=c(2,2,2,2),mar=c(4,4,4,4))
    cat("Autocorrelation metrics successfully computed.\n", file=file)

    #now we take the median of the genes' autocorrelation for each gene in each dataset and then look at the rank product over the different cancer types
    if (length(names_datasets) > 1){
      overall_rank_mat <- matrix(NA,nrow=length(unique(gene_sig)),ncol=length(names_datasets))
      row.names(overall_rank_mat) <- unique(gene_sig)
      colnames(overall_rank_mat) <- names_datasets
      for (i in 1:length(names_datasets)){
        autocors <- cor(t(na.omit(mRNA_expr_matrix[[names_datasets[i]]][intersect(unique(gene_sig),rownames(mRNA_expr_matrix[[names_datasets[i]]])),])),method='spearman')
        median_scores <- apply(autocors,2,function(x) median(na.omit(x)))
        overall_rank_mat[names(median_scores),i] <- median_scores
      }

      require(RankProd)
      RP.out <-RP(data = overall_rank_mat,cl = rep(1,times=length(names_datasets)),logged = F,gene.names=rownames(overall_rank_mat))
      plotRP(RP.out ,cutoff=0.05)
      topGene(RP.out,cutoff=0.05,method="pfp",logged=F, gene.names=intersect(gene_sig,rownames(mRNA_expr_matrix[[names_datasets[i]]])))
      cat("Autocorrelation rank product successfully computed.\n", file=file)

    }else{
        cat("Rank product not computed as there is only one dataset.\n", file=file)

    }
    if(showResults){
      dev.copy(pdf,file.path(out_dir,paste0('sig_autocor_rankProd_',names_sigs[k],'.pdf')),width=10,height=10)
    }
    dev.off()
  }
  radar_plot_values

}

make_radar_chart_loc <- function(radar_plot_values,showResults,names_sigs,names_datasets, out_dir = '~',file){
  radar_plot_mat <- c()
  # print(radar_plot_values)
  for(k in 1:length(names_sigs)){
    t <- lapply(radar_plot_values[[names_sigs[k]]],function(x) radar_plot_mat <<- rbind(radar_plot_mat,x))
  }
  radar_plot_mat <- rbind(rep(1,length(radar_plot_values[[1]][[1]])),rep(0,length(radar_plot_values[[1]][[1]])),radar_plot_mat)
  radar_plot_mat <- abs(radar_plot_mat)
  # print(radar_plot_mat)
  legend_labels <- c()
  legend_cols <- c()
  legend_lty <- c()
  for(k in 1:length(names_sigs) ){
    for(i in 1:length(names_datasets)){
      legend_labels <- c(legend_labels,paste0(names_datasets[i],' ',names_sigs[k]))
      legend_cols <- c(legend_cols,i)
      legend_lty <- c(legend_lty, k)
    }
  }
  # print(legend_labels)
  row.names(radar_plot_mat) <- c('max','min',legend_labels)#names(radar_plot_values))
  # print(radar_plot_mat)
   if (showResults){
      dev.new()
    }else{
      pdf(file.path(out_dir,'sig_radarplot.pdf'),width=10,height=10)
    }

  areas <- c()
  for (i in 3:dim(radar_plot_mat)[1]){
    areas<- c(areas,sum(sapply(1:length(radar_plot_mat[i,]),function(x) if(x < length(radar_plot_mat[i,])){radar_plot_mat[i,x] * radar_plot_mat[i,x+1]}else{radar_plot_mat[i,x]* radar_plot_mat[i,1]})))

  }
  areas <- areas /dim(radar_plot_mat)[2]
  # print(areas)
 
#  legend_labels <- cbind(1:length(names(radar_plot_values)),paste0(names(radar_plot_values),' (',format(areas,digits=2),')'))
  count <-1
  legend_labels <-c()
   for(k in 1:length(names_sigs) ){
    for(i in 1:length(names_datasets)){
      legend_labels <- c(legend_labels,paste0(names_datasets[i],' ',names_sigs[k],' (',format(areas[count],digits=2),')'))
      count <- count + 1
    }
  }
  layout(rbind(1,2), heights=c(7,1))  # put legend on bottom 1/8th of the chart

  radarchart(as.data.frame(radar_plot_mat),
    maxmin = T,axistype = 1,
    cglcol = 'grey',axislabcol = 'black',
    caxislabels = seq(0,1,length.out = 5),
    cglty = 1,cglwd = 1,calcex = 0.5,
    vlabels = c('Ratio of\nMed. SD','Skew Ratio','Prop in\ntop 10% var','Prop in\ntop 25% var','Prop in\ntop 50% var','Coef. of Var. \nRatio',
      'Med. non-NA Prop.','Med. Prop.\nExpressed',
      'Med. Autocor.','Mean, Med.\nScore Cor.',
      'PCA1, Med.\nScore Cor.','Mean, PCA1\nScore Cor.',
      'Med., Z_Med.\nScore Cor.'),
    vlcex = 0.8,
    title='Signature Summary',
    pty=16, plty=legend_lty,pcol=legend_cols,plwd = 2)
    legend_labels <- legend_labels[order(-areas)]
    legend_cols <- legend_cols[order(-areas)]
    legend_lty <- legend_lty[order(-areas)]
    par(mar=c(0, 0, 0, 0))

    plot.new()
    par(xpd=TRUE)

    legend('center', legend=legend_labels, seg.len=0.5, title="Datasets",lty=legend_lty,pch=1, 
       bty="n" ,lwd=1, col=legend_cols,cex=0.6,horiz=T)
  if(showResults){
    dev.copy(pdf,file.path(out_dir,'sig_radarplot.pdf'),width=10,height=10)
  }
  dev.off()
  cat('Radar chart made successfully.\n', file=file)

}

