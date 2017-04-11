# eval_struct_loc.R
#
# This function creates the plots to evaluate for signature structure. Produces plots for the
# expression of the signature with hierarchical clustering as well as biclustering plots on binarized data and on non-binarized data
# @param gene_sigs_list A list of genes representing the gene signature to be tested.
# @param names_sigs The names of the gene signatures (one name per gene signature, in gene_sigs_list)
# @param mRNA_expr_matrix A list of expression matrices
# @param names_datasets The names of the different datasets contained in mRNA_expr_matrix
# @param covariates A list containing a sub-list of 'annotations' and 'colors' which contains the annotation matrix for the given dataset and the associated colours with which to plot in the expression heatmap
# @param out_dir A path to the directory where the resulting output files are written
# @param file File representing the log file where errors can be written
# @param showResults Tells if open dialog boxes showing the computed results. Default is FALSE
# @param radar_plot_values A list of values that store computations that will be used in the final summary radarplot
# @keywords eval_struct_loc

eval_struct_loc <- function(gene_sigs_list,names_sigs, mRNA_expr_matrix,names_datasets,covariates, out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values){
  # library(gplots)
  # require(biclust)
  # require(ComplexHeatmap)
  # par(cex.main=0.7,cex.lab = 0.6,oma=c(2,2,3,2),mar=c(2,2,2,2))
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    hmaps <- sapply(1:length(names_datasets),function(i) {

      data.matrix = mRNA_expr_matrix[[names_datasets[i]]]
      inter <- intersect(gene_sig[,1], row.names(data.matrix))

      sig_scores <- (as.matrix(data.matrix[inter,]))
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

          ans_hmap <- ComplexHeatmap::Heatmap((sig_scores),show_column_dend = F,
                                              show_column_names = F,
                                              name=names_datasets[i],
                                              heatmap_legend_param = list(title = names_datasets[i], color_bar = "continuous",legend_direction='vertical'),
                                              column_title = paste0(names_datasets[i]),
                                              row_names_gp =  grid::gpar(fontsize = row_names.fontsize),
                                              row_title = 'Genes')#,
        }else{
          if (is.vector(covariates[[names_datasets[i]]][['annotations']])){
            ha1 = ComplexHeatmap::HeatmapAnnotation(df = as.data.frame(covariates[[names_datasets[i]]][['annotations']][intersect(names(covariates[[names_datasets[i]]][['annotations']]),colnames(sig_scores))]),
                                                    col=covariates[[names_datasets[i]]][['colors']],
                                                    na_col="grey")#,
            #show_annotation_name = TRUE)#, col = list(type = c("a" = "red", "b" = "blue")
          }else{
            ha1 = ComplexHeatmap::HeatmapAnnotation(df = as.data.frame(covariates[[names_datasets[i]]][['annotations']][intersect(rownames(covariates[[names_datasets[i]]][['annotations']]),colnames(sig_scores)),]),
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

          # row_names_gp = grid::gpar(fontsize = row_names.fontsize)

          ans_hmap <- ComplexHeatmap::Heatmap((sig_scores),show_column_dend = F,
                                              show_column_names = F,
                                              name=names_datasets[i],
                                              heatmap_legend_param = list(title = names_datasets[i], color_bar = "continuous",legend_direction='vertical'),
                                              column_title = paste0(names_datasets[i]),
                                              row_names_gp = grid::gpar(fontsize = row_names.fontsize),
                                              row_title = 'Genes', top_annotation=ha1)#,

        }
      },
      error=function(err){
        #print(paste0("There was an error, likely due to NA values in ",names[i] ," : ", err))
        cat(paste0('Error when creating expression heatmaps for ',names_datasets[i],' ', names_sigs[k],': ',err,'\n'), file=file)

        graphics::plot.new()
        graphics::title(paste0('\n \n \n',names_datasets[i],' ',names_sigs[k]))
      })
      ans_hmap

      #grab_grob()
    })
    for ( i in 1:length(names_datasets)){
      all_hmaps <- hmaps[[1]]

      if (length(hmaps) > 1){
        tmp <- lapply(hmaps[2:length(hmaps)],function(x) all_hmaps <<- ComplexHeatmap::add_heatmap(all_hmaps,x))
      }
      # str_to_eval <- 'ComplexHeatmap::draw(hmaps[[1]]'
      # if (length(names_datasets) > 1) {
      #   for (j in 2:length(names_datasets)){
      #     str_to_eval <- paste0(str_to_eval,' + hmaps[[',j,']]')
      #   }
      # }

      # str_to_eval <- paste0(str_to_eval,',heatmap_legend_side = \"left\",annotation_legend_side = \"left\", main_heatmap = \"',names_datasets[i],'\")')
      if (showResults){
        grDevices::dev.new()
      } else{
        grDevices::pdf(file.path(out_dir, paste0('sig_eval_struct_clustering_',names_datasets[i],'_',names_sigs[k],'.pdf')),width=10,height=10)#paste0(out_dir,'/sig_autocor_hmps.pdf'))
      }
      ComplexHeatmap::draw(all_hmaps)

      # print(str_to_eval)
      # eval(parse(text=str_to_eval))
      # ComplexHeatmap::draw(hmaps[[1]] + hmaps[[2]])
      if(showResults){
        grDevices::dev.copy(grDevices::pdf,file.path(out_dir, paste0('sig_eval_struct_clustering_',names_datasets[i],'_',names_sigs[k],'.pdf')),width=10,height=10)#paste0(out_dir,'/sig_autocor_hmps.pdf'))
      }
      cat('Expression heatmaps saved successfully.\n', file=file)
      grDevices::dev.off()
    }

  }

  grDevices::dev.off()

  # draw.heatmaps(hmaps,names)

  if (showResults){
    grDevices::dev.new()
  }

  graphics::par(cex.main=0.8,cex.lab = 0.8,oma=c(4,2,2,2),mar=c(4,4,4,4))

  hmaps <- lapply(1:(length(names_datasets)* length(names_sigs)),function(i) {

    dataset_ind <- i %% length(names_datasets)
    if (dataset_ind == 0 ){
      dataset_ind <- length(names_datasets)
    }

    sig_ind <- ceiling(i/length(names_datasets))
    gene_sig <- gene_sigs_list[[names_sigs[sig_ind]]]

    data.matrix = mRNA_expr_matrix[[names_datasets[dataset_ind]]]
    inter = intersect(gene_sig[,1], row.names(data.matrix))
    sig_scores <- as.matrix(data.matrix[inter,])
    sig_scores[!is.finite(sig_scores)] <- NA

    #need to standardize here the matrix
    for (gene in inter){
      sig_scores[gene,] <- (as.numeric(sig_scores[gene,]) - mean(as.numeric(sig_scores[gene,]),na.rm=T)) / stats::sd(as.numeric(sig_scores[gene,]),na.rm=T)
    }


    x <- biclust::binarize(stats::na.omit(t(sig_scores)))#discretize(stats::na.omit(t(sig_scores)),nof=10,quant=F)
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
    Xmotif <- biclust::biclust(x, method=biclust::BCCC(), delta=1,alpha=1.5, number=50)# alpha=0.5, nd=floor(dim(sig_scores)/num_cols_chosen), ns=num_cols_chosen, sd=floor(dim(sig_scores)/num_cols_chosen))

    if(Xmotif@Number > 1){
      biclust::heatmapBC(stats::na.omit(t(sig_scores)),bicResult=Xmotif,col = gplots::colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')
      #	biclust::heatmapBC(x,bicResult=Xmotif,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')

      graphics::title(paste0('\n \n \nBivariate clustering\n',names_datasets[dataset_ind],' ',names_sigs[sig_ind]))
      graphics::axis(1,at=1:length(rownames(sig_scores)), labels=rownames(sig_scores),las=2,tck=0,cex.axis=0.6)
      #mtext(rownames(sig_scores))
      # }else if (Xmotif@Number == 1){
      # 	biclust::heatmapBC(stats::na.omit(t(sig_scores)),bicResult=Xmotif,number=1,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')
    }else{
      #print(paste0("Zero or one co-clusters found: ", names[i]))
      graphics::plot.new()
      graphics::title(paste0('\n\n\n <=1 bivariate clusters for\n',names_datasets[dataset_ind],' ',names_sigs[sig_ind]))
      cat(paste0('<= 1 bi-clusters found for: ', names_datasets[dataset_ind],' ',names_sigs[sig_ind],'\n'), file=file)

    }
    grab_grob()
  })

  draw.heatmaps(hmaps,names_datasets,names_sigs)

  grDevices::dev.copy(grDevices::pdf,file.path(out_dir,'sig_eval_bivariate_clustering.pdf'),width=10,height=10)
  cat('Bi-clustering completed successfully\n', file=file)

  grDevices::dev.off()
  #same thing, but with binarized heatmaps underneath
  if (showResults){
    grDevices::dev.new()
  }
  # pdf(file.path(out_dir,'sig_eval_bivariate_clustering_binarized_maps.pdf'),width=10,height=10)

  graphics::par(cex.main=0.8,cex.lab = 0.8,oma=c(4,2,2,2),mar=c(4,4,4,4))
  hmaps <- lapply(1:(length(names_datasets) * length(names_sigs)),function(i) {
    dataset_ind <- i %% length(names_datasets)
    if (dataset_ind == 0 ){
      dataset_ind <- length(names_datasets)
    }
    sig_ind <- ceiling(i/length(names_datasets))

    gene_sig <- gene_sigs_list[[names_sigs[sig_ind]]]

    data.matrix = mRNA_expr_matrix[[names_datasets[dataset_ind]]]
    inter = intersect(gene_sig[,1], row.names(data.matrix))

    sig_scores <- as.matrix(data.matrix[inter,])
    sig_scores[!is.finite(sig_scores)] <- NA
    #standardise by z-transform
    for (gene in gene_sig){
      sig_scores[gene,] <- (as.numeric(sig_scores[gene,]) - mean(as.numeric(sig_scores[gene,]),na.rm=T)) / stats::sd(as.numeric(sig_scores[gene,]),na.rm=T)
    }
    x <- biclust::binarize(stats::na.omit(t(sig_scores)))#discretize(stats::na.omit(t(sig_scores)),nof=10,quant=F)
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
    Xmotif <- biclust::biclust(x, method=biclust::BCCC(), delta=1,alpha=1.5, number=50)
    if(Xmotif@Number > 1){
      #biclust::heatmapBC(stats::na.omit(t(sig_scores)),bicResult=Xmotif,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')
      biclust::heatmapBC(x,bicResult=Xmotif,col = gplots::colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')

      graphics::title(paste0('\n \n \nBivariate clustering\n',names_datasets[dataset_ind],' ',names_sigs[sig_ind]))
      graphics::axis(1,at=1:length(rownames(sig_scores)), labels=rownames(sig_scores),las=2,tck=0,cex.axis=0.6)
      #mtext(rownames(sig_scores))
      # }else if (Xmotif@Number == 1){
      # 	biclust::heatmapBC(stats::na.omit(t(sig_scores)),bicResult=Xmotif,number=1,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')
    }else{
      # print(paste0("Zero or one co-clusters found: ", names[i]))
      graphics::plot.new()
      graphics::title(paste0('\n\n\n <=1 bivariate clusters for\n',names_datasets[dataset_ind],' ',names_sigs[sig_ind]))
    }
    grab_grob()
  })

  draw.heatmaps(hmaps,names_datasets,names_sigs)

  grDevices::dev.copy(grDevices::pdf,file.path(out_dir,'sig_eval_bivariate_clustering_binarized_maps.pdf'),width=10,height=10)
  grDevices::dev.off()
  radar_plot_values
}
