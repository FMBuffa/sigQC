# eval_compactness_loc.R
#
# This function creates the plots of autocorrelation, as well as the rank product computation. Specifically, it creates
# plots for heatmap of autocorrelation, density plot of autocorrelation values, and if there is more than one dataset, for
# each individual gene signature, will do a rank product analysis on the lists of median autocorrelation coefficients to
# identify genes consistently low ranking in autocorrelation among more than one dataset.
# @param gene_sigs_list A list of genes representing the gene signature to be tested.
# @param names_sigs The names of the gene signatures (one name per gene signature, in gene_sigs_list)
# @param mRNA_expr_matrix A list of expression matrices
# @param names_datasets The names of the different datasets contained in mRNA_expr_matrix
# @param out_dir A path to the directory where the resulting output files are written
# @param file File representing the log file where errors can be written
# @param showResults Tells if open dialog boxes showing the computed results. Default is FALSE
# @param radar_plot_values A list of values that store computations that will be used in the final summary radarplot
# @keywords eval_compactness_loc

eval_compactness_loc <- function(gene_sigs_list,names_sigs, mRNA_expr_matrix, names_datasets, out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values){
  # require(gplots)
  if (showResults){
    grDevices::dev.new()
  }
  # else{
  #   grDevices::pdf(file.path(out_dir, 'sig_autocor_hmps.pdf'),width=10,height=10)
  # }
  # pdf(file.path(out_dir,'sig_autocor_hmps.pdf'),width=10,height=10)

  graphics::par(cex.main=0.8,cex.lab = 0.6,oma=c(2,0,0,0),mar=c(0,0,0,0))
  hmaps <- lapply(1:(length(names_sigs) *length(names_datasets)),function(i) {
    dataset_ind <- i %% length(names_datasets)
    if (dataset_ind == 0 ){
      dataset_ind <- length(names_datasets)
    }
    sig_ind <- ceiling(i/length(names_datasets))
    gene_sig <- gene_sigs_list[[names_sigs[sig_ind]]]
    data.matrix = mRNA_expr_matrix[[names_datasets[dataset_ind]]]
    inter <- intersect(gene_sig[,1], row.names(data.matrix))
    autocors <- stats::cor(t(stats::na.omit(data.matrix[inter,])),method='spearman')

    tryCatch({
      gplots::heatmap.2( stats::na.omit(autocors),
                         col = gplots::colorpanel(100,"blue","white","red"), #redgreen(100),#colorpanel(100,"red","yellow","green"),
                         trace = "none",
                         xlab = "Gene ID",
                         ylab="Gene ID",
                         na.color="grey",
                         labRow=rownames(autocors),
                         labCol=colnames(autocors),#gene_sig,
                         main = paste0("\n\nAutocorrelation\n", names_datasets[[dataset_ind]] ,' ',names_sigs[[sig_ind]]),
                         dendrogram = "col",
                         symbreaks = T,
                         Rowv = T,Colv=T ,key.xlab='Rho',key.ylab=NA,  key.title=NA,cexRow=0.5,cexCol=0.5,margins=c(4,4))
    },
    error=function(err){
      cat(paste0("There was an error, likely due to NA values in data, for dataset: ",names_datasets[dataset_ind]," ", names_sigs[sig_ind]," ", err,'\n'), file=file)

    })
    grab_grob()
  })
  draw.heatmaps(hmaps,names_datasets,names_sigs)
  # if(showResults){
  grDevices::dev.copy(grDevices::pdf,file.path(out_dir,'sig_autocor_hmps.pdf'),width=4*(length(names_datasets)),height=4*(length(names_sigs)))#width=10,height=10)
  # }
  grDevices::dev.off()

  if (showResults){
    grDevices::dev.new()
  }else{
    grDevices::pdf(file.path(out_dir,'sig_autocor_dens.pdf'),width=5,height=5)
  }

  graphics::par(cex.main=0.8,cex.lab = 0.6,oma=c(2,2,2,2),mar=c(2,2,2,2),mfrow=c(1,1))
  max_dens <- -9999
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets) ){
      data.matrix = mRNA_expr_matrix[[names_datasets[i]]]
      inter <- intersect(gene_sig[,1], row.names(data.matrix))
      autocors <- stats::cor(stats::na.omit(t(data.matrix[inter,])),method='spearman')
      cur_max <- max(stats::density(unlist(stats::na.omit(autocors)))$y)

      if (max_dens  < cur_max){
        max_dens <- cur_max
      }
    }
  }

  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets) ){
      data.matrix = mRNA_expr_matrix[[names_datasets[i]]]
      inter <- intersect(gene_sig[,1], row.names(data.matrix))
      autocors <- stats::cor(t(stats::na.omit(data.matrix[inter,])),method='spearman')
      if ((i ==1) && (k==1)){
        graphics::plot(stats::density(unlist(stats::na.omit(autocors))),ylim=c(0,ceiling(max_dens)),xlim=c(-1,1),col=i,main=NA,lwd=2,lty=k)
      }else{
        graphics::lines(stats::density(unlist(stats::na.omit(autocors))),ylim=c(0,ceiling(max_dens)),xlim=c(-1,1),col=i,main=NA,lwd=2,lty=k)
      }
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['autocor_median'] <- stats::median(stats::na.omit(autocors))

    }
  }

  graphics::mtext(side = 2, line = 2, 'Density',cex=0.8)
  graphics::mtext(side = 1, line = 2, 'Rho',cex=0.8)
  graphics::mtext(side=3,line=2,'Autocorrelation Density')

  op <- graphics::par(cex=0.6)
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
  graphics::legend("topright",legend_names,col=legend_cols,lty=legend_lty,lwd=rep(2,times=(length(names_datasets) * length(names_sigs))))

  if(showResults){
    grDevices::dev.copy(grDevices::pdf,file.path(out_dir,'sig_autocor_dens.pdf'),width=5,height=5)
  }
  grDevices::dev.off()
  if (length(names_datasets) > 1){

  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (showResults){
      grDevices::dev.new()
    }else{
      grDevices::pdf(file.path(out_dir,paste0('sig_autocor_rankProd_',names_sigs[k],'.pdf')),width=10,height=10)
    }

    graphics::par(cex.main=0.8,cex.lab = 0.6,oma=c(2,2,2,2),mar=c(4,4,4,4))
    cat("Autocorrelation metrics successfully computed.\n", file=file)

    #now we take the median of the genes' autocorrelation for each gene in each dataset and then look at the rank product over the different cancer types
    #note that the rank product analysis is only done if there is more than one dataset (otherwise not done, and is doen separately for each gene signature)
      overall_rank_mat <- matrix(NA,nrow=length(unique(gene_sig[,1])),ncol=length(names_datasets))
      row.names(overall_rank_mat) <- unique(gene_sig[,1])
      colnames(overall_rank_mat) <- names_datasets
      for (i in 1:length(names_datasets)){
        data.matrix = mRNA_expr_matrix[[names_datasets[i]]]
        inter = intersect(unique(gene_sig[,1]),rownames(data.matrix))
        autocors <- stats::cor(t(stats::na.omit(data.matrix[inter,])),method='spearman')
        median_scores <- as.matrix(apply(autocors,2,function(x) {stats::median(stats::na.omit(x))}))
        overall_rank_mat[rownames(median_scores),i] <- median_scores[,1]
      }

      # require(RankProd)
      RP.out <-RankProd::RP(data = overall_rank_mat,cl = rep(1,times=length(names_datasets)),logged = F,gene.names=rownames(overall_rank_mat))
      RankProd::plotRP(RP.out ,cutoff=0.05)
     
      table_rank_prod <- RankProd::topGene(RP.out,cutoff=0.05,method="pfp",logged=F, gene.names=rownames(overall_rank_mat))#intersect(gene_sig[,1],rownames(mRNA_expr_matrix[[names_datasets[i]]])))
      # output the rank product table to file
      if( (!is.null(table_rank_prod$Table1)) & (!is.null(table_rank_prod$Table2))){
          dir.create(file.path(out_dir,'rank_prod'))
          utils::write.csv(table_rank_prod$Table1,file=file.path(out_dir, 'rank_prod',paste0('rank_product_table1_',names_sigs[k],'.txt')),quote=F,sep='\t')
          utils::write.csv(table_rank_prod$Table2,file=file.path(out_dir, 'rank_prod',paste0('rank_product_table2_',names_sigs[k],'.txt')),quote=F,sep='\t')
      }
      
      cat("Autocorrelation rank product successfully computed.\n", file=file)

    
    if(showResults){
      grDevices::dev.copy(grDevices::pdf,file.path(out_dir,paste0('sig_autocor_rankProd_',names_sigs[k],'.pdf')),width=10,height=10)
    }
    grDevices::dev.off()
  }
  }else{
      cat("Rank product not computed as there is only one dataset.\n", file=file)

    }

  radar_plot_values

}
