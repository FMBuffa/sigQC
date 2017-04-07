eval_compactness_loc <- function(gene_sigs_list,names_sigs, mRNA_expr_matrix, names_datasets, out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values){
  # require(gplots)
  if (showResults){
    grDevices::dev.new()
  }
  # pdf(file.path(out_dir,'sig_autocor_hmps.pdf'),width=10,height=10)

  graphics::par(cex.main=0.8,cex.lab = 0.6,oma=c(2,0,0,0),mar=c(0,0,0,0))
  hmaps <- lapply(1:(length(names_sigs) *length(names_datasets)),function(i) {
    dataset_ind <- i %% length(names_datasets)
    if (dataset_ind == 0 ){
      dataset_ind <- length(names_datasets)
    }
    sig_ind <- ceiling(i/length(names_datasets))
    gene_sig <- gene_sigs_list[[names_sigs[sig_ind]]]
    autocors <- stats::cor(t(stats::na.omit(mRNA_expr_matrix[[names_datasets[dataset_ind]]][intersect(gene_sig,rownames(mRNA_expr_matrix[[names_datasets[dataset_ind]]])),])),method='spearman')

    tryCatch({
      gplots::heatmap.2( stats::na.omit(autocors),
                         col = gplots::colorpanel(100,"blue","white","red"), #redgreen(100),#colorpanel(100,"red","yellow","green"),
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
  grDevices::dev.copy(grDevices::pdf,file.path(out_dir,'sig_autocor_hmps.pdf'),width=10,height=10)
  grDevices::dev.off()

  if (showResults){
    grDevices::dev.new()
  }else{
    grDevices::pdf(file.path(out_dir,'sig_autocor_dens.pdf'),width=10,height=10)
  }

  graphics::par(cex.main=0.8,cex.lab = 0.6,oma=c(2,2,2,2),mar=c(2,2,2,2),mfrow=c(1,1))
  max_dens <- -9999
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets) ){
      autocors <- stats::cor(stats::na.omit(t(mRNA_expr_matrix[[names_datasets[i]]][intersect(gene_sig,rownames(mRNA_expr_matrix[[names_datasets[i]]])),])),method='spearman')
      cur_max <- max(stats::density(unlist(stats::na.omit(autocors)))$y)
      if (max_dens  < cur_max){
        max_dens <- cur_max
      }
    }
  }

  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets) ){
      autocors <- stats::cor(t(stats::na.omit(mRNA_expr_matrix[[names_datasets[i]]][gene_sig,])),method='spearman')
      if ((i ==1) && (k==1)){
        graphics::plot(stats::density(unlist(stats::na.omit(autocors))),ylim=c(0,ceiling(max_dens)),col=i,main=NA,lwd=2,lty=k)
      }else{
        graphics::lines(stats::density(unlist(stats::na.omit(autocors))),col=i,main=NA,lwd=2,lty=k)
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
    grDevices::dev.copy(grDevices::pdf,file.path(out_dir,'sig_autocor_dens.pdf'),width=10,height=10)
  }
  grDevices::dev.off()
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
    if (length(names_datasets) > 1){
      overall_rank_mat <- matrix(NA,nrow=length(unique(gene_sig)),ncol=length(names_datasets))
      row.names(overall_rank_mat) <- unique(gene_sig)
      colnames(overall_rank_mat) <- names_datasets
      for (i in 1:length(names_datasets)){
        autocors <- stats::cor(t(stats::na.omit(mRNA_expr_matrix[[names_datasets[i]]][intersect(unique(gene_sig),rownames(mRNA_expr_matrix[[names_datasets[i]]])),])),method='spearman')
        median_scores <- apply(autocors,2,function(x) {stats::median(stats::na.omit(x))})
        overall_rank_mat[names(median_scores),i] <- median_scores
      }

      # require(RankProd)
      RP.out <-RankProd::RP(data = overall_rank_mat,cl = rep(1,times=length(names_datasets)),logged = F,gene.names=rownames(overall_rank_mat))
      RankProd::plotRP(RP.out ,cutoff=0.05)
      RankProd::topGene(RP.out,cutoff=0.05,method="pfp",logged=F, gene.names=intersect(gene_sig,rownames(mRNA_expr_matrix[[names_datasets[i]]])))
      cat("Autocorrelation rank product successfully computed.\n", file=file)

    }else{
      cat("Rank product not computed as there is only one dataset.\n", file=file)

    }
    if(showResults){
      grDevices::dev.copy(grDevices::pdf,file.path(out_dir,paste0('sig_autocor_rankProd_',names_sigs[k],'.pdf')),width=10,height=10)
    }
    grDevices::dev.off()
  }
  radar_plot_values

}
