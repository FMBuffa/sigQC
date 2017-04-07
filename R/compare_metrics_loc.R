compare_metrics_loc <- function(gene_sigs_list,names_sigs, mRNA_expr_matrix, names_datasets, out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values){
  # require(gplots)
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    hmaps <- list()
    if (showResults){
      grDevices::dev.new()
    }else{
      grDevices::pdf(file.path(out_dir,paste0('sig_compare_metrics_',names_sigs[k],'.pdf')),width=10,height=10)
    }

    graphics::par(mfcol = c(4,length(names_datasets)),mar=c(4,4,4,4))
    for ( i in 1:length(names_datasets)){
      mRNA_expr_matrix[[names_datasets[i]]][!(is.finite(as.matrix(mRNA_expr_matrix[[names_datasets[i]]])))] <- NA

      med_scores <- apply(mRNA_expr_matrix[[names_datasets[i]]][gene_sig,],2,function(x){stats::median(stats::na.omit(x))})
      mean_scores <- apply(mRNA_expr_matrix[[names_datasets[i]]][gene_sig,],2,function(x){mean(stats::na.omit(x))})
      pca1_scores <- stats::prcomp(stats::na.omit(t(mRNA_expr_matrix[[names_datasets[i]]][gene_sig,])),retx=T)
      props_of_variances <- pca1_scores$sdev^2/(sum(pca1_scores$sdev^2))
      pca1_scores <- pca1_scores$x[,1]

      common_score_cols <- intersect(names(med_scores),intersect(names(mean_scores),names(pca1_scores)))
      jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      graphics::smoothScatter(med_scores[common_score_cols],mean_scores[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='Median vs Mean')
      graphics::points(med_scores[common_score_cols],mean_scores[common_score_cols],pch='.')
      graphics::par(new=T)#,srt=45)
      graphics::plot(stats::density(med_scores[common_score_cols]), axes=F, xlab=NA, ylab=NA,  col='red',main=NA)
      graphics::axis(side = 4)
      graphics::mtext(side = 4, line = 2, 'Density',cex=0.8)
      graphics::mtext(side = 2, line = 2, 'Mean',cex=0.8)
      graphics::mtext(side = 1, line = 2, 'Median',cex=0.8)
      graphics::mtext(side=3,line=2.5,paste0(names_datasets[i],' ',names_sigs[k]))

      rho <- stats::cor(med_scores[common_score_cols],mean_scores[common_score_cols],method='spearman')
      rho_mean_med <- rho
      graphics::mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(med_scores[common_score_cols]))
      graphics::smoothScatter(mean_scores[common_score_cols],pca1_scores[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='Mean vs PCA1')
      graphics::points(mean_scores[common_score_cols],pca1_scores[common_score_cols],pch='.')
      graphics::par(new=T)
      graphics::plot(stats::density(mean_scores[common_score_cols]), axes=F, xlab=NA, ylab=NA, col='red',main=NA)
      graphics::axis(side = 4)
      graphics::mtext(side = 4, line = 2, 'Density',cex=0.8)
      rho <- stats::cor(mean_scores[common_score_cols],pca1_scores[common_score_cols],method='spearman')
      rho_mean_pca1 <- rho
      graphics::mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(mean_scores[common_score_cols]))
      graphics::mtext(side = 2, line = 2, 'PCA1',cex=0.8)
      graphics::mtext(side = 1, line = 2, 'Mean',cex=0.8)

      graphics::smoothScatter(pca1_scores[common_score_cols],med_scores[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='PCA1 vs Median')
      graphics::points(pca1_scores[common_score_cols],med_scores[common_score_cols],pch='.')
      graphics::par(new=T)
      graphics::plot(stats::density(pca1_scores[common_score_cols]), axes=F, xlab=NA, ylab=NA, col='red',main=NA)
      graphics::axis(side = 4)
      graphics::mtext(side = 4, line = 2, 'Density',cex=0.8)
      rho <- stats::cor(pca1_scores[common_score_cols],med_scores[common_score_cols],method='spearman')
      rho_pca1_med <- rho
      graphics::mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(pca1_scores[common_score_cols]))
      graphics::mtext(side = 2, line = 2, 'Median',cex=0.8)
      graphics::mtext(side = 1, line = 2, 'PCA1',cex=0.8)

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


      # hmaps[[i]] <- gplots::heatmap.2(autocors_mat,
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
      graphics::barplot(bars_plot,main="PCA vs proportion\n of variance") #ylim= c(0,1),
      graphics::mtext(side = 1, line = 2, 'PCA',cex=0.8)
    }
    if(showResults){
      grDevices::dev.copy(grDevices::pdf,file.path(out_dir,paste0('sig_compare_metrics_',names_sigs[k],'.pdf')),width=10,height=10)
    }
    grDevices::dev.off()
  }
  cat('Metrics compared successfully.\n', file=file)
  radar_plot_values
}
