#' eval_stan_loc.R
#'
#' This function creates the plots of effects of standardisation. It checks the correlation between the median and the z-transformed
#' median of the signature genes in the dataset. 
#' @param gene_sigs_list A list of genes representing the gene signature to be tested.
#' @param names_sigs The names of the gene signatures (one name per gene signature, in gene_sigs_list)
#' @param mRNA_expr_matrix A list of expression matrices
#' @param names_datasets The names of the different datasets contained in mRNA_expr_matrix
#' @param out_dir A path to the directory where the resulting output files are written
#' @param file File representing the log file where errors can be written
#' @param showResults Tells if open dialog boxes showing the computed results. Default is FALSE
#' @param radar_plot_values A list of values that store computations that will be used in the final summary radarplot
#' @keywords eval_stan_loc

eval_stan_loc <- function(gene_sigs_list, names_sigs,mRNA_expr_matrix, names_datasets,out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values ){
  # defining the size of the plotting area
  num_rows <- length(names_sigs)#ceiling(sqrt(length(names)))
  num_cols <- length(names_datasets)#ceiling(length(names)/num_rows)

  if (showResults){
    grDevices::dev.new()
  }else{
    grDevices::pdf(file.path(out_dir, 'sig_standardisation_comp.pdf'),width=10,height=10)
  }

  graphics::par(mfrow=c(num_rows,num_cols),oma=c(2,2,2,2),mar=c(4,4,4,4))
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets)){

      z_transf_mRNA <- mRNA_expr_matrix[[names_datasets[i]]][gene_sig,]
      for (j in 1:length(gene_sig)){
        z_transf_mRNA[gene_sig[j],] <- (as.numeric(z_transf_mRNA[gene_sig[j],]) - mean(as.numeric(z_transf_mRNA[gene_sig[j],]),na.rm=T)) / stats::sd(as.numeric(z_transf_mRNA[gene_sig[j],]),na.rm=T)
      }
      z_transf_scores <- apply(z_transf_mRNA[gene_sig,],2,function(x) {stats::median(stats::na.omit(x))})

      med_scores <- apply(mRNA_expr_matrix[[names_datasets[i]]][gene_sig,],2,function(x){ stats::median(stats::na.omit(x))})

      combined_scores <- as.data.frame(cbind(med_scores,z_transf_scores))
      colnames(combined_scores) <- c('Median',"Z_Median")

      jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      graphics::smoothScatter(med_scores,z_transf_scores,colramp = jet.colors,xlab=NA,ylab=NA,main='Median vs Z-median')
      graphics::points(med_scores,z_transf_scores,pch='.')
      graphics::par(new=T,srt=45)
      graphics::plot(stats::density(med_scores), axes=F, xlab=NA, ylab=NA,  col='red',main=NA)
      graphics::axis(side = 4)
      graphics::mtext(side = 4, line = 2, 'Density',cex=0.8)
      graphics::mtext(side = 2, line = 2, 'Z-median',cex=0.8)
      graphics::mtext(side = 1, line = 2, 'Median',cex=0.8)
      graphics::mtext(side=3,line=2.5,paste0(names_datasets[i],', ',names_sigs[k] ))

      rho <- stats::cor(med_scores,z_transf_scores,method='spearman')
      rho_mean_med <- rho
      graphics::mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(med_scores))
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['standardization_comp'] <- rho
    }
  }
  cat('Standardisation compared successfully.\n', file=file)
  if(showResults){
    grDevices::dev.copy(grDevices::pdf,file.path(out_dir, 'sig_standardisation_comp.pdf'),width=10,height=10)
  }
  grDevices::dev.off()
  radar_plot_values
}
