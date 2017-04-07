eval_expr_loc <- function(gene_sigs_list,names_sigs, mRNA_expr_matrix,names_datasets, thresholds = NULL, out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values){
  #calculate the number of rows and columns in the image
  num_rows <- length(names_sigs)#ceiling(sqrt(length(names)))
  num_cols <- length(names_datasets)#ceiling(length(names)/num_rows)

  if (showResults){
    grDevices::dev.new()
  }else{
    grDevices::pdf(file.path(out_dir,'sig_expr_barcharts_NA_values.pdf'),width=10,height=10)

  }

  graphics::par(mfrow=c(num_rows,num_cols),cex=0.7, cex.axis=0.5)
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets)){
      #calculate the porportion of non-NA expression data in the matrix
      genes_expr <- mRNA_expr_matrix[[names_datasets[i]]][gene_sig,]
      gene_expr_vals <- (rowSums(is.na(genes_expr)) / dim(genes_expr)[2])
      gene_expr_vals <- -sort(-gene_expr_vals)
      if(max(gene_expr_vals)==0){
        bar_expr <- graphics::barplot(gene_expr_vals,xlab="Signature Gene IDs",ylab="Proportion of NA expression ",main=paste0("Signature gene expression across samples\n",names_datasets[i],' ',names_sigs[k]), axisnames=F,axis=F,ylim=c(0,1))
      }else{
        bar_expr <- graphics::barplot(gene_expr_vals,xlab="Signature Gene IDs",ylab="Proportion of NA expression ",main=paste0("Signature gene expression across samples\n",names_datasets[i],' ',names_sigs[k]), axisnames=F,axis=F)

      }
      graphics::text(bar_expr, graphics::par("usr")[3], labels = names(gene_expr_vals), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.5)
      graphics::axis(2)
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['med_prop_na'] <- stats::median(1-gene_expr_vals)

    }
  }
  if(showResults){
    grDevices::dev.copy(grDevices::pdf,file.path(out_dir,'sig_expr_barcharts_NA_values.pdf'),width=10,height=10)
  }
  grDevices::dev.off()

  if (showResults){
    grDevices::dev.new()
  }else{
    grDevices::pdf(file.path(out_dir,'sig_expr_barcharts.pdf'),width=10,height=10)
  }

  graphics::par(mfrow=c(num_rows,num_cols),cex=0.7, cex.axis=0.5)
  if (length(thresholds) == 0) {
    thresholds <- rep(0,length(names_datasets))
    for ( i in 1:length(names_datasets)){
      thresholds[i] <- stats::median(unlist(stats::na.omit(mRNA_expr_matrix[[names_datasets[i]]])))
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
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['med_prop_above_med'] <- stats::median(gene_expr_vals)
      bar_expr <- graphics::barplot(gene_expr_vals,xlab="Signature Gene IDs",ylab="Proportion with expression above threshold",main=paste0("Signature gene expression across samples\n",names_datasets[i],' ',names_sigs[k]), axisnames=F,axis=F)
      graphics::text(bar_expr, graphics::par("usr")[3], labels = names(gene_expr_vals), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.5)
      graphics::axis(2)
    }
  }
  if(showResults){
    grDevices::dev.copy(grDevices::pdf,file.path(out_dir,'sig_expr_barcharts.pdf'),width=10,height=10)
  }
  grDevices::dev.off()


  if (showResults){
    grDevices::dev.new()
  }else{
    grDevices::pdf(file.path(out_dir,'sig_expr_density_plots.pdf'),width=10,height=10)
  }

  graphics::par(mfrow=c(num_rows,num_cols))
  for (k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    for (i in 1:length(names_datasets)){
      #calculate the porportion of nonzero expression data in the matrix
      genes_expr <- mRNA_expr_matrix[[names_datasets[i]]][gene_sig,]
      gene_expr_vals <- 1 - (rowSums(genes_expr < thresholds[i]) / (dim(genes_expr)[2]))
      graphics::plot(stats::density(stats::na.omit(gene_expr_vals),adjust=0.25),main=paste0("Signature gene expression across samples\n",names_datasets[i], ' ',names_sigs[k]),ylab="Density")

    }
  }
  if(showResults){
    grDevices::dev.copy(grDevices::pdf,file.path(out_dir,'sig_expr_density_plots.pdf'),width=10,height=10)
  }
  grDevices::dev.off()
  cat('Expression and density graphs created successfully.\n', file=file)

  # print(paste0("Min expression occurs for gene ID ",names(gene_expr_vals)[1], ", with ", round(gene_expr_vals[1]*100,2),"% non-zero expression."))
  radar_plot_values
}
