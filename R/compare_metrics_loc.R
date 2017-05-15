# compare_metrics_loc.R
#
# This function creates the plots to compare various signature summary statistic metrics against each other
# That is, compares the mean, median and first principal component with each other, and also produces PCA vs variance
# plot for the first 10 principal components of the dataset
# @param gene_sigs_list A list of genes representing the gene signature to be tested.
# @param names_sigs The names of the gene signatures (one name per gene signature, in gene_sigs_list)
# @param mRNA_expr_matrix A list of expression matrices
# @param names_datasets The names of the different datasets contained in mRNA_expr_matrix
# @param out_dir A path to the directory where the resulting output files are written
# @param file File representing the log file where errors can be written
# @param showResults Tells if open dialog boxes showing the computed results. Default is FALSE
# @param radar_plot_values A list of values that store computations that will be used in the final summary radarplot
# @keywords compare_metrics_loc

compare_metrics_loc <- function(gene_sigs_list,names_sigs, mRNA_expr_matrix, names_datasets, out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values){
  # require(gplots)
  dir.create(file.path(out_dir,'metrics_tables'))
  for(k in 1:length(names_sigs)){ 
    #for each signature we will make a separate file comparing the datasets
    gene_sig <- gene_sigs_list[[names_sigs[k]]] # load the gene signature
    #set up canvas for plotting
    if (showResults){
      grDevices::dev.new()
    }else{
      grDevices::pdf(file.path(out_dir,paste0('sig_compare_metrics_',names_sigs[k],'.pdf')),width=3*length(names_datasets),height=10)
    }

    #find the max number of characters in the title
    max_title_length <- -999
    for( i in 1:length(names_datasets)){
      if(max_title_length < nchar(paste0(names_datasets[i],' ',names_sigs[k]))){
        max_title_length <- nchar(paste0(names_datasets[i],' ',names_sigs[k]))
      }
    }
    #set up the canvas
    graphics::par(mfcol = c(4,length(names_datasets)),mar=c(4,4,4,4))

    for ( i in 1:length(names_datasets)){
      # now we can loop over the datasets for the plot and generate the metrics for every dataset with this signature
      data.matrix = mRNA_expr_matrix[[names_datasets[i]]] #load the data
      data.matrix[!(is.finite(as.matrix(data.matrix)))] <- NA #ensure that the data is not infintie
      inter = intersect(gene_sig[,1],rownames(data.matrix)) #consider only the genes actually present in the data

      med_scores <- apply(data.matrix[inter,],2,function(x){stats::median(stats::na.omit(x))}) #compute median
      mean_scores <- apply(data.matrix[inter,],2,function(x){mean(stats::na.omit(x))}) #compute mean
      
      pca1_scores <- NULL

      tryCatch({
        pca1_scores <- stats::prcomp(stats::na.omit(t(data.matrix[inter,])),retx=T) #compute PCA1
        
        },error=function(e){
          pca1_scores <<- NULL
          cat(paste0("There was an error:  ",names_datasets[i]," ", names_sigs[k]," ", e,'\n'), file=file)

     #     print(paste0("error: ", e))
      })


      # print(paste0("test ",pca1_scores ))
      if(length(pca1_scores) > 1){#!is.null(pca1_scores)){
        props_of_variances <- pca1_scores$sdev^2/(sum(pca1_scores$sdev^2)) #for the scree plot
        pca1_scores <- pca1_scores$x[,1] #takes the actual PCA1 scores
        common_score_cols <- intersect(names(med_scores),intersect(names(mean_scores),names(pca1_scores))) #ensures we have the same samples for each plot
      }else{
        common_score_cols <- intersect(names(med_scores),names(mean_scores))
      }

      if(length(common_score_cols) > 1){
        #the following is the colourmap for the 2D scatter
        jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        #plotting commands for the 2D scatterplot for mean-median correlation as follows
        graphics::smoothScatter(med_scores[common_score_cols],mean_scores[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='Median vs Mean')
        graphics::points(med_scores[common_score_cols],mean_scores[common_score_cols],pch='.') #draw the points on top of it
        graphics::par(new=T)#,srt=45)
        graphics::plot(stats::density(med_scores[common_score_cols]), axes=F, xlab=NA, ylab=NA,  col='red',main=NA,lwd=2) #draw the density plot behind it
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density',cex=0.8) #labels for the axes
        graphics::mtext(side = 2, line = 2, 'Mean',cex=0.8)
        graphics::mtext(side = 1, line = 2, 'Median',cex=0.8)
        graphics::mtext(side=3,line=2.5,paste0(names_datasets[i],' ',names_sigs[k]),cex=min(1,3*10/max_title_length)) #title
        rho <- stats::cor(med_scores[common_score_cols],mean_scores[common_score_cols],method='spearman') 
        rho_mean_med <- rho
        graphics::mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(med_scores[common_score_cols]))
      }else{
        graphics::plot.new()
        graphics::mtext(side=3,line=2.5,paste0(names_datasets[i],' ',names_sigs[k]),cex=min(1,3*10/max_title_length)) #title

        graphics::title(paste0('\n\nToo many NA values for Mean/Median in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
        rho_mean_med <- 0
      }
      if (length(pca1_scores) > 1){#(!is.null(pca1_scores)){
        #plotting for the mean-pca1
        graphics::smoothScatter(mean_scores[common_score_cols],pca1_scores[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='Mean vs PCA1')
        graphics::points(mean_scores[common_score_cols],pca1_scores[common_score_cols],pch='.')
        graphics::par(new=T)
        graphics::plot(stats::density(mean_scores[common_score_cols]), axes=F, xlab=NA, ylab=NA, col='red',main=NA,lwd=2)
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density',cex=0.8)
        rho <- stats::cor(mean_scores[common_score_cols],pca1_scores[common_score_cols],method='spearman')
        rho_mean_pca1 <- rho
        graphics::mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(mean_scores[common_score_cols]))
        graphics::mtext(side = 2, line = 2, 'PCA1',cex=0.8)
        graphics::mtext(side = 1, line = 2, 'Mean',cex=0.8)

        #plotting for the pca1-median
        graphics::smoothScatter(pca1_scores[common_score_cols],med_scores[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='PCA1 vs Median')
        graphics::points(pca1_scores[common_score_cols],med_scores[common_score_cols],pch='.')
        graphics::par(new=T)
        graphics::plot(stats::density(pca1_scores[common_score_cols]), axes=F, xlab=NA, ylab=NA, col='red',main=NA,lwd=2)
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density',cex=0.8)
        rho <- stats::cor(pca1_scores[common_score_cols],med_scores[common_score_cols],method='spearman')
        rho_pca1_med <- rho
        graphics::mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(pca1_scores[common_score_cols]))
        graphics::mtext(side = 2, line = 2, 'Median',cex=0.8)
        graphics::mtext(side = 1, line = 2, 'PCA1',cex=0.8)
      }else{

#        graphics::par(new=T)
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for PCA1/Mean in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
        rho_mean_pca1 <- 0

 #       graphics::par(new=T)
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for Median/PCA1 in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
        rho_pca1_med <- 0

      }
      #here we put the cor plot of the mean median and pca1 (if we wanted an autocorrelation heatmap of these metrics)

      # autocors_mat <- matrix(0,nrow=3,ncol=3)
      # row.names(autocors_mat) <- c('Mean' , 'Median','PCA1')
      # colnames(autocors_mat) <- c('Mean' , 'Median','PCA1')
      # autocors_mat[1,1] <- 1
      # autocors_mat[1,2] <- rho_mean_med
      # autocors_mat[1,3] <- rho_mean_pca1

      # autocors_mat[2,1] <- rho_mean_med
      # autocors_mat[2,2] <- 1
      # autocors_mat[2,3] <- rho_pca1_med

      # autocors_mat[3,1] <- rho_mean_pca1
      # autocors_mat[3,2] <- rho_pca1_med
      # autocors_mat[3,3] <- 1

      #here we output the table of mean, median and pca1 scores for each sample to a table
      if(length(pca1_scores) > 1){#(!is.null(pca1_scores)){
        output_mat <- cbind(mean_scores[common_score_cols],cbind(med_scores[common_score_cols],pca1_scores[common_score_cols]))
        colnames(output_mat) <- c('Mean Scores','Median Scores','PCA1 Scores')
      }else{
        output_mat <- cbind(mean_scores[common_score_cols],med_scores[common_score_cols])
        colnames(output_mat) <- c('Mean Scores','Median Scores')

      }

      utils::write.table(output_mat,file=file.path(out_dir,'metrics_tables', paste0('metrics_table_',names_sigs[k],'_',names_datasets[i],'.txt')),quote=F,sep='\t')

      #stores values that will be used in the radarplot
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['rho_mean_med'] <- rho_mean_med
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['rho_pca1_med'] <- rho_pca1_med
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['rho_mean_pca1'] <- rho_mean_pca1
      
      if(length(pca1_scores) > 1){#(!is.null(pca1_scores)){
        #draws the scree plot
        bars_plot <- props_of_variances[1:min(10,length(props_of_variances))]
        graphics::barplot(bars_plot,main="PCA vs proportion\n of variance") #ylim= c(0,1),
        graphics::mtext(side = 1, line = 2, 'PCA',cex=0.8)
      }else{
       # graphics::par(new=T)
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for PCA1 in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
      }
    }
    #saves file
    if(showResults){
      grDevices::dev.copy(grDevices::pdf,file.path(out_dir,paste0('sig_compare_metrics_',names_sigs[k],'.pdf')),width=3*length(names_datasets),height=10)
    }
    grDevices::dev.off()
  }
  cat('Metrics compared successfully.\n', file=file) #output to log
  radar_plot_values #returns the radarplot values
}
