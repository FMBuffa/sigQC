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
    grDevices::dev.new()
  }else{
    grDevices::pdf(file.path(out_dir,'sig_radarplot.pdf'),width=10,height=10)
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
  graphics::layout(rbind(1,2), heights=c(7,1))  # put legend on bottom 1/8th of the chart

  fmsb::radarchart(as.data.frame(radar_plot_mat),
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
  graphics::par(mar=c(0, 0, 0, 0))

  graphics::plot.new()
  graphics::par(xpd=TRUE)

  graphics::legend('center', legend=legend_labels, seg.len=0.5, title="Datasets",lty=legend_lty,pch=1,
                   bty="n" ,lwd=1, col=legend_cols,cex=0.6,horiz=T)
  if(showResults){
    grDevices::dev.copy(grDevices::pdf,file.path(out_dir,'sig_radarplot.pdf'),width=10,height=10)
  }
  grDevices::dev.off()
  cat('Radar chart made successfully.\n', file=file)

}
