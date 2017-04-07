# draw.heatmaps.R
#
# Makes the plots in each of the subfunctions.
# @param hmaps A list of genes representing the gene signature to be tested.
# @param names_datasets A list of expression matrices
#
draw.heatmaps <- function(hmaps,names_datasets,names_sigs){
  grid::grid.newpage()
  num_rows <- length(names_sigs)#ceiling(sqrt(length(names)))
  num_cols <- length(names_datasets)#ceiling(length(names)/num_rows)

  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = num_rows, ncol=num_cols)))
  count <- 1
  for (i in 1:num_rows){
    for(j in 1:num_cols){
      grid::grid.draw(grid::editGrob(hmaps[[count]], vp=grid::viewport(layout.pos.row = i, layout.pos.col = j , clip=T)))
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
