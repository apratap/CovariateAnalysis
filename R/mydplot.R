library(ggdendro)
library(ggplot2)
library(reshape)
library(grid)

mydplot <- function(ddata, row=!col, col=!row, labels=col) {
  ## plot a dendrogram
  yrange <- range(ddata$segments$y)
  yd <- yrange[2] - yrange[1]
  nc <- max(nchar(as.character(ddata$labels$label)))
  tangle <- if(row) { 0 } else { 90 }
  tshow <- col
  p <- ggplot() +
    geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = NULL, y = NULL) + theme_dendro()
  if(row) {
    p <- p +
      scale_x_continuous(expand=c(0.5/length(ddata$labels$x),0)) +
      coord_flip()
  } else {
    p <- p +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  return(p)
}

