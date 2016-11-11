mono <- function(lfdr){
  library(sva)
  .Call("monotone",as.numeric(lfdr),PACKAGE="sva")
}