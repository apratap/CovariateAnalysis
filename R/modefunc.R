modefunc <- function(x) {
  return(as.numeric(names(sort(-table(x)))[1]))
}