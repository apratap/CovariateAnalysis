edge.lfdr <- function(p, trunc=TRUE, monotone=TRUE, transf=c("probit", "logit"), adj=1.5, eps=10^-8, lambda=0.8, ...) {
  pi0 <- mean(p >= lambda)/(1 - lambda)
  pi0 <- min(pi0, 1)
  
  n = length(p)
  transf = match.arg(transf)
  
  if(transf=="probit") {
    p = pmax(p, eps)
    p = pmin(p, 1-eps)
    x = qnorm(p)
    myd = density(x, adjust=adj)
    mys = smooth.spline(x=myd$x, y=myd$y)
    y = predict(mys, x)$y
    lfdr = pi0*dnorm(x)/y
  }
  
  if(transf=="logit") {
    x = log((p+eps)/(1-p+eps))
    myd = density(x, adjust=adj)
    mys = smooth.spline(x=myd$x, y=myd$y)
    y = predict(mys, x)$y
    dx = exp(x)/(1+exp(x))^2
    lfdr = pi0*dx/y
  }
  
  if(trunc) {lfdr[lfdr > 1] = 1}
  if(monotone) {	
    lfdr = lfdr[order(p)]
    lfdr = mono(lfdr)
    lfdr = lfdr[rank(p)]
  }
  
  return(lfdr)
}