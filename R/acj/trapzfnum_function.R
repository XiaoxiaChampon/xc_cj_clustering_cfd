
# library("pracma")

#Function to find the L2 distance between two latent curves
#' @param yy 1D vector (true curve)
#' @param yy2 1D vector (estimated curve)
#' @return scalar-L2 distance
trapzfnum <- function(yy,yy2,timestamps01)
{
  st=timestamps01[1]
  et=tail(timestamps01,n=1)
  x=seq(st,et,length=5000)
  xx=seq(st,et,length=length(yy))
  y1 <- cubicspline(xx, yy,x)
  y2 <- cubicspline(xx, yy2,x)
  out=sqrt(trapz(x, (y1-y2)^2) )
  return(out)
}

#'Function to find the Hellinger distance between two probability curves
#' @param yy 1D vector (true curve)
#' @param yy2 1D vector (estimated curve)
#' @param timestamps01 1D vector 
#' @return scalar-Hellinger distance
trapzfnump <- function(yy,yy2,timestamps01)
{
  st=timestamps01[1]
  et=tail(timestamps01,n=1)
  x=seq(st,et,length=5000)
  xx=seq(st,et,length=length(yy))
  y1 <- cubicspline(xx, yy,x)
  y1[y1<0]=0
  y2 <- cubicspline(xx, yy2,x)
  y2[y2<0]=0
  out=sqrt(trapz(x, (sqrt(y1)-sqrt(y2))^2) )
  return(out)
}
