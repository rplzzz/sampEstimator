#' Bounded secant solver for a 1-D function of integer values
#' 
#' This solver is specialized to the case where the function being solved is
#' only defined for integer values.  We still pass noninteger values to the
#' function, but we expect that they will be rounded.  This affects our stopping
#' criteria.
#' 
#' @param x1 lower bound
#' @param x2 upper bound
#' @param targfun function to solve
#' @keywords internal
bssolv <- function(x1, x2, targfun)
{
  fx1 <- targfun(x1)
  fx2 <- targfun(x2)

  while(fx1*fx2 > 0) {
    x2 <- x2*1.25
    fx2 <- targfun(x2)
  }
  
  ## Two stopping conditions.  
  while(x2-x1 > 0.25 && floor(x1) != floor(x2)) {
    m <- (fx2-fx1)/(x2-x1)
    xmid <- x1 - fx1/m
    if(xmid <= x1 || xmid >= x2) {
      xmid <- 0.5*(x1+x2)
    }
    fxmid <- targfun(xmid)
    
    if(fxmid*fx1 > 0) {
      x1 <- xmid
      fx1 <- fxmid
    }
    else {
      x2 <- xmid
      fx2 <- fxmid
    }
  }
  ceiling(x2)
}
