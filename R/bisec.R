#' @title Bisection algorithm in Beta distribution
#' 
#' @description Running bisection algorithm to search for a2, the minimizer of (log((a2)^dea1*(1-a2)^dea2-delta))^2
#'
#' @param delta Test statistic in Han et al. (2012), 
#'               delta = (ttot1/(ttot1+ttot2))^dea1*(ttot2/(ttot1+ttot2))^dea2;  
#' @param dea1 first parameter in Beta distribution (number of events from the first arm)
#' @param dea2 second parameter in Beta distribution (number of events from the second arm)
#' @param upbd upper bound of a2
#' @param lowbd lower bound of a2
#' 
#' 
#' @return
#' a2
#' 
#' @export
#'
#' @examples
#' bisec(-74.4824, 33, 98, 1, 0.252)
bisec <- function(delta,dea1,dea2,upbd,lowbd){
  a2init=0.5*(lowbd+upbd)
  lowvalue=dea1*log(lowbd)+dea2*log(1-lowbd)-delta
  a2initvalue=dea1*log(a2init)+dea2*log(1-a2init)-delta
  a=a2initvalue
  l=lowvalue
  if (l<a)
  {
    #monotone increasing
    while(abs(upbd-lowbd)> 0.00000001)
    {
      n=a2initvalue
      if (0<n) {
          upbd = a2init
          a2init = 0.5*(lowbd+upbd)
          a2initvalue = log(a2init)*dea1+log(1-a2init)*dea2-delta
          #lowbd
          #upbd
      } else if (0>n) {
          lowbd = a2init
          lowvalue = a2initvalue
          a2init = 0.5*(lowbd+upbd)
          a2initvalue = log(a2init)*dea1+log(1-a2init)*dea2-delta
          #lowbd
          #upbd
      } else {
          lowbd = upbd
          a2 = a2init
      }
    }
    a2=0.5*(lowbd+upbd)
  }
  if (l>a)
  {
     #monotone decreasing
     while (abs(upbd-lowbd)>0.00000001)
     {
       n=a2initvalue
      if (0>n)
      {
          upbd = a2init
          a2init = 0.5*(lowbd+upbd)
          a2initvalue = log(a2init)*dea1+log(1-a2init)*dea2-delta
           #lowbd
           #upbd
      } else if (0<n){
          lowbd = a2init
          lowvalue = a2initvalue
          a2init = 0.5*(lowbd+upbd)
          a2initvalue = log(a2init)*dea1+log(1-a2init)*dea2-delta
           #lowbd
           #upbd
      } else {
          lowbd = upbd
          a2 = a2init
      }
    }
    a2=0.5*lowbd+0.5*upbd
  }
  if (l==a)
    {
  #no change, means lowvalue == a2initvalue == 0;
  #                 a2 = 1 if upbd = 1;
  #                 a2 = 0 if upbd < 1.
      if (upbd == 1)
        a2 = upbd
      else
        a2 = lowbd
    }
    return(a2)
}