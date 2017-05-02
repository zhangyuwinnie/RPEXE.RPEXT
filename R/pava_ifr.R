#' @title PAVA order restriction under decreasing failure rate (DFR)
#' 
#' @description This function imposes the PAVA DFR order restriction by eliminating change-points violating the restriction 
#' 
#' @usage pava_ifr(time_die,ttot,deaths)
#' 
#' @param time_die event times
#' @param ttot the total time on test (ttot) corresponding to the event times
#' @param deaths the number of deaths at each event time
#' 
#' @return
#' time2 the event times after PAVA
#' ttot2 the corresponding ttot after PAVA
#' deaths2 the corresponding number of deaths after PAVA
#' 
#' @export
#'
#' @examples
#' data(pava_dfrd)
#' t_d = pava_dfrd[,1]
#' t = pava_dfrd[,2]
#' d = pava_dfrd[,3]
#' pava_ifr(t_d, t, d)
#'
pava_ifr <- function(time_die,ttot,deaths)
{
    ttotrev= (-1)*ttot
    returnval=pava_dfr(time_die,ttotrev,deaths)
    m=dim(returnval)[2]/3
    time2=returnval[,1:m]
    ttot3=returnval[,(m+1):2*m]
    deaths2=returnval[,(2*m+1):3*m]
    #[time2, ttot3, deaths2] = retuenval
    ttot2=-ttot3
    returnval_t=cbind(time2,ttot2,deaths2)
    return(returnval_t)
}
