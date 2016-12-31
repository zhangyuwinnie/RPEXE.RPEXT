# A function computing the log likelihood from the gamma distribution under
#       an order restriction reduction
# Inputs
#       structtime   == times under restriction
#       structttot   == ttots between each time point and the previous
#                           time point (or 0) under restriction
#       structdeaths == number of deaths corresponding to struct.ttot
#       time_die      == all possible times to make the cuts
#       ttot          == ttots corresponding to time_die
#       deaths        == dealthes corresponding to ttot
# Outputs
#       loglik        == log of the likelihood

#' @title Computing the log likelihood from the gamma distribution
#' 
#' @description A function computing the log likelihood from the gamma distribution under
#'       an order restriction reduction
#'
#' @param structtime   times under restriction
#' @param structttot  ttots between each time point and the previous
#'                           time point (or 0) under restriction
#' @param structdeaths number of deaths corresponding to struct.ttot
#' @param time_die all possible times to make the cuts
#' @param ttot ttots corresponding to time_die
#' @param deaths dealthes corresponding to ttot
#' 
#' @usage gamllik(structtime,structttot,structdeaths,time_die,ttot,deaths)
#' 
#' @return
#' log of the likelihood
#' 
#' @export
#'
#' @examples
#' 
gamllik=function(structtime,structttot,structdeaths,time_die,ttot,deaths)
{

  #compute the Gamma parameter
  structgamindi=array(0,c(length(structtime),1))
  for(j in 1:length(structtime))
      structgamindi[j]= structttot[j]/structdeaths[j]
  #the likelihood
  #get the indices of the cut time:
  structindi=array(0,c(length(structtime),1))
  for(j in 1:length(structtime))
    for (jj in 1:length(time_die))
       if (time_die[jj]==structtime[j])
           structindi[j] =jj

  #set the scale parameter for Gamma distribution of the ttot
  structgampar=array(0,c(length(time_die),1))
  for (ii in 1:length(time_die))
      for (j in 1:length(structtime))
          if (ii<=structindi[j])
            {
              if (j==1)
                  structgampar[ii]= structgamindi[j]
              else
                  if (ii>structindi[j-1])
                      structgampar[ii] = structgamindi[j]
             }
  loglik = 0
  for (ii in 1:length(structgampar))
      loglik = loglik + log(dgamma(ttot[ii]/deaths[ii],shape=deaths[ii],scale=structgampar[ii],log=FALSE))
  return(loglik)
}