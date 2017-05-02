#' @title Log likelihood from the gamma distribution
#' 
#' @description A function computing the log likelihood from the gamma distribution under an order restriction reduction
#'
#' @param structtime change-point times to be used to compute the likelihood value
#' @param structttot total time on test (ttot) between each time point and the previous time point (or 0) corresponding to structtime
#' @param structdeaths number of deaths corresponding to structttot
#' @param time_die all event and censoring times from small to large
#' @param ttot total time on test corresponding to time_die
#' @param deaths the number of deaths corresponding to "ttot"
#' 
#' @usage gamllik(structtime,structttot,structdeaths,time_die,ttot,deaths)
#' 
#' @return
#' log of the likelihood
#' 
#' @export
#'
#' @examples
#' data(gamllik_data)
#' structtime = c(0.08493151, 1.89315068, 2.38630137, 12.70958904, 18.59452055, 23.46027397, 24.70958904, 28.03013699, 47.85479452, 51.07671233)
#' structttot = c(15.08493, 302.03836,  73.34795, 1090.68493, 345.04110, 199.20000, 41.22740, 99.71507, 369.89589,  82.45479)
#' structdeaths = c(2, 25, 6, 76, 22, 10, 2, 4, 6, 1)
#' time_die = gamllik_data[,1]
#' ttot = gamllik_data[,2]
#' deaths = gamllik_data[,3]
#' gamllik(structtime,structttot,structdeaths,time_die,ttot,deaths)
gamllik=function(structtime,structttot,structdeaths,time_die,ttot,deaths)
{
  #compute the Gamma parameter
  structgamindi=array(0,c(length(structtime),1))
  for(j in 1:length(structtime))
      structgamindi[j]= structttot[j]/structdeaths[j]
  #the likelihood
  #get the indices of the cut times:
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