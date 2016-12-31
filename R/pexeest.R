# a function produces the the quantile estimate at tx, 
# when a piecewise exponential distribution is fitted to 
# (times,cens) cens = 0 for censored, cens = 1 for uncensored.
# the change point is tchange
# lamest is the estimated parameters

#' pexeest
#' 
#' a function produces the the quantile estimate at tx, 
#' when a piecewise exponential distribution is fitted to 
#' (times,cens) cens = 0 for censored, cens = 1 for uncensored.
#'
#' @param times 
#' @param cens 
#' @param tchange 
#' @param tx 
#' 
#' @usage pexeest(times, cens, tchange, tx)
#'
#' @return
#' quan
#' lamest
#' 
#' @export
#'
#' @examples
pexeest <- function(times, cens, tchange, tx){
  
  quan = vector()
  tchange = sort(tchange)
  nchange = length(tchange)
  returnv=totaltest(times,cens) #returnv=list(time_die,ttot,deaths),the variables owns the same length
  m=dim(returnv)[2]/3
  time_die=returnv[,1:m]
  ttot=returnv[,(m+1):(2*m)]
  deaths=returnv[,(2*m+1):3*m]
  ntime = length(time_die)
  
  if (nchange >= 1)
  {
    # find the index
    indchange = rep(0, nchange)
    for (j in 1:nchange)
    {
      for (i in 1:ntime)
      {
        if (abs(time_die[i] - tchange[j]) < 0.00001) # for round off error
        { indchange[j] = i
        }
      }
    }
    if (length(unique(indchange)) < length(indchange))
    {
      # there is a need to look for unique variables
      a = order(indchange)[!duplicated(sort(indchange))]
      tchange = tchange[a]
      nchange = length(tchange)
    }
    
    # estimate the piecewise exponential parameter lambda1-lambda_nchange
    # compute all the changepoint quantile
    lamest = rep(0, nchange+1)
    E = rep(0, nchange)
    for (j in 1:nchange)
    {
      if (j == 1)
      {
        lamest[j] = sum(ttot[1:indchange[j]])/sum(deaths[1:indchange[j]])
        E[j] = exp(-tchange[j]/lamest[j])
      } else {
        a = ttot[(indchange[j-1]+1):indchange[j]]
        b = deaths[(indchange[j-1]+1):indchange[j]]
        lamest[j] = sum(a)/sum(b)
        E[j] = E[j-1]*exp(-(tchange[j]-tchange[j-1])/lamest[j])
      }
    }
    lamest[nchange+1] = sum(ttot[(indchange[nchange]+1):ntime])/sum(deaths[(indchange[nchange]+1):ntime])
    for (k in 1:length(tx))
    {
      if (tx[k] < tchange[1])
      {
        quan[k] = exp(-tx[k]/lamest[1])
      } else if (tx[k] < tchange[nchange]){
        for (j in 2:nchange)
        {
          if (tx[k] >= tchange[j-1])
          {
            if (tx[k] < tchange[j])
            {                
              quan[k] = E[j-1]*exp(-(tx[k]-tchange[j-1])/lamest[j])
            }
          }
        }
      } else{
          #piece nchange+1
          c = -(tx[k]-tchange[nchange])/lamest[(nchange+1)]
          quan[k] = E[nchange]* exp(c)
      }
    }
    
    
  }else{
    # nchange < 1
    lamest = sum(ttot)/sum(deaths)
    for (k in 1:length(tx))
    {
      quan[k] = exp(-tx[k]/lamest)
    }
  }
  # 
  # tchange
  # E
  # lamest
  pexeout=list("quan"=quan, "lamest"=lamest)
  return (pexeout)
  
}