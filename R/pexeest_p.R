# a function produces the the quantile estimate at tx, 
# when a piecewise exponential distribution is fitted to 
# (times,cens) cens = 0 for censored, cens = 1 for uncensored.
# the change point is tchange
# lamest is the estimated parameters

#' pexeest_p
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
#' @usage pexeest_p(times, cens, tchange, tx)
#'
#' @return
#' pchange
#' 
#' @export
#'
#' @examples
pexeest_p <- function(times, cens, tchange, tx){
  
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
        ttot1 = sum(ttot[1:indchange[1]])
        d1 = sum(deaths[1:indchange[1]])
        ttot2 = sum(ttot[(indchange[1]+1):indchange[2]])
        d2 = sum(deaths[(indchange[1]+1):indchange[2]])
        reout = exact_pvalue(ttot1,ttot2,d1,d2,mono)
        pchange[1] = reout[2]
      } else if (j==nchange) {
        ttot1 = sum(ttot[(indchange[j-1]+1):indchange[j]])
        d1 = sum(deaths[(indchange[j-1]+1):indchange[j]])
        ttot2 = sum(ttot[(indchange[j]+1):length(ttot)])
        d2 = sum(deaths[(indchange[j]+1):length(ttot)])
        reout = exact_pvalue(ttot1,ttot2,d1,d2,mono)
        pchange[j] = reout[2]
      } else{
        ttot1 = sum(ttot[(indchange[j-1]+1):indchange[j]])
        d1 = sum(deaths[(indchange[j-1]+1):indchange[j]])
        ttot2 = sum(ttot[(indchange[j]+1):indchange[j+1]])
        d2 = sum(deaths[(indchange[j]+1):indchange[j+1]])
        reout = exact_pvalue(ttot1,ttot2,d1,d2,mono)
        pchange[j] = reout[2]
      }
    }

  }
  # 
  # tchange
  # E
  # lamest
  return (pchange)
  
}