#' Plot a Kaplan Meier curve in log scale
#'
#' @description The function plots a Kaplan Meier curve in log scale
#' @param time time of observed event
#' @param censor a vector indicating censored or not at the given times, 0 = censored; 1 = uncensored
#' @param plotcens 0: add censored data to the output curve
#' 
#'                 1: don't add censored data to the output curve
#'
#' @usage km_log(time, censor, plotcens)
#' @return
#' A Kaplan Meier curve in log scale
#' 
#' @export
#'
#' @examples
#' t1 <- c(2,3,4,5.5,7,10,12,15)
#' c1 <- c(0,0,1,0,0,1,0,0)
#' km_log(t1,c1,0)
km_log <- function(time, censor, plotcens){
  # compute realt and deaths
  tmptime = time
  #cat(tmptime)
  for (i in 1:length(censor))
  {
    if (censor[i] == 0)
    {
      tmptime[i] = 0
    }
  }
  timesort = sort(tmptime)
  realt = unique(timesort)
  # input #
  # realt:  sorted times when death occur,
  # deaths: corresponding number of deaths at the sorted times
  # time:   all the times, including censor and death,
  # censor: a vector indicating censored or not at the given times in time
  #         censored: censor() = 0; noncensored: censor() =1;
  
  # output #
  # the label on the kaplan-mejer curve, which is corresponding to each
  # point on realt and the plot of the kaplan-meier curve.
  
  ldea = length(realt)
  pos_km = rep(0, ldea)
  # at each time point, the key is to compute ni and di,
  # the kaplan meier estimate has the form
  #       est_s(t) = prod_{ti<t}((ni-di)/ni);
  #
  # get the values of the first item in pos_km;
  
  # compute two sequences ni and di
  niseq = rep(0, ldea)
  diseq = rep(0, ldea)
  difnidi = rep(0, ldea)
  for (i in 1:ldea)
  {
    numberp = rep(1,length(censor))
    for (ii in 1:length(censor))
    {
      if (time[ii] < realt[i])
      {
        numberp[ii] = 0
      }
      if (time[ii] <= realt[i] && censor[ii] == 1)
      {
        diseq[i] = diseq[i] + 1
      }
      niseq[i] = sum(numberp)
    }
  }
  
  #substract the cumulative death
  for (j in 2:ldea)
  {
    diseq[ldea+2-j] = diseq[ldea+2-j] - diseq[ldea+2-j-1]
  }
  
  difnidi = niseq - diseq
  elekm = rep(0, ldea)
  
  for (i in 1:ldea)
  {
    elekm[i] = difnidi[i]/niseq[i]
    pos_km[i] = prod(elekm[1:i])
  }
  
  # given pos_km, plot kaplan-meier curve
  upperx              = rep(0, ldea+1)
  lowerx              = rep(0, ldea+1)
  upperx[1]           = 0
  upperx[2:(ldea+1)]  = realt
  lowerx[1:ldea]      = realt
  lowerx[1+ldea]      = max(time)
  uppery              = rep(0, ldea+1)
  lowery              = rep(0, ldea+1)
  uppery[1:2]         = 1
  uppery[2:ldea+1]    = pos_km[1:ldea-1]
  
  lowery[1:ldea]      = pos_km
  lowery[1+ldea]      = pos_km[ldea]
  
  xpart               = rep(0, 2*ldea+2)
  ypart               = rep(0, 2*ldea+2)
  xpartl              = lowerx[1:ldea]
  xpartu              = upperx[2:ldea+1]
  ypartl              = lowery[1:ldea]
  ypartu              = uppery[2:ldea+1]
  ypart
  # save the data in xpart and ypart
  xpart[1] = upperx[1]
  xpart[2] = upperx[2]
  for (i in 2:ldea)
  {
    xpart[2*(i-1)+1] = lowerx[i-1]
    xpart[2*(i-1)+2] = upperx[i+1]
  }
  xpart[2*ldea+1] = lowerx[ldea]
  xpart[2*ldea+2] = lowerx[ldea+1]
  
  ypart[1] = uppery[1]
  ypart[2] = uppery[2]
  for (i in 2:ldea)
  {
    ypart[2*(i-1)+1]= lowery[i-1]
    ypart[2*(i-1)+2]= uppery[i+1]
  }
  ypart[2*ldea+1] = lowery[ldea]
  ypart[2*ldea+2] = lowery[ldea+1]
  
  #create the log of the y axis and plot it
  logypart = log(ypart)
  # plot the main survival curve
  plot(xpart, logypart, type = "s",xlab="Years", ylab="Log of survival probability")
  
  # add the points of censor data
  if (plotcens == 1){
    for (i in 1:length(censor))
    {
      tcen = 0
      if (censor[i] == 0)
      {
        tcen = time[i]
      }
      # find the value
      k = 1
      while (xpart[k] < tcen)
      {
        k = k+1
      }
      if (tcen != 0)
      {
        points(tcen, logypart[k],cex=.8,pch=1)
      }
    }
  }
  
  
}
