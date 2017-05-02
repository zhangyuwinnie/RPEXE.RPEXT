#' @title RPEXE main function
#' 
#' @description This is the RPEXE main function taking inputs including time, censoring, 
#' change-point candidates, order restriction, criticl value, and display position. This function 
#' produces the RPEXE estimate. The prediction of the survival probability will be made on 100 equally
#' spaced time points within the range of the event times based on the piecewise exponential 
#' estimate determined by all the changepoints. 
#'
#' @param times A sequence of times where the events occur
#' @param censoring A sequence of dichotomous values indicating censored or not (0=censored and 1=not censored)
#' @param cuttimes A vector of unique, sorted, possible times to make the cuts. When it's set to NULL, it's the Default value, 
#'                 which is sorted event times from small to large. 
#' @param monotone An input having indicating the monotonicity assumption
#'                -- 0: no monotonic assumption (default)
#'                -- 1: failure rate is decreasing over time
#'                -- 2: failure rate is increasing over time
#'                -- 3: monotonic failure rate
#'                -- 4: failure rate is increasing and then decreasing
#'                -- 5: failure rate is decreasing and then increasing
#'                -- 6: failure rate is increasing and then decreasing with
#'                      the peak removed first
#'                -- 7: failure rate is decreasing and then increasing with
#'                      the peak removed first
#'                                      
#' @param criticalp  The critical (naive) p-value cutoff where all p-values 
#'                in the backward elimination that are lower than this 
#'                will be regarded as being significant. For example, at type I error rate 0.05, 
#'                the critical p-value was 0.004 in the real example of Han et al. (2014).
#'                Default == -1 (equivalent to NA).
#' @param pos The position of the legend. Can be 0 or 1. The legend will be 
#'    on the topright if set to 0. The legend will be on the bottomleft if set to 1. Default is 0.
#'
#' @usage RPEXEv1_2(times,censoring,cuttimes=NULL, monotone=0, criticalp=-1, pos = 0)
#'
#' @return
#' times: event/censoring times taking out from the backward elimination
#' pvalues: p-values corresponding to "times" 
#' times_c: significant change-points
#' pvalues_c: critical p-values that are smaller than the critical p-value 
#' trend: trend information
#' struct:  structure information for multiple order restrictions
#' changet:  change-point time of trend for umbrella alternatives.
#' 
#' 
#' @export
#'
#' @examples
#' t1 <- c(2,3,4,5.5,7,10,12,15)
#' c1 <- c(0,0,1,0,0,1,0,0)
#' RPEXEv1_2(t1, c1, monotone = 1,criticalp=0.05, pos = 0)
RPEXEv1_2 <- function(times, censoring, cuttimes=NULL, monotone = 0, criticalp=-1, pos = 0)
{
  # initialize the output list
  pexeout=list(times=as.null(),pvalues=as.null(), times_c = as.null(), pvalues_c = as.null(), trend=as.null(),struct=as.null(),changet=as.null(),
               plotdatakme_times=as.null(), plotdatakme_censoring=as.null(), plotdatapexe_t100=as.null(),
               plotdatapexe_pred100=as.null(), plotdatapexe_tchange=as.null(), plotdatapexe_predc=as.null())
  
  #reset parameter inputs
    if (!is.null(cuttimes))
      cuttimes_default = 1
    else
      cuttimes_default = 0
  
  
  # Compute the time of the death(in increasing order), the total time on
  # test, and the number of deaths corresponding to ttot.
  # These quantities will be used later.
  returnv=totaltest(times,censoring) #returnv=list(time_die,ttot,deaths),the variables owns the same length
  m=dim(returnv)[2]/3
  time_die=returnv[,1:m]
  ttot=returnv[,(m+1):(2*m)]
  deaths=returnv[,(2*m+1):3*m]
 
  # Set the default values
  if (cuttimes_default == 0)
     cuttimes = time_die[1:(length(time_die)-1)]

  if (monotone== 0)
     {
     returnva_l=loopcuts(times,censoring,cuttimes,monotone)
     n=dim(returnva_l)[2]/2
     ts=returnva_l[,1:n]
     pvalues=returnva_l[,(n+1):2*n]
     pexeout$trend="No order restriction"
     }
  if (monotone == 1)
     {
     returnva_p=pava_dfr(time_die,ttot,deaths)
     n=dim(returnva_p)[2]/3
     time2=returnva_p[,1:n]
     ttot2=returnva_p[,(n+1):(2*n)]
     deaths2=returnva_p[,(2*n+1):3*n]
     cuttimes_trend= time2[1:(length(time2)-1)]
     cuttimes= sort(intersect(cuttimes_trend, cuttimes))
     returnva_l=loopcuts(times,censoring,cuttimes,monotone)
     n=dim(returnva_l)[2]/2
     ts=returnva_l[,1:n]
     pvalues=returnva_l[,(n+1):2*n]
     pexeout$trend="Decreasing falilure rate"
     }
  if (monotone == 2)
     {
     returnva_p= pava_ifr(time_die,ttot,deaths)
     n=dim(returnva_p)[2]/3
     time2=returnva_p[,1:n]
     ttot2=returnva_p[,(n+1):(2*n)]
     deaths2=returnva_p[,(2*n+1):3*n]
     cuttimes_trend=time2[1:(length(time2)-1)]
     cuttimes= sort(intersect(cuttimes_trend, cuttimes))
     returnva_l=loopcuts(times,censoring,cuttimes,monotone)
     n=dim(returnva_l)/2
     ts=returnva_l[,1:n]
     pvalues=returnva_l[,(n+1):2*n]
     pexeout$trend="Increasing falilure rate"
     }
  if (monotone == 3)
     {
      returnva_u=umbrella(time_die,ttot,deaths,monotone-3)
      time2=returnva_u[1]
      struct=returnva_u[2]
      label=returnva_u[3]
      indx=returnva_u[4]
      time2=time2[[1]][]
      cuttimes_trend  = time2[1:(length(time2)-1)]
      cuttimes =sort(intersect(cuttimes_trend, cuttimes))
      if (struct[[1]][[4]][1]>struct[[1]][[4]][2])
         monotone = 1
      else
         monotone = 2
      returnva_l=loopcuts(times,censoring,cuttimes,monotone)
      n=dim(returnva_l)[2]/2
      ts=returnva_l[,1:n]
      pvalues=returnva_l[,(n+1):(2*n)]
      pexeout$struct= struct
      pexeout$trend="Monotone failure rate"
      }
  if (monotone==4) #increasing then decreasing failure rate
      {
      returnva_u=umbrella(time_die,ttot,deaths,monotone-3)
      time2=returnva_u[1]
      struct=returnva_u[2]
      label=returnva_u[3]
      indx=returnva_u[4]
      #indx
      #label
      time2=time2[[1]][]
      cuttimes_trend  = time2[1:(length(time2)-1)]
      cuttimes        = sort(intersect(cuttimes_trend,cuttimes))
      indx=indx[[1]][1]
      changetime      = time_die[indx]
      cuttimes        = c(cuttimes[cuttimes< changetime],changetime,cuttimes[cuttimes>changetime])
      mono            = c(2*matrix(1,nrow=length(cuttimes[cuttimes<changetime]),ncol=1),0,matrix(1,nrow=length(cuttimes[cuttimes>changetime]),ncol=1))
      returnva_l=loopcuts_umbrella(times,censoring,cuttimes,mono)
      n=dim(returnva_l)[2]/2
      ts=returnva_l[,1:n]
      pvalues=returnva_l[,(n+1):(2*n)]
      pexeout$struct= struct
      pexeout$changet = changetime
      pexeout$trend   = "Increasing-decreasing failure rate"
  }
  
  if  (monotone==5) #decreasing then increasing failure rate
     {
      returnva_u=umbrella(time_die,ttot,deaths,monotone-3)
      time2=returnva_u[1]
      struct=returnva_u[2]
      label=returnva_u[3]
      indx=returnva_u[4]
      time2=time2[[1]][]
      cuttimes_trend  = time2[1:(length(time2)-1)]
      cuttimes =sort(intersect(cuttimes_trend, cuttimes))
      indx=indx[[1]][1]
      changetime      = time_die[indx]
      cuttimes        = c(cuttimes[cuttimes< changetime],changetime,cuttimes[cuttimes>changetime])
      mono            = c(matrix(1,nrow=length(cuttimes[cuttimes< changetime]),ncol=1),0,2*matrix(1,nrow=length(cuttimes[cuttimes>changetime]),ncol=1))
      returnva_lu=loopcuts_umbrella(times,censoring,cuttimes,mono)
      n=dim(returnva_lu)[2]/2
      ts=returnva_lu[,1:n]
      pvalues=returnva_lu[,(n+1):2*n]
      pexeout$struct  = struct
      pexeout$changet = changetime
      pexeout$trend   = "Decreasing-increasing failure rate"
      }
  if (monotone==6) # increasing then decreasing failure rate, no peak
     {
      returnva_u=umbrella(time_die,ttot,deaths,monotone-3)
      time2=returnva_u[1]
      struct=returnva_u[2] #struct=data.frame(time,ttot,deaths,loglik)
      label=returnva_u[3]
      indx=returnva_u[4]
      time2=time2[[1]][]
      cuttimes_trend  = time2[1:(length(time2)-1)]
      cuttimes =sort(intersect(cuttimes_trend, cuttimes))
      indx=indx[[1]][1]
      changetime      = time_die[indx]
      cuttimes        = c(cuttimes[cuttimes< changetime],cuttimes[cuttimes>changetime])
      mono            = c(2*matrix(1,nrow=length(cuttimes[cuttimes< changetime]),ncol=1),matrix(1,nrow=length(cuttimes[cuttimes>changetime]),ncol=1))
      returnva_lu=loopcuts_umbrella(times,censoring,cuttimes,mono)
      n=dim(returnva_lu)[2]/2
      ts=returnva_lu[,1:n]
      pvalues=returnva_lu[,(n+1):2*n]
      pexeout$struct  = struct
      pexeout$changet = changetime
      pexeout$trend   ="Increasing-decreasing failure rate"
     }

  if (monotone == 7) #decreasing then increasing failure rate, no peak
      {
      returnva_u=umbrella(time_die,ttot,deaths,monotone-3)
      time2=returnva_u[1]
      struct=returnva_u[2]
      label=returnva_u[3]
      indx=returnva_u[4]
      time2=time2[[1]][]
      cuttimes_trend  = time2[1:(length(time2)-1)]
      cuttimes        = sort(intersect(cuttimes_trend,cuttimes))
  #   cuttimes_trend
  #   cuttimes
  #   changetime
      indx=indx[[1]][1]
      changetime = time_die[indx]
  #   cuttimes
      cuttimes        = c(cuttimes[cuttimes<changetime],cuttimes[cuttimes>changetime])
  #   mono
      mono=c(matrix(1,nrow=length(cuttimes[cuttimes<changetime]),ncol=1),2*matrix(1,nrow=length(cuttimes[cuttimes>changetime]),ncol=1))
      returnva_lu=loopcuts_umbrella(times,censoring,cuttimes,mono)
      n=dim(returnva_lu)[2]/2
      ts=returnva_lu[,1:n]
      pvalues=returnva_lu[,(n+1):2*n]
      pexeout$struct  = struct
      pexeout$changet = changetime
      pexeout$trend   ="Decreasing-increasing failure rate"
      }
  pexeout$times       = unique(ts[pvalues<1])
  #print("pexeout$times")
  #print(pexeout$times)
  pexeout$pvalues     = unique(pvalues[pvalues<1])
  
  # The following part is new in the RPEXE version 2.
  # 1. In the new version the program will determine location of the 
  # change point given a naive p-value.
  # 2. The program is able to plot the estimated survival function 
  # using the change point determined by the naive p-value and overlay 
  # the estimated piecewise exponential estimate with the Kaplan-Meier curve.
  
  if (criticalp != -1)
  {
    # compute the grid of 100 times that is equally spaced between 0 and the
    # maximum of the event time; 
    tmax = max(times)
    t100 = seq(from = 0.01, to = 1, by = 0.01)*tmax
    
    # Compute the predicted value using the grid and use it as 
    # critical sentence: when any p-value is less than the critical 
    #     value, keep all the times at that point.
    
    # if the change point is found
    if (min(pexeout$pvalues) < criticalp)
    {
      # exclude 1:min()-1 from times
      last = min(which(pexeout$pvalues<criticalp))
      print("last")
      print(last)
      if(last == 1)
      {
        tchange = pexeout$times
      } else{
        tchange = pexeout$times[-(1:last-1)]
      }
      pexeest1 = pexeest(times, censoring, tchange, t100)
      pred100 = pexeest1$quan
      lamest = pexeest1$lamest
      # save critical time and p-values;
      tchange = sort(tchange)
      # pexeout.pvalues_c = ...
      # pexeout.pvalues(min(find(pexeout.pvalues < criticalp)):end);
      # compute the critical value at the changetime of tchange;
      pexeest2 = pexeest(times, censoring, tchange, tchange)
      predc = pexeest2$quan
      lamest = pexeest2$lamest
      # pchange = pexeest_p(times, censoring, tchange, mono);  
      # program mono
      if (monotone < 4)
      {
        mono = monotone * rep(1,length(tchange))
      } else{
        indsmall = which(tchange < changetime)
        indlarge = which(tchange > changetime)
        if (monotone == 4)
        {
          if (is.element(changetime, tchange))
          { 
            #order restriction changepoint in the critical times set
            tchange = c(tchange[indsmall], changetime, tchange[indlarge])
            mono = c(2*rep(1,length(indsmall)),0,rep(1,length(indlarge)))
          } else{
            #order restriction changepoint not in the critical times set
            tchange = c(tchange[indsmall], tchange[indlarge])
            mono = c(2*rep(1,length(indsmall)),rep(1,length(indlarge)))
          }
        } else if (monotone == 5){
          if (is.element(changetime, tchange))
          { 
            #order restriction changepoint in the critical times set
            tchange = c(tchange[indsmall], changetime, tchange[indlarge])
            mono = c(rep(1,length(indsmall)),0,2*rep(1,length(indlarge)))
          } else{
            #order restriction changepoint not in the critical times set
            tchange = c(tchange[indsmall], tchange[indlarge])
            mono = c(rep(1,length(indsmall)),2*rep(1,length(indlarge)))
          }
        } else if (monotone == 6){
          tchange = c(tchange[indsmall], tchange[indlarge])
          mono = c(2*rep(1,length(indsmall)), rep(1,length(indlarge)))
        } else if (monotone == 7){
          tchange = c(tchange[indsmall], tchange[indlarge])
          mono = c(rep(1,length(indsmall)),2*rep(1,length(indlarge)))
        }
      }
      pexeout$times_c = tchange
      pexeout$lamest = lamest
      pvalall = loopcuts_onestep(times,censoring,tchange,mono)
      pexeout$pvalues_c = pvalall
      # If an output has p-value 
        
    }else {
      # Note that if we set the changepoint to be the maximum time in the
      # record, then the function will return the exponential estimate;
      # [pred100 lamest] = pexeest(times, censoring, max(times), t100);
      expar1 = sum(times)/sum(censoring)
      pred100 = exp(-t100/expar1)
      pexeout$lamest  = expar1 
    }
    
    # Plot the Kaplan-Meier curve overlaid with the estimates;
    # original scale;
    # figure(11)
    km(times, censoring,1)
    lines(t100,pred100, type = "l", lty = 2, lwd = 2, col = "red")
    if (min(pexeout$pvalues)<criticalp)
    {
      
      points(tchange,predc,pch = 24,lwd = 4,col="red")
    }
    #legend("bottomleft",legend=c("Kaplan-Meier estimate","Exponential estimate","Significant failure rate change"),
    #       col=c("black", "red","red"), lty=c(1,2,NA),pch = c(NA,NA,24), cex=0.8)
    
    # figure(12)
    km_log(times, censoring,1)
    lines(t100,log(pred100), type = "l", lty = 2, lwd = 2, col = "red")
    if (min(pexeout$pvalues)<criticalp)
    {
     
      points(tchange,log(predc),pch = 24,lwd = 4, col="red")
    }
    #legend("topright",legend=c("Kaplan-Meier estimate","Exponential estimate","Significant failure rate change"),
    #       col=c("black", "red","red"), lty=c(1,2,NA),pch = c(NA,NA,24), cex=0.8)
    
    # figure(13)
    km(times, censoring,0)
    lines(t100,pred100, type = "l", lty = 2, lwd = 2, col = "red")
    if (min(pexeout$pvalues)<criticalp)
    {
     
      points(tchange,predc,pch = 24,lwd = 4,col="red")
    }
    if (pos == 0)
      position = "topright"
    if (pos == 1)
      position = "bottomleft"
    legend(position,legend=c("Kaplan-Meier estimate","Exponential estimate"),
           col=c("black", "red"), lty=c(1,2),cex=0.8)
  }
  pexeout$plotdatakme_times = times
  pexeout$plotdatakme_censoring = censoring
  pexeout$plotdatapexe_t100 = t100
  pexeout$plotdatapexe_pred100 = pred100
  if (min(pexeout$pvalues) < criticalp)
  {
    pexeout$plotdatapexe_tchange = tchange
    pexeout$plotdatapexe_predc = predc
  }
  return(pexeout)
}