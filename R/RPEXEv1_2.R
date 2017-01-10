############################################################
#Input:
# 'EventTime' == A sequence of times where the events occur
# 'Censor'    == A sequence of dichotomous values indicating
#                censored or not (0=censored and 1=not censored)
# 'CutTimes'  == A vector of unique, sorted, possible times to
#                make the cuts. Default is sorted (from small to
#                large) event times
#                Default == 'EventTime'
# 'Trend'     == An input having indicating the monotonicity assumption
#                -- 0: no monotonic assumption
#                -- 1: failure rate is decreasing over time
#                -- 2: failure rate is increasing over time
#                -- 3: monotonic failure rate
#                -- 4: failure rate is increasing and then decreasing
#                -- 5: failure rate is decreasing and then increasing
#                -- 6: failure rate is increasing and then decreasing with
#                      the peak removed first
#                -- 7: failure rate is decreasing and then increasing with
#                      the peak removed first
#                Default == 0
# 'Criticalp'  == The critical (naive) p-value cutoff where all p-values 
#                in the backward elimination that are lower than this 
#                will be regarded as being significant. The prediction of 
#                the survival probability will be made on 100 equally
#                spaced time points within the range of the event times 
#                based on the piecewise exponential estimate determined by 
#                all the changepoints. 
#                Default == -1 (equivalent to NA).
# Output:
# pexeout.times   ==  times to make the cuts
# pexeout.pvalues ==  pvalues correspond to the times
# pexeout.times_c ==  critical times to make the cuts
# pexeout.pvalues_c ==  critical p-values that are smaller than the 
# pexeout.trend   ==  trend information
# pexeout.struct  ==  structure information for multiple order restrictions
# pexeout.changet ==  change point in time for umbrella alternatives.
#############################################################################

#' RPEXE
#'
#' @description  RPEXE
#'
#' @param eventtime A sequence of times where the events occur
#' @param censor A sequence of dichotomous values indicating censored or not (0=censored and 1=not censored)
#' @param cuttime A vector of unique, sorted, possible times to make the cuts. Default is sorted (from small to
#'                large) event times
#'                Default = 'EventTime'
#' @param trend An input having indicating the monotonicity assumption
#'                -- 0: no monotonic assumption
#'                
#'                -- 1: failure rate is decreasing over time
#'                
#'                -- 2: failure rate is increasing over time
#'                
#'                -- 3: monotonic failure rate
#'                
#'                -- 4: failure rate is increasing and then decreasing
#'                
#'                -- 5: failure rate is decreasing and then increasing
#'                
#'                -- 6: failure rate is increasing and then decreasing with
#'                      the peak removed first
#'                      
#'                -- 7: failure rate is decreasing and then increasing with
#'                      the peak removed first
#'                      
#'                Default == 0
#' @param criticalps  The critical (naive) p-value cutoff where all p-values 
#'                in the backward elimination that are lower than this 
#'                will be regarded as being significant. The prediction of 
#'                the survival probability will be made on 100 equally
#'                spaced time points within the range of the event times 
#'                based on the piecewise exponential estimate determined by 
#'                all the changepoints. 
#'                Default == -1 (equivalent to NA).
#'
#' @usage RPEXEv1_2(eventtime,censor,cuttime, trend,criticalps)
#'
#' @return
#' times: times to make the cuts
#' pvalues: pvalues correspond to the times
#' times_c: critical times to make the cuts
#' pvalues_c: critical p-values that are smaller than the 
#' trend: trend information
#' struct:  structure information for multiple order restrictions
#' changet:  change point in time for umbrella alternatives.
#' 
#' 
#' @export
#'
#' @examples
RPEXEv1_2 <- function(eventtime=NULL,censor=NULL,cuttime=NULL,trend=NULL, criticalps=NULL)
{


  cuttimes_default = 0
  monotone_default = 0
  criticalp_default = 0
  pexeout=list(times=as.null(),pvalues=as.null(), times_c = as.null(), pvalues_c = as.null(), trend=as.null(),struct=as.null(),changet=as.null(),
               plotdatakme_times=as.null(), plotdatakme_censoring=as.null(), plotdatapexe_t100=as.null(),
               plotdatapexe_pred100=as.null(), plotdatapexe_tchange=as.null(), plotdatapexe_predc=as.null())
  

  #reset parameter inputs
    if (!is.null(eventtime))
       times = eventtime
    if (!is.null(censor))
       censoring= censor
    if (!is.null(cuttime))
      {
       cuttimes    = cuttime
       cuttimes_default = 1
       }
    if (!is.null(trend))
      {
       monotone= trend
       monotone_default   = 1
      }
    if (!is.null(criticalps))
     {
       criticalp=criticalps
       criticalp_default = 1
     }
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
     cuttimes   = time_die[1:(length(time_die)-1)]
  #print(cuttimes)
  if (monotone_default == 0)
     monotone = 0
  if (criticalp_default == 0)
     criticalp = -1
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
  #     [time2,struct,label,indx] = umbrella(time_die,ttot,deaths,monotone-3);
  #     cuttimes_trend  = time2(1:(length(time2)-1));
  #     cuttimes        = intersect(cuttimes_trend, cuttimes);
  #     changetime      = time_die(indx);
  #     indsmall        = find(cuttimes< changetime);
  #     indlarge        = find(cuttimes>=changetime);
  #     [ts1,pvalues1]  = loopcuts(times,censoring,cuttimes(indsmall),1);
  #                         % decreasing Failure Rate
  #     [ts2,pvalues2]  = loopcuts(times,censoring,cuttimes(indlarge),2);
  #                         % increasing Failure Rate
  #     ts              = [ts1;ts2];
  #     pvalues         = [pvalues1;pvalues2];

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

  #     [time2,struct,label,indx] = umbrella(time_die,ttot,deaths,monotone-3);
  #     cuttimes_trend  = time2(1:(length(time2)-1));
  #     cuttimes        = intersect(cuttimes_trend, cuttimes);
  #     changetime      = time_die(indx);
  #     indsmall        = find(cuttimes< changetime);
  #     indlarge        = find(cuttimes>changetime);
  #     cuttimes        = [cuttimes(indsmall);changetime;cuttimes(indlarge)];
  #     mono            = [2*ones(length(indsmall),1);0;ones(length(indlarge),1)];
  #     [ts,pvalues]    = loopcuts_umbrella(times,censoring,cuttimes,mono);
  #     pexeout.struct  = struct;
  #     pexeout.changet = changetime;
  #     pexeout.trend   = 'Increasing-decreasing failure rate';

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
  
  if (criticalp_default == 1)
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
    legend("bottomleft",legend=c("Kaplan-Meier estimate","Exponential estimate"),
           col=c("black", "red"), lty=c(1,2),cex=0.8)
    
    # Save the data for makng the
    # not done
  }
  
  #times = pexeout$times 
  #values = pexeout$pvalues 
  #times_c =  pexeout$times_c 
  #pvalues_c =  pexeout$pvalues_c
  #trend   =  pexeout$trend 
  #struct  =  pexeout$structure 
  #changet = pexeout$changet
  
  # initialize the output vector
  #plotdatakme_times = vector()
  #plotdatakme_censoring = vector()
  #plotdatapexe_t100 = vector()
  #plotdatapexe_pred100 = vector()
  #plotdatapexe_tchange = vector()
  #plotdatapexe_predc = vector()
  pexeout$plotdatakme_times = times
  pexeout$plotdatakme_censoring = censoring
  pexeout$plotdatapexe_t100 = t100
  pexeout$plotdatapexe_pred100 = pred100
  if (min(pexeout$pvalues) < criticalp)
  {
    pexeout$plotdatapexe_tchange = tchange
    pexeout$plotdatapexe_predc = predc
  }
  #pexeout = list("times" = times, "pvalues"=pvalues, "times_c" =times_c, "pvalues_c" = pvalues_c, "trend" = trend, "struct"=struct, "changet" = changet,  "plotdatakme_times" = plotdatakme_times, "plotdatakme_censoring" = plotdatakme_censoring, "plotdatapexe_t100" = plotdatapexe_t100,
  #                 "plotdatapexe_pred100" = plotdatapexe_pred100, "plotdatapexe_tchange" = plotdatapexe_tchange, "plotdatapexe_predc" = plotdatapexe_predc)
  return(pexeout)
}