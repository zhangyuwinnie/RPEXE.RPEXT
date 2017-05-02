#' @title Umbrella alternative.
#' 
#' @description Using the umbrella alternative to merge certain entries to make 
#' the sequence of ttot/deaths to increase then decrease or to decrease then increase. 
#' Note that the pava function imposes non-decreasing or non-increasing order. 
#' This function directly uses function pava().
#'
#' @param time_die a sequence of times where deaths happened.
#' @param ttot the total time on test between each time point
#'                     and the previous time point (or 0).
#' @param deaths the number of deaths at each time point.
#' @param indi an indicator
#'                indi == 0: monotonic failure rate (either decrease or increase)
#'                indi == 1: denoting the failure rate increase then decrease
#'                indi == 2: denoting the failure rate decrease then increase
#'
#' @usage umbrella(time_die,ttot,deaths,indi)
#' 
#' @return
#' time2  == the merged time_die after the umbrealla alternative order restriction;
#' struct  == a structure saves the partition information;
#' label  == a note about how the failure rate varies;
#' indx   == the position where the change point value is.
#' 
#' @export
#'
#' @examples
#' data(pava_dfrd)
#' t_d = pava_dfrd[,1]
#' t = pava_dfrd[,2]
#' d = pava_dfrd[,3]
#' umbrella(t_d, t, d, 2)
umbrella <- function(time_die,ttot,deaths,indi)
{
  n=length(time_die)
  time_all = matrix(0,nrow=n+1,ncol=1)
  time_all[1]=0
  time_r=array(0,c(n,n))
  ttot_r=array(0,c(n,n))
  deaths_r=array(0,c(n,n))
  loglik_r=array(0,c(1,n))

  for (i in 2:(n+1))
     time_all[i]= time_die[i-1]
  if (indi==0) # monotonic failure rate.
    {
     # in the structure, the first element is decreasing rate; the second element is increasing rate
     returnva_p1=pava_dfr(time_die,ttot,deaths)#returnva_p=list(struct$time[1],struct$ttot[1],struct$deaths[1])
     m=dim(returnva_p1)[2]/3
     w1=length(returnva_p1)/3
     time_r[1:w1,1]=returnva_p1[,1:m]
     ttot_r[1:w1,1]=returnva_p1[,(m+1):2*m]
     deaths_r[1:w1,1]=returnva_p1[,(2*m+1):3*m]
     returnva_p2=pava_ifr(time_die,ttot,deaths)#returnva_p2=list(struct$time[2],struct$ttot[2],struct$deaths[2])
     m=dim(returnva_p2)[2]/3
     w2=length(returnva_p2)/3
     time_r[1:w2,2]=returnva_p2[,1:m]
     ttot_r[1:w2,2]=returnva_p2[,(m+1):2*m]
     deaths_r[1:w2,2]=returnva_p2[,(2*m+1):3*m]
     #computing the likelihood
     for (i in 1:2)
       {
         t=time_r[,i]
         t=t[t!=0]
         tt=ttot_r[,i]
         tt=tt[tt!=0]
         d=deaths_r[1:length(tt),i]
       loglik_r[1,i]=gamllik(t,tt,d,time_die,ttot,deaths)
       }
     #select the times using the larger loglikelihood
     loglik=loglik_r[1,1]
     time2=time_r[,1]
     i=2
     if (loglik<loglik_r[1,i])
       {
        time2  = time_r[,i]
        loglik = loglik_r[1,i]
        label  = "Decreasing failure rate."
        indx   = 1
        }
     else
        {
        label  ="Increasing failure rate."
        indx   = length(time2)
        }
    }
  if (indi==1)#failure rate increases then decreases. ttot/death U shape
    {
    for (j in 1:length(time_die))
      {
        if (j==1)
          {
          returnva_p3=pava_dfr(time_die,ttot,deaths)#returnva_p3=cbind(struct$time[j],struct$ttot[j],struct$deaths[j])
           m=dim(returnva_p3)[2]/3
           w1=length(returnva_p3)/3
           time_r[1:w1,j]=returnva_p3[,1:m]
           ttot_r[1:w1,j]=returnva_p3[,(m+1):2*m]
           deaths_r[1:w1,j]=returnva_p3[,(2*m+1):3*m]
          }
        if (j==length(time_die))
          {
           returnva_p4=pava_ifr(time_die,ttot,deaths)#returnva_p4=list(struct$time[j],struct$ttot[j],struct$deaths[j])
           m=dim(returnva_p4)[2]/3
           w2=length(returnva_p4)/3
           time_r[1:w2,j]=returnva_p4[,1:m]
           ttot_r[1:w2,j]=returnva_p4[,(m+1):2*m]
           deaths_r[1:w2,j]=returnva_p4[,(2*m+1):3*m]
           }
        if ((j!=1)&&(j<length(time_die)))
           {
           returnva_p5=pava_ifr(time_die[1:j],ttot[1:j],deaths[1:j])
           m=dim(returnva_p5)[2]/3
           w3=length(returnva_p5)/3
           tmp_time1=returnva_p5[,1:m]
           tmp_ttot1=returnva_p5[,(m+1):2*m]
           tmp_deaths1=returnva_p5[,(2*m+1):3*m]
           returnva_p6=pava_dfr(time_die[-(1:j)],ttot[-(1:j)],deaths[-(1:j)])
           m=dim(returnva_p6)[2]/3
           w4=length(returnva_p6)/3
           tmp_time2=returnva_p6[,1:m]
           tmp_ttot2=returnva_p6[,(m+1):2*m]
           tmp_deaths2=returnva_p6[,(2*m+1):3*m]
           time_r[1:(w3+w4),j] = c(tmp_time1,tmp_time2)
           ttot_r[1:(w3+w4),j]= c(tmp_ttot1,tmp_ttot2)
           deaths_r[1:(w3+w4),j] = c(tmp_deaths1,tmp_deaths2)
           }
         #compute the loglikelihood for struct(j)
         t=time_r[,j]
         t=t[t!=0]
         tt=ttot_r[,j]
         tt=tt[tt!=0]
         d=deaths_r[1:length(tt),j]
         loglik_r[1,j]= gamllik(t,tt,d,time_die,ttot,deaths)
       }
     label="increasing then decreasing failure rate"
     #select the times using the largest loglikelihood
     loglik=loglik_r[1,1]
     time2=time_r[,1]
     indx    = 1
     for (j in 1:length(time_die))
        if (loglik < loglik_r[1,j])
           {
            time2  = time_r[,j]
            loglik = loglik_r[1,j]
            indx   = j
            }
     }
  if (indi==2)#failure rate decreases then increases. ttot/death arc shape
    {
     for (j in  1:length(time_die))
       {
       if (j==1)
         {
           returnva_p7=pava_ifr(time_die,ttot,deaths)#returnva_p4=list(struct$time[j],struct$ttot[j],struct$deaths[j])
           m=dim(returnva_p7)[2]/3
           w1=length(returnva_p7)/3
           time_r[1:w1,j]=returnva_p7[,1:m]
           ttot_r[1:w1,j]=returnva_p7[,(m+1):2*m]
           deaths_r[1:w1,j]=returnva_p7[,(2*m+1):3*m]
          }
       if (j == length(time_die))
          {
          returnva_p=pava_dfr(time_die,ttot,deaths)#returnva_p4=list(struct$time[j],struct$ttot[j],struct$deaths[j])
          m=dim(returnva_p)[2]/3
          w2=length(returnva_p)/3
          time_r[1:w2,j]=returnva_p[,1:m]
          ttot_r[1:w2,j]=returnva_p[,(m+1):2*m]
          deaths_r[1:w2,j]=returnva_p[,(2*m+1):3*m]
          }
       if ((j!=1)&&j < length(time_die))
          {

           returnva_p8=pava_dfr(time_die[1:j],ttot[1:j],deaths[1:j])
           m=dim(returnva_p8)[2]/3
           w3=length(returnva_p8)/3
           tmp_time1=returnva_p8[,1:m]
           tmp_ttot1=returnva_p8[,(m+1):2*m]
           tmp_deaths1=returnva_p8[,(2*m+1):3*m]
           returnva_p9=pava_ifr(time_die[-(1:j)],ttot[-(1:j)],deaths[-(1:j)])
           m=dim(returnva_p9)[2]/3
           w4=length(returnva_p9)/3
           tmp_time2=returnva_p9[,1:m]
           tmp_ttot2=returnva_p9[,(m+1):2*m]
           tmp_deaths2=returnva_p9[,(2*m+1):3*m]
           time_r[1:(w3+w4),j] = c(tmp_time1,tmp_time2)
           ttot_r[1:(w3+w4),j]= c(tmp_ttot1,tmp_ttot2)
           deaths_r[1:(w3+w4),j] = c(tmp_deaths1,tmp_deaths2)
           }
        # compute the loglikelihood for struct(j);
         t=time_r[,j]
         t=t[t!=0]
         tt=ttot_r[,j]
         tt=tt[tt!=0]
         d=deaths_r[1:length(tt),j]
         loglik_r[1,j]= gamllik(t,tt,d,time_die, ttot, deaths)
      }
      label = "decreasing then increasing failure rate"
      # select the times using the largest loglikelihood
      loglik  = loglik_r[1,1]
      time2   = time_r[,1]
      indx    = 1
      for (j in 1:length(time_die))
         if (loglik < loglik_r[1,j])
           {
            time2  = time_r[,j]
            loglik = loglik_r[1,j]
            indx   = j
            }
    }
  #  time2
  #  loglik
  #  indx

  if (indi == 3)#failure rate increases then decreases.ttot/death U shape, peak point removed
    {
    for (j in 1:length(time_die))
      {
      if (j==1)
         {
          returnva_p7=pava_dfr(time_die,ttot,deaths)
          m=dim(returnva_p7)[2]/3
          w1=length(returnva_p7)
          time_r[1:w1,j]=returnva_p7[,1:m]
          ttot_r[1:w1,j]=returnva_p7[,(m+1):2*m]
          deaths_r[1:w1,j]=returnva_p7[,(2*m+1):3*m]
          }
       if (j == length(time_die))
          {
          returnva_p=pava_ifr(time_die,ttot,deaths)#returnva_p=list(struct$time[j],struct$ttot[j],struct$deaths[j])
          m=dim(returnva_p)[2]/3
          w2=length(returnva_p)
          time_r[1:w2,j]=returnva_p[,1:m]
          ttot_r[1:w2,j]=returnva_p[,(m+1):2*m]
          deaths_r[1:w2,j]=returnva_p[,(2*m+1):3*m]
          }
       if ((j!=1)&&(j < length(time_die)))
          {
           returnva_p8=pava_ifr(time_die[1:j],ttot[1:j],deaths[1:j])
           m=dim(returnva_p8)[2]/3
           w3=length(returnva_p8)/3
           tmp_time1=returnva_p8[,1:m]
           tmp_ttot1=returnva_p8[,(m+1):2*m]
           tmp_deaths1=returnva_p8[,(2*m+1):3*m]
           returnva_p9=pava_dfr(time_die[-(1:j)],ttot[-(1:j)],deaths[-(1:j)])
           m=dim(returnva_p9)[2]/3
           w4=length(returnva_p9)/3
           tmp_time2=returnva_p9[,1:m]
           tmp_ttot2=returnva_p9[,(m+1):2*m]
           tmp_deaths2=returnva_p9[,(2*m+1):3*m]
           time_r[1:(w3+w4),j] = c(tmp_time1,tmp_time2)
           ttot_r[1:(w3+w4),j]= c(tmp_ttot1,tmp_ttot2)
           deaths_r[1:(w3+w4),j] = c(tmp_deaths1,tmp_deaths2)
           }
       # compute the loglikelihood for struct(j)
       t=time_r[,j]
       t=t[t!=0]
       tt=ttot_r[,j]
       tt=tt[tt!=0]
       d=deaths_r[1:length(tt),j]
       loglik_r[1,j]= gamllik(t,tt,d,time_die, ttot, deaths)

       if ((j < length(time_die))&& (j>1))
         {
           time_r[length(tmp_time1),j]=NA
           ttot_r[-length(tmp_time1),j]=NA
           deaths_r[length(tmp_time1),j]=NA
          }
     }
      label = "increasing then decreasing failure rate"
      # select the times using the largest loglikelihood
      loglik  = loglik_r[1,1]
      time2   = time_r[,1]
      indx    = 1
      for (j in 1:length(time_die))
         if (loglik < loglik_r[1,j])
            {
            time2  = time_r[,j]
            loglik = loglik_r[1,j]
            indx   = j
            }
     }
  if (indi == 4)# failure rate decreases then increases ttot/death arc shape; peak point removed.
    {
    for (j in 1:length(time_die))
      {
      if (j==1)
         {
          returnva_p7=pava_ifr(time_die,ttot,deaths)
          m=dim(returnva_p7)[2]/3
          w1=length(returnva_p7)/3
          time_r[1:w1,j]=returnva_p7[,1:m]
          ttot_r[1:w1,j]=returnva_p7[,(m+1):2*m]
          deaths_r[1:w1,j]=returnva_p7[,(2*m+1):3*m]
          }
       if (j == length(time_die))
          {
          returnva_p=pava_dfr(time_die,ttot,deaths)#returnva_p=list(struct$time[j],struct$ttot[j],struct$deaths[j])
          m=dim(returnva_p)[2]/3
          w2=length(returnva_p)/3
          time_r[1:w2,j]=returnva_p[,1:m]
          ttot_r[1:w2,j]=returnva_p[,(m+1):2*m]
          deaths_r[1:w2,j]=returnva_p[,(2*m+1):3*m]
          }
       if (j<length(time_die)&&j!=1)
          {
           returnva_p8=pava_dfr(time_die[1:j],ttot[1:j],deaths[1:j])
           m=dim(returnva_p8)[2]/3
           w3=length(returnva_p8)/3
           tmp_time1=returnva_p8[,1:m]
           tmp_ttot1=returnva_p8[,(m+1):2*m]
           tmp_deaths1=returnva_p8[,(2*m+1):3*m]
           returnva_p9=pava_ifr(time_die[-(1:j)],ttot[-(1:j)],deaths[-(1:j)])
           m=dim(returnva_p9)[2]/3
           w4=length(returnva_p9)/3
           tmp_time2=returnva_p9[,1:m]
           tmp_ttot2=returnva_p9[,(m+1):2*m]
           tmp_deaths2=returnva_p9[,(2*m+1):3*m]
           time_r[1:(w3+w4),j] = c(tmp_time1,tmp_time2)
           ttot_r[1:(w3+w4),j]= c(tmp_ttot1,tmp_ttot2)
           deaths_r[1:(w3+w4),j] = c(tmp_deaths1,tmp_deaths2)
           }
       # compute the loglikelihood for struct(j)
       t=time_r[,j]
       t=t[t!=0]
       tt=ttot_r[,j]
       tt=tt[tt!=0]
       d=deaths_r[1:length(tt),j]
       loglik_r[1,j]= gamllik(t,tt,d,time_die, ttot, deaths)
       if (j < length(time_die) && j > 1)
         {
           time_r[length(tmp_time1),j]=NA
           ttot_r[-length(tmp_time1),j]=NA
           deaths_r[length(tmp_time1),j]=NA
          }
      }
      label = "decreasing then increasing failure rate"
      # select the times using the largest loglikelihood
      loglik  = loglik_r[1,1]
      time2   = time_r[,1]
      indx    = 1
      for (j in 1:length(time_die))
         if (loglik < loglik_r[1,j])
            {
             time2  = time_r[,j]
             loglik = loglik_r[1,j]
             indx   = j
             }
   }
  #  loglik
  #  indx
  time_r=time_r[!is.na(time_r)&time_r!=0]
  ttot_r=ttot_r[!is.na(ttot_r)&ttot_r!=0]
  deaths_r=deaths_r[!is.na(deaths_r)&deaths_r!=0]
  time2=time2[time2!=0&!is.na(time2)]
  struct=list(time_r,ttot_r,deaths_r,loglik_r)
  returnc=list(time2, struct, label, indx)
  return(returnc)
}