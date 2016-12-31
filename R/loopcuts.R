#using a loop format to find out the times to make the cuts
#from large pvalues to small pvalues and the lambdas in between
#
#Input:
#        time: A sequence of time
#        censor: and censor (0=censored and 1=not censored)
#        cuttimes: unique, sorted, possible times to make the cuts,
#             including 0 and the ending time
#       mono: indicate the type of the test
#            mono == 0: 2-sided hypothesis: H0:lam1=lam2; H1:lam1 \ne lam2
#                 == 1: 1-sides hypothesis: H0:lam1>=lam2; H1:lam1 < lam2
#                         decreasing failure rate constraint
#                 == 2: 1-sides hypothesis: H0:lam1<=lam2; H1:lam1 > lam2
#                         increasing failure rate constraint
#Output:
#        ts: the times where the cuts shall be made
#        pvalues: the p values for deleting each cutting times

#' Loopcuts
#'
#' The function uses a loop format to find out the times to make the cuts
#' from large pvalues to small pvalues and the lambdas in between
#' 
#' @usage loopcuts(time,censor,cuttimes,mono)
#' 
#' @param time a sequence of time
#' @param censor a vector indicating censored or not at the given times, 0 = censored; 1 = uncensored
#' @param cuttimes unique, sorted, possible times to make the cuts, including 0 and the ending time
#' @param mono 0: 2-sided hypothesis: H0: lam1 is equal to lam2; H1: lam1 is not equal to lam2
#' 
#'             1: 1-sided hypothesis: H0: lam1 is greater than or equal to lam2; H1: lam1 is less than lam2
#'             
#'             2: 1-sided hypothesis: H0: lam1 is less than or equal to lam2; H1: lam1 is greater than lam2

#'
#' @return
#' the times where the cuts shall be made
#' 
#' the p values for deleting each cutting times
#' @export
#'
#' @examples
loopcuts <- function(time,censor,cuttimes,mono)
{
  #Sort cuttimes and find cutlen
  cuttimes=unique(sort(cuttimes))
  cutlen=length(cuttimes)

  #prepare the death time, ttot, and the number of deaths.returnval_1 is the four returned valure from function totaltest
  returnval_1=totaltest(time,censor)
  m=dim(returnval_1)[2]/3
  time_die=returnval_1[,1:m]
  ttot=returnval_1[,(m+1):(2*m)]
  deaths=returnval_1[,(2*m+1):3*m]
  #Based on cuttimes, find out the corresponding totaltime on test and the number of deaths.

  #Set up ts and pvalues
  ts=array(0,c(cutlen,1))
  pvalues=array(0,c(cutlen,1))
  if (cutlen==1)
  {
    for(j in 1:length(time_die))
       if (time_die[j]==cuttimes)
        {
          ttot1 = sum(ttot[1:j])
          ttot2 = sum(ttot[-(1:j)])
          death1 = sum(deaths[1:j])
          death2 = sum(deaths[-(1:j)])
          returnval_2= exact_pvalue(ttot1,ttot2,death1,death2,mono)
          a2=returnval_2[1]
          p=returnval_2[2]
         }
    ts= cuttimes
    pvalues = p
    }
  if (cutlen!=1)
     {
     for(i in 1:cutlen)
      {
       if (i==1)
         {
          ###initialize the matrix storing ttots and deaths and p values
          ttotvec= array(0,c(cutlen+1,1))
          deavec= array(0,c(cutlen+1,1))
          allt= cuttimes
          pvalall=array(0,c(cutlen,1))
          for(ii in 1:cutlen)
             for(j in 1:length(time_die))
                 if (time_die[j]==cuttimes[ii])
                  {
                   if (ii==1)
                      {
                       ttotvec[ii] = sum(ttot[1:j])
                       deavec[ii]  = sum(deaths[1:j])
                       }
                    if(ii!=1)
                       {
                       ttotvec[ii] = sum(ttot[1:j])-sum(ttotvec[1:(ii-1)])
                       deavec[ii]=sum(deaths[1:j])-sum(deavec[1:(ii-1)])
                       }
                    }
         #do the last item in ttotvec and deavec
         ttotvec[cutlen+1] = sum(ttot)-sum(ttotvec[1:cutlen])
         deavec[cutlen+1]  = sum(deaths)-sum(deavec[1:cutlen])
         #Compute ts(1) and pvalues(1)
         for (k in 1:length(allt))
           {
           ttot1 = ttotvec[k]
           ttot2 = ttotvec[k+1]
           dea1  = deavec[k]
           dea2  = deavec[k+1]
           returnval_3=exact_pvalue(ttot1,ttot2,dea1,dea2,mono)
           a2=returnval_3[1]
           pvalall[k]=returnval_3[2]
            }
         maxp = max(pvalall)
         maxplab = 0
        for (k1 in 1:length(pvalall))
          if (maxp == pvalall[k1])
             maxplab = k1
        #pvalall
        #save ts(i) and pvalues(i)
        ts[i] = allt[maxplab]
        pvalues[i]= pvalall[maxplab]
        allt=allt[-maxplab]
        pvalall=pvalall[-maxplab]
        #merge centain items on ttotvec and deavec
        ttotvec[maxplab] = ttotvec[maxplab]+ ttotvec[maxplab+1]
        ttotvec=ttotvec[-(maxplab+1)]
        deavec[maxplab]= deavec[maxplab]+ deavec[maxplab+1]
        deavec=deavec[-(maxplab+1)]
       }
      if (i!=1)
       {
        #Use position maxplab, two vectors:(allt,pvalall), and ttotvec, deavec to compute ts(i) and pvalues(i)
        maxplab1=maxplab
        n=length(allt)
        if(maxplab1==1)
         {
           #p = lrtpvalue(ttotvec(1),ttotvec(2),deavec(1),deavec(2),n)
           returnval_4= exact_pvalue(ttotvec[1],ttotvec[2],deavec[1],deavec[2],mono)
           a2=returnval_4[1]
           p=returnval_4[2]
           if (length(deavec)==2)
             {
               pvalues[i]=p
               ts[i]=allt[maxplab]
              }
           else
              {
               pvalall[1]= p
               maxp=max(pvalall)
               maxplab= 0
               for (k1 in 1:length(pvalall))
                   if (maxp == pvalall[k1])
                        maxplab= k1
                #save ts(i) and pvalues(i)
                ts[i] = allt[maxplab]
                pvalues[i] = pvalall[maxplab]
                allt=allt[-maxplab]
                #pvalall
                pvalall=pvalall[-maxplab]
                #merge centain items on ttotvec and deavec
                ttotvec[maxplab] = ttotvec[maxplab] + ttotvec[maxplab+1]
                ttotvec=ttotvec[-(maxplab+1)]
                deavec[maxplab]  = deavec[maxplab] + deavec[maxplab+1]
                deavec=deavec[-(maxplab+1)]
                }
           }
         if((maxplab1!=1)&(maxplab1>n))
          {
            returnval_5=exact_pvalue(ttotvec[maxplab-1],ttotvec[maxplab],deavec[maxplab-1],deavec[maxplab],mono)
            a2=returnval_5[1]
            p=returnval_5[2]
            pvalall[length(allt)]= p
            maxp= max(pvalall)
            maxplab=0
            for (k1 in 1:length(pvalall))
              if (maxp == pvalall[k1])
                  maxplab= k1
            #save ts(i) and pvalues(i)
            ts[i] = allt[maxplab]
            pvalues[i] = pvalall[maxplab]
            allt=allt[-maxplab]
            #pvalall
            pvalall=pvalall[-maxplab]
            #merge centain items on ttotvec and deavec
            ttotvec[maxplab] = ttotvec[maxplab]+ ttotvec[maxplab+1]
            ttotvec=ttotvec[-(maxplab+1)]
            deavec[maxplab] = deavec[maxplab] + deavec[maxplab+1]
            deavec=deavec[-(maxplab+1)]
            }
         if ((maxplab1!=1)&&(maxplab1<=n))
          {
          returnval_6=exact_pvalue(ttotvec[maxplab-1],ttotvec[maxplab],deavec[maxplab-1],deavec[maxplab],mono)
          a2=returnval_6[1]
          pfront=returnval_6[2]
          returnval_7=exact_pvalue(ttotvec[maxplab],ttotvec[maxplab+1],deavec[maxplab],deavec[maxplab+1],mono)
          a2=returnval_7[1]
          pback=returnval_7[2]
          pvalall[maxplab-1]= pfront
          pvalall[maxplab]= pback
          maxp= max(pvalall)
          maxplab= 0
          for (k1 in 1:length(pvalall))
             if (maxp==pvalall[k1])
                maxplab=k1
          #save ts(i) and pvalues(i)
          ts[i]= allt[maxplab]
          pvalues[i]= pvalall[maxplab]
          allt=allt[-maxplab]
          #pvalall
          pvalall=pvalall[-maxplab]
          #merge centain items on ttotvec and deavec
          ttotvec[maxplab]= ttotvec[maxplab]+ ttotvec[maxplab+1]
          ttotvec=ttotvec[-(maxplab+1)]
          deavec[maxplab]= deavec[maxplab]+ deavec[maxplab+1]
          deavec=deavec[-(maxplab+1)]
         }
       }
      }
     }
    returnv=cbind(ts,pvalues)
    return(returnv)
}

