#' @title PAVA order restriction under decreasing failure rate (DFR)
#' 
#' @description This function imposes the PAVA DFR order restriction by eliminating change-points violating the restriction 
#' 
#' @usage pava_dfr(time_die,ttot,deaths)
#' 
#' @param time_die event times
#' @param ttot the total time on test (ttot) corresponding to the event times
#' @param deaths the number of deaths at each event time
#' 
#' @return
#' time2: the event times after PAVA
#' ttot2: the corresponding ttot
#' deaths2 the corresponding number of deaths
#' 
#' @export
#'
#' @examples
#' data(pava_dfrd)
#' t_d = pava_dfrd[,1]
#' t = pava_dfrd[,2]
#' d = pava_dfrd[,3]
#' pava_dfr(t_d, t, d)
pava_dfr <- function(time_die,ttot,deaths)
{
  len=length(ttot)
  if (len == 1)
  {
     time2   = time_die
     ttot2   = ttot
     deaths2 = deaths
  }
  if (len!=1)
  {
   for (j in 1:len)
   {
      #indicate deleted items
     n=length(ttot)
     if(n!=1)
     {
      for (i in 1:(n-1))
         #indicate items to be deleted
         if ((ttot[i]/deaths[i])>(ttot[i+1]/deaths[i+1]))
         {
             deaths[i+1] = deaths[i]+deaths[i+1]
             deaths[i]= 0
             ttot[i+1]=ttot[i]+ttot[i+1]
             ttot[i]=0
         }
       #delete the indicated items
       k=as.null()
       for (i in 1:length(ttot))
          if(ttot[i]==0)
            k[length(k)+1]=i
      if (!is.null(k))
        {
         ttot=ttot[-k]
         deaths=deaths[-k]
         time_die=time_die[-k]
         }
      }
     }
   }
   time2   = time_die
   ttot2   = ttot
   deaths2 = deaths
   returnval=cbind(time2,ttot2,deaths2)
   return(returnval)
}
