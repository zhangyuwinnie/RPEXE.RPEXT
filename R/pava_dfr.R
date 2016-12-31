#merge certain entries to make the sequence of ttot to be non decreasing
#Input:
#      time_die == a sequence of times where deaths happened.
#      ttot     == the total time on test at each time point.
#      deaths   == the number of deaths at each time point.
#Output:
#      time2    == the merged time_die
#      ttot2    == .......... ttot
#      deaths2  == .......... deaths

#' Pava_dfr
#'
#' Merge certain entries to make the sequence of ttot to be non decreasing
#'
#'
#' @param time_die 
#' @param ttot 
#' @param deaths 
#'
#' @usage pava_dfr(time_die,ttot,deaths)
#' 
#' @return
#' time2: the merged time_die
#' ttot2: the merged ttot
#' deaths2 the merged deaths
#' 
#' @export
#'
#' @examples
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
