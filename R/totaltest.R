#Function 'totaltest' computes total time on test
#Inputs: times and censor (0 = censored; 1 = uncensored)
#Return: time_die: times that events occur (in ascending order)
#        ttot: total time on test at each time point in time_die
#        deaths: number of death at each time point in time_die


#' Computes total time on test
#'
#' @param time 
#' @param censor 
#'
#' @usage totaltest(time,censor)
#' @return
#' time_die
#' @export
#'
#' @examples
totaltest <- function(time,censor)
{
  #save the data in the orginal structure
  sew=array(0,c(length(time),1))
  for(i in 1:length(time))
     sew[i]=i
  tmpdata=array(cbind(sew,censor,time),c(length(sew),10))

  #sort the tmp data
  tmp2 =array(0,c(dim(tmpdata)))
  for (i in 1:dim(tmp2)[1])
     {
      if (i!=dim(tmp2)[1])
      {
      tmp2[i,3]=min(tmpdata[,3])
      j=which.min(tmpdata[,3])
      tmp2[i,1:2]=tmpdata[j,1:2]
      tmpdata=tmpdata[-j,]
       }
     else
       tmp2[i,]=tmpdata
     }
  #Compute alpha's for the sequence
  for (i in 1:dim(tmp2)[1])
     {
     if (tmp2[i,2]==0)
        tmp2[i,4]=0
     else
        tmp2[i,4]=1
      }
  #Deal with alpha > 1
  for (i in 1:(dim(tmp2)[1]-1))
      if (tmp2[dim(tmp2)[1]+1-i,3]== tmp2[dim(tmp2)[1]-i,3])
        {
         tmp2[dim(tmp2)[1]-i,4]= tmp2[dim(tmp2)[1]-i,4] + tmp2[dim(tmp2)[1]+1-i,4]
         tmp2[dim(tmp2)[1]+1-i,4]=0
        }
  # Delete the repeats
  k=as.null()
  for (i in 1:dim(tmp2)[1])
      if (tmp2[i,2] == 1&tmp2[i,4]==0)
          k[length(k)+1]=i
  tmp3 = tmp2
  if (length(k)!=0)
     tmp3=tmp3[-k,]

  #Compute the number of patients in the study
  for(i in 1:dim(tmp3)[1])
    {
     if (tmp3[i,4]==0)
         tmp3[i,5]= 1
     else
         tmp3[i,5]= tmp3[i,4]
     }
  for(i in 1:dim(tmp3)[1]-1)
    {
     tmp3[i,6]= sum(tmp3[,5])-sum(tmp3[1:i,5])
     tmp3[dim(tmp3)[1],6]= 0
    }

  #Compute the survival time of this cohort
  for(i in 1:dim(tmp3)[1])
     {
     if (i==1)
        tmp3[i,7]= sum(tmp3[,5])*tmp3[i,3]
     else
        ###survival time == patient number * time difference
        tmp3[i,7] = tmp3[i-1,6]* (tmp3[i,3]-tmp3[i-1,3])
     }
  tmp3[,8] = tmp3[,7]
  for (i in 1:dim(tmp3)[1])
     if (tmp3[i,2]==0)
       {
         if (t(tmp3[i:dim(tmp3)[1],2])%*%tmp3[i:dim(tmp3)[1],2]>0)
             tmp3[i+1,8] = tmp3[i,8]+tmp3[i+1,8]
         if (t(tmp3[i:dim(tmp3)[1],2])%*%tmp3[i:dim(tmp3)[1],2]==0&tmp3[i-1,2]!=0)
               {
              ### put all the credit to the last noncensered data
               k = length(tmp3[i:dim(tmp3)[1],2])
               for (j in 1:k)
                  tmp3[i-1,8] = tmp3[i-1,8]+tmp3[i-1+j,8]
                }
         }

  #Build the survival reaction
  tmp3[,9] = tmp3[,6]
  for (i in 2:length(tmp3[,9]))
     tmp3[i,9]= tmp3[i-1,9]-tmp3[i,4]

  #plot (tmp3[,3],tmp3[,9])

  ###delete all the censered items
  k=as.null()
  for (i in 1:length(tmp3[,1]))
     if (tmp3[i,2]== 0)
        k[length(k)+1]=i
  tmp4 = tmp3
  if (length(k)!=0)
     tmp4=tmp4[-k,]

  time_die=tmp4[,3]
  ttot=  tmp4[,8]
  deaths=  tmp4[,5]
  returnv=cbind(time_die,ttot,deaths)
  return(returnv)
}
