#' @title total time on test
#' 
#' @description Function 'totaltest' computes total-time-on-test.
#'
#' @param time event/censoring times
#' @param censor censoring status
#'
#' @usage totaltest(time,censor)
#' @return
#' time_die time points where events occur (in ascending order)
#' ttot total time on test corresponding to each time point in "time_die"
#' deaths number of death corresponding to each time point in "time_die"
#' @export
#'
#' @examples
#' t1 <- c(2,3,4,5.5,7,10,12,15)
#' c1 <- c(0,0,1,0,0,1,0,0)
#' totaltest(t1,c1)
totaltest <- function(time,censor)
{
  #save the data in the orginal structure
  sew=rep(0, length(time))
  for(i in 1:length(time))
    sew[i]=i
  tmpdata=cbind(sew,censor,time)
  
  # sort the tmp data
  tmp2 = tmpdata[order(time),]
  tmp2 = cbind(tmp2,rep(0,dim(tmp2)[1]))
  
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
  tmp3 = cbind(tmp3,rep(0,dim(tmp3)[1]))
  
  #Compute the number of patients in the study
  for(i in 1:dim(tmp3)[1])
  {
    if (tmp3[i,4]==0)
      tmp3[i,5]= 1
    else
      tmp3[i,5]= tmp3[i,4]
  }
  tmp3 = cbind(tmp3,rep(0,dim(tmp3)[1]))
  for(i in 1:dim(tmp3)[1]-1)
  {
    tmp3[i,6]= sum(tmp3[,5])-sum(tmp3[1:i,5])
    tmp3[dim(tmp3)[1],6]= 0
  }
  
  #Compute the survival time of this cohort
  tmp3 = cbind(tmp3,rep(0,dim(tmp3)[1]))
  for(i in 1:dim(tmp3)[1])
  {
    if (i==1)
      tmp3[i,7]= sum(tmp3[,5])*tmp3[i,3]
    else
      ###survival time == patient number * time difference
      tmp3[i,7] = tmp3[i-1,6]* (tmp3[i,3]-tmp3[i-1,3])
  }
  tmp3 = cbind(tmp3,rep(0,dim(tmp3)[1]))
  tmp3[,8] = tmp3[,7]
  for (i in 1:dim(tmp3)[1])
    if (tmp3[i,2]==0)
    {
      if (t(tmp3[i:dim(tmp3)[1],2])%*%tmp3[i:dim(tmp3)[1],2]>0)
        tmp3[i+1,8] = tmp3[i,8]+tmp3[i+1,8]
      if (t(tmp3[i:dim(tmp3)[1],2])%*%tmp3[i:dim(tmp3)[1],2]==0 && tmp3[i-1,2]!=0)
      {
        ### put all the credit to the last noncensered data
        k = length(tmp3[i:dim(tmp3)[1],2])
        for (j in 1:k)
          tmp3[i-1,8] = tmp3[i-1,8]+tmp3[i-1+j,8]
      }
    }
  
  #Build the survival reaction
  tmp3 = cbind(tmp3,rep(0,dim(tmp3)[1]))
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
