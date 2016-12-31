# KMPLOT Plot the Kaplan-Meier estimation of the survival function
# Survival times are data that measure follow-up time from a defined
# starting point to the occurrence of a given event, for example the time
# from the beginning to the end of a remission period or the time from the
# diagnosis of a disease to death. Standard statistical techniques cannot
# usually be applied because the underlying distribution is rarely Normal
# and the data are often "censored". A survival time is described as
# censored when there is a follow-up time but the event has not yet
# occurred or is not known to have occurred. For example, if remission time
# is being studied and the patient is still in remission at the end of the
# study, then that patients remission time would be censored. If a patient
# for some reason drops out of a study before the end of the study period,
# then that patients follow-up time would also be considered to be
# censored. The survival function S(t) is defined as the probability of
# surviving at least to time t. The graph of S(t) against t is called the
# survival curve. The Kaplan meier method can be used to estimate this
# curve from the observed survival times without the assumption of an
# underlying probability distribution.
#
# Syntax: 	kmplot(x,alpha)
#      
#     Inputs:
#           X (mandatory)- Nx2 data matrix:
#                          (X:,1) = survival time of the i-th subject
#                          (X:,2) = censored flag 
#                                   (0 if not censored; 1 if censored)
#           ALPHA  - significance level (default 0.05) 
#           
#     Outputs:
#           Kaplan-Meier plot
  
#' Obtain values for Kaplan-Meier plotting
#'
#' @param x Nx2 data matrix,first columen represents survival time of the i-th subject, second column represents censored flag (0 if not censored, 1 if censored)
#' @param alpha  significance level
#'
#' @usage kmplot(x, alpha)
#' @return
#' Values used for Kaplan-Meier plotting
#' @export
#'
#' @examples
kmplot <- function(x, alpha){
  # string for LEGEND function
  str1 = paste(toString((1-alpha)*100),"% confidence interval",sep = "")
  # sort data by survival time
  x = x[order(x[,1]),]
  # table of patients observed for each survival time
  # the TABULATE function sets up this matrix:
  # table1=[time count]
  # get frequent list for survival time
  frequent <- table(x[,1])
  # convert frequency list to matrix
  y <- matrix(c(as.numeric(names(frequent)), frequent), ncol=2, byrow=FALSE, 
                  dimnames=NULL)
  # add (0,size) to the first row of matrix y
  table1 = rbind(c(0,sum(y[,2])),y)
  
  # Table of censored data
  # list of censored data, "1" represent censored
  cenData = x[which(x[,2]==1)]
  freqCen <- table(cenData)
  # convert censor data frequency list to matrix
  table12 <- matrix(c(as.numeric(names(freqCen)), freqCen), ncol=2, byrow=FALSE, 
              dimnames=NULL)
  # setup the vector of the censored data (function is.element() return T or F whether elements in first 
  # vector contained in second vector
  cens <- as.numeric(is.element(table1[,1],table12[,1]))
  # loc return index in vector 2 of matched element  
  loc <- match(table1[,1],table12[,1])
  # replace NA with 0
  loc[is.na(loc)] <- 0
  
  # place in the third column how many subjects are still alive at the
  # beginning of the i-th interval.
  a1 = c(table1[1,2], -1*table1[-1,2])
  table1 <- cbind(table1,cumsum(a1))
  end1 <- dim(table1)[1]
  table1[-1,3] = table1[-end1,3]
  # number of deaths in the intervals (don't take in account the censored
  # data)
  table1[which(cens==1),2] = table1[which(cens==1),2]-table12[loc[which(cens == 1)],2]
  # finally, delete the first row that is now useless
  table1<- table1[-1,]
  
  # this is the x variable (time);
  t1 = c(0,table1[,1])
  # this is the y variable (survival function)
  T1 = c(1, cumprod(( 1-(table1[,2]/table1[,3])   )))
 # censored data plotting
  
 # if there are censored data after max(t1), add a new cell into the t1,T1 
 end12 <- dim(table12)[1]
 if (table12[end12, 1] >= t1[length(t1)]){
   t1[length(t1)+1] = table12[end12,1]+1
   T1[length(T1)+1] = T1[length(T1)]
 }

 # data of censored points
 # vectors preallocation
  xcg = rep(0, sum(table12[,2]))
  ycg = xcg
  J = 1
  # for each censored data into the i-th time interval...
  for (I in 1:end12)
  { 
    # compute how many position into the array they must occupy
    JJ = J + table12[I,2]-1
    #find the correct time interval in which censored data must be placed
    B = min(which(t1>table12[I,1]))
    A = B-1
    # equally divide this interval
    int = seq(from = table12[I,1], to = t1[B], len = table12[I,2]+2)
    xcg[J:JJ] = int[2:(length(int)-1)]
    ycg[J:JJ] = T1[A]
    # update the counter
    J = JJ + 1
  }

 
 # compute the hazard rate
 c1 = T1*dim(x)[1]*dim(x)[2]
 c2 = -(diff(log(c1[-length(c1)]))/diff(t1[-length(t1)]))
 lambda = mean(c2[which(c2 != 0)])
 
 # output
 kmout = list("table1"=table1, "table12" = table12, "t1" = t1, "T1" = T1, "xcg" = xcg, "ycg" = ycg, 
              "lambda" = lambda)
 return (kmout)
}
#x = cbind(time_p,cens_p)
