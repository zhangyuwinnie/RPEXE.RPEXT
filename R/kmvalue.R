#' Obtain values for Kaplan-Meier plotting
#'
#' @param x Nx2 data matrix,first columen represents survival time of the i-th subject, second column represents censored flag (0 if not censored, 1 if censored)
#'
#' @usage kmvalue(x)
#' @return
#' Values used for Kaplan-Meier plotting
#' @export
#'
#' @examples
#' t1 <- c(2,3,4,5.5,7,10,12,15)
#' c1 <- c(0,0,1,0,0,1,0,0)
#' x1<-cbind(t1,c1)
#' kmvalue(x1)
kmvalue <- function(x){
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
  # list of censored data, "0" represent censored
  cenData = x[which(x[,2]==0)]
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
    B = min(which(t1>=table12[I,1]))
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
