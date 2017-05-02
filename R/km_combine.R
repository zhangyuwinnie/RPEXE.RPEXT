#' Comparing two Kaplan Meier curves in one plot
#' 
#' @description The function compares two Kaplan Meier curves in one plot
#'
#' @param x1 Nx2 data matrix,first columen represents survival time of the i-th subject, second column represents censored flag (0 if not censored, 1 if censored)
#' @param x2 Nx2 data matrix,first columen represents survival time of the i-th subject, second column represents censored flag (0 if not censored, 1 if censored)
#' @param pos The position of the legend. Can be 0 or 1. The legend will be 
#'    on the topright if set to 0. The legend will be on the bottomleft if set to 1. Default is 0.
#'
#' @return
#' A combined Kaplan Meier curve 
#' 
#' @export
#'
#' @examples
#' t1 <- c(2,3,4,5.5,7,10,12,15)
#' c1 <- c(0,0,1,0,0,1,0,0)
#' t2 <- c(1,3,5,4,8,10,9,11)
#' c2 <- c(0,0,0,0,1,0,0,0)
#' x1<-cbind(t1,c1)
#' x2<-cbind(t2,c2)
#' km_combine(x1,x2)
#' km_combine(x1,x2,pos=1)
km_combine <- function(x1, x2, pos = 0){

  
  kmout1 <- kmvalue(x1)
  kmout2 <- kmvalue(x2)
  table1 = kmout1$table1
  table12 = kmout1$table12
  t1 = kmout1$t1
  T1 = kmout1$T1
  xcg1 = kmout1$xcg
  ycg1 = kmout1$ycg
  t2 = kmout2$t1
  T2 = kmout2$T1
  xcg2 = kmout2$xcg
  ycg2 = kmout2$ycg

  
  # plot Kaplan_Meier
  plot(t1, T1, type = "s",xlab="Time", ylab="Estimated survival functions", lty = 2, 
       col = "blue",ylim = c(0,1))
  points(xcg1, ycg1,cex=.8,pch=3, col = "black")
  
  lines(t2, T2, type = "s",xlab="Time", ylab="Estimated survival functions", lty = 2, 
        col = "red",ylim = c(0,1))
  points(xcg2, ycg2,cex=.8,pch=3, col = "black")
  title('Kaplan-Meier estimate of survival functions')
  if (pos == 0)
    position = "topright"
  if (pos == 1)
    position = "bottomleft"
  legend(position,legend=c("Treatment 1","Treatment 2","Censored"),col=c("blue", "red","black"), lty=c(2,2,NA),pch = c(NA,NA,3), cex=0.8)
  
  
  
}
