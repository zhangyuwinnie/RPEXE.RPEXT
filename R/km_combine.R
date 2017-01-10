#' Comparing two Kaplan Meier curves in one plot
#' 
#' @description The function...
#'
#' @param x1 
#' @param x2 
#'
#' @return
#' A combined Kaplan Meier curve
#' 
#' @export
#'
#' @examples
km_combine <- function(x1,x2){

  
  kmout1 <- kmplot(x1, 0.05)
  kmout2 <- kmplot(x2, 0.05)
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
  legend("bottomleft",legend=c("Treatment 1","Treatment 2","Censored"),col=c("blue", "red","black"), lty=c(2,2,NA),pch = c(NA,NA,3), cex=0.8)
  
  
  
}
