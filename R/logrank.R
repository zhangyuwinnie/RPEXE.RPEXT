# LOGRANK Comparing survival curves of two groups using the log rank test
# Comparison of two survival curves can be done using a statistical
# hypothesis test called the log rank test. It is used to test the null
# hypothesis that there is no difference between the population survival
# curves (i.e. the probability of an event occurring at any time point is
# the same for each population). This function use the Kaplan-Meier
# procedure to estimate the survival function, so it is mandatory to download
# KMPLOT (http://www.mathworks.com/matlabcentral/fileexchange/22293).
#
# Syntax: 	logrank(x1,x2,alpha,censflag)
#      
#     Inputs:
#           X1 and X2 (mandatory)- Nx2 data matrix:
#                     (X:,1) = survival time of the i-th subject
#                     (X:,2) = censored flag 
#                             (0 if not censored; 1 if censored)
#           
#      Outputs:
#           Kaplan-Meier plot
#           Log-rank statistics
#
#      Example: 
#           load logrankdata x1 x2
#           logrank(x1,x2)
#
# LOG-RANK TEST FOR KAPLAN-MEIER SURVIVAL FUNCTIONS
#
# --------------------------------------------------------------------------------
# UL				S.E.			z				p-value (2-tailed test) 	alpha
# --------------------------------------------------------------------------------
# 6.57226		2.80788			2.16258			0.03057                     0.050
# --------------------------------------------------------------------------------
#    		The survival functions are statistically different
#
#            Created by Giuseppe Cardillo
#           giuseppe.cardillo-edta@poste.it
#
# To cite this file, this would be an appropriate format:
# Cardillo G. (2008). LogRank: Comparing survival curves of two groups
# using the log rank test
# http://www.mathworks.com/matlabcentral/fileexchange/22317

#' @title Obtain the logrank test 
#' 
#' @description 
#' Comparing survival curves of two groups using the log rank test
#' Comparison of two survival curves can be done using a statistical
#' hypothesis test called the log rank test. It is used to test the null
#' hypothesis that there is no difference between the population survival
#' curves (i.e. the probability of an event occurring at any time point is
#' the same for each population). 
#' 
#' @usage logrank(x1, x2, Alpha)
#'
#' @param x1 Nx2 data matrix,first columen represents survival time of the i-th subject, second column represents censored flag (0 if not censored, 1 if censored)
#' @param x2 Nx2 data matrix,first columen represents survival time of the i-th subject, second column represents censored flag (0 if not censored, 1 if censored)
#' @param Alpha  Significance level(optional, default is 0.05)
#'
#' @return Log-rank statistics
#' @export
#'
#' @examples
#' logrank_Yates228(x1,x2,0.05)
#' logrank_Yates228(x1,x2)
logrank <- function(x1, x2, Alpha=NA,alpha0){
  #     Inputs:
  #           X1 and X2 (mandatory)- Nx2 data matrix:
  #                     (X:,1) = survival time of the i-th subject
  #                     (X:,2) = censored flag 
  #                             (0 if not censored; 1 if censored)
  
  #reset parameter inputs
  if (!is.na(Alpha))
  {
    alpha = alpha0
  } else{
    alpha = 0.05
  }
  
  kmout1 <- kmplot(x1, 0.05)
  kmout2 <- kmplot(x2, 0.05)
  table1 = kmout1$table1
  table12 = kmout1$table12
  t1 = kmout1$t1
  T1 = kmout1$T1
  xcg1 = kmout1$xcg
  ycg1 = kmout1$ycg
  lambda1 = kmout1$lambda
  table2 = kmout2$table1
  table22 = kmout2$table12
  t2 = kmout2$t1
  T2 = kmout2$T1
  xcg2 = kmout2$xcg
  ycg2 = kmout2$ycg
  lambda2 = kmout2$lambda
  
  # plot Kaplan_Meier
  plot(t1, T1, type = "s",xlab="Time", ylab="Estimated survival functions", lty = 2, 
       col = "blue",ylim = c(0,1))
  points(xcg1, ycg1,cex=.8,pch=3, col = "black")
  
  lines(t2, T2, type = "s",xlab="Time", ylab="Estimated survival functions", lty = 2, 
       col = "red",ylim = c(0,1))
  points(xcg2, ycg2,cex=.8,pch=3, col = "black")
  title('Kaplan-Meier estimate of survival functions')
  legend("bottomleft",legend=c("Treatment 1","Treatment 2","Censored"),col=c("blue", "red","black"), lty=c(2,2,NA),pch = c(NA,NA,3), cex=0.8)
  
  # Full-blown LOGRANK procedure
  # Merge the first columns of Table1 and Table2 (time intervals)
  # and pick-up unique values
  A = sort(unique(c((table1[,1]),(table2[,1]))))
  table = matrix(0, length(A), 9)
  # Out in the first column the time intervals
  table[,1]=A 
  # Put in the columns 2 and 3 and in the proper rows the deaths and alive
  # taken from table1 columns 2 and 3
  ib = match(table1[,1] ,A)
  # remove NA
  ib = ib[!is.na(ib)]
  ia = match(A, table1[,1])
  ia = ia[!is.na(ia)]
  table[ib,2:3] = table1[ia,2:3]
  # Put in the columns 4 and 5 and in the proper rows the deaths and alive
  # taken from table2 columns 2 and 3
  ib = match(table2[,1] ,A)
  ib = ib[!is.na(ib)]
  ia = match(A, table2[,1])
  ia = ia[!is.na(ia)]
  table[ib,4:5] = table2[ia,2:3]
  # remove the rows where there arent't deaths in both treatments
  table = table[-(which(table[,2]==0 & table[,4]==0)),]
  # fill the "pigeon-holes"
  c= which(table[,3]==0) # find the "pigeon-holes" of treatment 1
  for (I in 1:length(c))
  {
    if (c[I] != 1)
    {
      # find the first interval time before the hole where there is almost 1 death
      J = max(which(table[(1:(c[I]-1)),3]>0))
      table[c[I],3] = table[J,3] - table[J,2]
      # find eventually censored data
      m = which(table12[,1]<table[c[I],1] & table12[,1]>=table[J,1])
      if (length(m) != 0){
        K = max(m)
      } else{
        K = NA
      }
      # Put in the hole how many subject were alive before the interval time of the hole
      if (!is.na(K))
      {
        table[c[I],3] = table[c[I],3] - sum(table12[K,2])
      }
    } else {
      table[1,3] = dim(x1)[1]
    }
  }
  # Do the same for tratment 2
  c= which(table[,5]==0) # find the "pigeon-holes" of treatment 1
  for (I in 1:length(c))
  {
    if (c[I] != 1)
    {
      # find the first interval time before the hole where there is almost 1 death
      J = max(which(table[(1:(c[I]-1)),5]>0))
      table[c[I],5] = table[J,5] - table[J,4]
      # find eventually censored data
      m = which(table22[,1]<table[c[I],1] & table22[,1]>=table[J,1])
      if (length(m) != 0){
        K = max(m)
      } else{
        K = NA
      }
      # Put in the hole how many subject were alive before the interval time of the hole
      if (!is.na(K))
      {
        table[c[I],5] = table[c[I],5] - sum(table22[K,2])
      }
    } else {
      table[1,5] = dim(x2)[1]
    }
  }
  
  # Fill the table and compute the statistic variable
  # Compute the total deaths and alive before the i-th time interval
  table[,6]=table[,2]+table[,4]
  table[,7]=table[,3]+table[,5]
  # Compute the difference between observed deaths for treatment 1 and
  # expected deaths in the hyphthesis that the treatments are similar
  table[,8]=table[,2]-table[,3]*table[,6]/table[,7]
  # Log-rank statistic is the sum of column 8 values
  UL=abs(sum(table[,8]))
  # Compute the contribute to the standard error
  table[,9]=table[,3]*table[,5]*table[,6]*(table[,7]-table[,6])/(table[,7]^2*(table[,7]-rep(1,dim(table)[1])))
  # find if there is some NaN (i.e. 0/0)
  loc = is.na(table[,9])
  if (any(loc))
  {
    table[loc,9] = 0
  }
  SUL = sqrt(sum(table[,9])) #Compute the total standard error
  z = abs(UL/SUL) # normalized UL with Yates'es correction
  p= 2*(1-pnorm(z,0,1))
  
  # display results
  line1 = " LOG-RANK TEST FOR KAPLAN-MEIER SURVIVAL FUNCTIONS\n"
  line2 = "------------------------------------------------------------------------------------------\n"
  line3 = "HAZARD RATE IS AN EXPERIMENTAL FUNCTION!!!!\n"
  line4 = paste("Treatment 1: Hazard rate: ", round(lambda1,4),"\n")
  line5 = paste("Treatment 2: Hazard rate: ", round(lambda2,4),"\n")
  line6 = paste('Hazard ratio: ', round(lambda1/lambda2,4),"\n")
  line7 = 'UL\t\tS.E.\t\tz\t\tp-value (2-tailed test)\t\talpha\n'
  line8 = paste(round(UL,5),'\t',round(SUL,5),'\t',round(z,5),'\t',round(p,5),'\t\t\t',round(alpha,3),'\n')
  cat(line1, line2, line3, line4, line5, line6, line2, line7, line2,line8,line2)
  
  if (p<alpha){
    cat('\t\tThe survival functions are statistically different\n')
  } else{
    cat('\t\tThe survival functions are not statistically different\n')
  }
}
#x1= cbind(time_p, cens_p)
#x2 = cbind(time_n, cens_n)
#logrank(x1,x2,Alpha = "alpha",0.01)


