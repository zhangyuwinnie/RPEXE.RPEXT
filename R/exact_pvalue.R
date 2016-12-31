#Compute the exact P value from the likelihood ratio test
#Input:
#       ttot1, ttot2 = total time on test 1 and 2
#       dea1, dea2   = number of death 1 and 2
#       mono: indicate the type of the test
#           mono == 0: 2-sided hypothesis: H0:lam1=lam2; H1:lam1 \ne lam2
#                == 1: 1-sided hypothesis: H0:lam1>=lam2; H1:lam1 < lam2
#                == 2: 1-sided hypothesis: H0:lam1<=lam2; H1:lam1 > lam2
# Output:
#       pval         = the exact Pvalue for testing H0: lambda1 = lambda2
# Assumption:
#       ttot1 and ttot2 are from Gamma(dea1,lambda1) and
#       Gamma(dea2,lambda2)
# Compute the test statistic, delta = (ttot1/(ttot1+ttot2))^dea1*...
#                                     (ttot2/(ttot1+ttot2))^dea2;

#' @title Compute the exact P value from the likelihood ratio test
#' 
#' @description The function...
#'
#' @param ttot1 total time on test 1
#' @param ttot2 total time on test 2
#' @param dea1 number of death 1
#' @param dea2 number of death 2
#' @param mono 0: 2-sided hypothesis: H0: lam1 is equal to lam2; H1: lam1 is not equal to lam2
#' 
#'             1: 1-sided hypothesis: H0: lam1 is greater than or equal to lam2; H1: lam1 is less than lam2
#'             
#'             2: 1-sided hypothesis: H0: lam1 is less than or equal to lam2; H1: lam1 is greater than lam2
#'
#' @usage exact_pvalue(ttot1,ttot2,dea1,dea2,mono)
#' 
#' @return
#' a2
#' pval
#' 
#' @export
#'
#' @examples
exact_pvalue <- function(ttot1,ttot2,dea1,dea2,mono)
{


  a1=ttot1/(ttot1+ttot2)
  delta=(log(a1)*dea1)+(log(1-a1)*dea2)

  #compute the p-value
  if(dea1==1&dea2==1)
  {
    # if dea1 = dea2 = 1
    if (mono==0)
    {
       if(a1<0.5)
          pval = 2*pbeta(a1,dea1,dea2)
       else
          pval = 2*(1-pbeta(a1,dea1,dea2))
    } else if (mono==1) {
       pval = pbeta(a1,dea1,dea2)
       pval = pval/pbeta(dea1/(dea1+dea2),dea1,dea2)
    } else if (mono==2) {
       pval = 1-pbeta(a1,dea1,dea2)
       pval = pval/(1-pbeta(dea1/(dea1+dea2),dea1,dea2))
    }
    a2=1-a1
  }
  if (dea1!=1||dea2!=1)
   {
    #search for a2
    if (a1<(dea1/(dea1+dea2)))
      {
        upbd = 1
        lowbd = dea1/(dea1+dea2)
        a2=bisec(delta,dea1,dea2,upbd,lowbd)
      }
    if (a1>=(dea1/(dea1+dea2)))
       {
        upbd = dea1/(dea1+dea2)
        lowbd =0.0000000000000000001
        a2=bisec(delta,dea1,dea2,upbd,lowbd)
        }
      #compute the probability in the beta distribution:
      #P(bigger than max(a1,a2),smaller than min(a1,a2))
    if (mono==0)
          pval = pbeta(min(a1,a2),dea1,dea2)+(1-pbeta(max(a1,a2),dea1,dea2))
    if (mono==1)
        {
          pval = pbeta(a1,dea1,dea2)
          pval = pval/pbeta(dea1/(dea1+dea2),dea1,dea2)
         }
    if (mono==2)
         {
          pval = 1-pbeta(a1,dea1,dea2)
          pval = pval/(1-pbeta(dea1/(dea1+dea2),dea1,dea2))
         }
    }
  returnval=c(a2,pval)
  return(returnval)
}