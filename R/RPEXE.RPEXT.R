#' @title RPEXE.RPEXT
#' @name RPEXE.RPEXT
#' @docType package
#' @description The package allows you to .....
#'
#' @details This reduced piecewise exponential survival software implements the likelihood 
#' ratio test procedure in Han, Schell, and Kim (2009). Inputs to the program can be either 
#' times when events/censoring occur or the vectors of total time on test and the number of 
#' events. Outputs of the programs are times of events and the corresponding p- values. The 
#' order of times and p-values is determined by a backward elimination procedure. Details 
#' about the model and implementation are given in Han, Schell, and Kim (2009). This program
#' can run in R version 12 and above.
#' 
#' @references 
#' Han, G., Schell, M. J., and Kim, J. (2009) Improved survival modeling using a piece- wise exponential approach
#' 
NULL

#' RPEXE_fitting
#' 
#' A dataset containing ....
#' 
#' \itemize{
#'   \item first column. times
#'   \item second column. censor
#'   \item group
#' }
#' @docType data
#' @keywords datasets
#' @name data2
#' @usage data(data2)
NULL

#' JAMA Breast cancer 
#' 
#' A dataset containing ....
#' 
#' \itemize{
#'   \item validate: ... 
#'   \item drfs: .... 
#'   \item drfs.time: ...
#'   \item group
#' }
#' @docType data
#' @keywords datasets
#' @name df
#' @usage data(df)
NULL

#' None Small Cell Lung cancer data
#' 
#' A dataset containing ....
#' 
#' \itemize{
#'   \item validate: ... 
#'   \item drfs: .... 
#'   \item drfs.time: ...
#'   \item group
#' }
#' @docType data
#' @keywords datasets
#' @name simple
#' @usage data(simple)
NULL
