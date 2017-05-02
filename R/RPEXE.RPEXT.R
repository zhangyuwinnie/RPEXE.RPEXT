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

#' Example data for loopcuts_cuttimes
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name loopcuts_cut
#' @usage data(loopcuts_cut)
NULL

#' Example data for loopcut_times_censoring
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name loopcuts_t_c
#' @usage data(loopcuts_t_c)
NULL

#' Example data for pava
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name pava_dfrd
#' @usage data(pava_dfrd)
NULL

#' Example data for pexeest_times_censoring
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name pexeest_times_censoring
#' @usage data(pexeest_times_censoring)
NULL

#' Example data for pexeest_tx
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name t100
#' @usage data(t100)
NULL

#' Example data for gamllik
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name gamllik_data
#' @usage data(gamllik_data)
NULL

#' Example data for loopcut_onestep
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name loopcut_onestep_data
#' @usage data(loopcut_onestep_data)
NULL

#' Example data for loopcut_umbrella
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name loopcuts_umbrella_cuttimes_mono
#' @usage data(loopcuts_umbrella_cuttimes_mono)
NULL
