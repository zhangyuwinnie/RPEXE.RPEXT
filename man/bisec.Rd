% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bisec.R
\name{bisec}
\alias{bisec}
\title{Bisection algorithm in Beta distribution}
\usage{
bisec(delta, dea1, dea2, upbd, lowbd)
}
\arguments{
\item{delta}{Test statistic in Han et al. (2012), 
delta = (ttot1/(ttot1+ttot2))^dea1*(ttot2/(ttot1+ttot2))^dea2;}

\item{dea1}{first parameter in Beta distribution (number of events from the first arm)}

\item{dea2}{second parameter in Beta distribution (number of events from the second arm)}

\item{upbd}{upper bound of a2}

\item{lowbd}{lower bound of a2}
}
\value{
a2
}
\description{
Running bisection algorithm to search for a2, the minimizer of (log((a2)^dea1*(1-a2)^dea2-delta))^2
}
\examples{
bisec(-74.4824, 33, 98, 1, 0.252)
}
