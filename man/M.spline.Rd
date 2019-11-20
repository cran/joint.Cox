\name{M.spline}
\alias{M.spline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ M-spline basis function }
\description{ Calculate the M-spline basis functions (a M-spline basis is a B-spline basis normalized so that the integral is 1).}
\usage{
M.spline(time, xi1, xi3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{ a vector of time points }
  \item{xi1}{ lower bound of time points }
  \item{xi3}{ upper bound of time points }
}

\details{
  The output shows the values of the 5 basis functions at "time",
 giving a matrix with nrow=length(time) and ncol=5.
 The five basis functions were originally given in the Supplementary Material of Emura et al. (2017).
 More details can be found in Emura and Chen (2018), Emura et al. (2019), and Shih and Emura (2020-). 
 The "time" argument should be a vector satisfying the contraints xi1<=time<=xi3.
 If "time" does not meet the constraints, error messages are shown.
}
\value{
\item{NULL }{A matrix with nrow=length(time) and ncol=5, containing the values of the 5 spline basis functions at "time".}
}
\references{
Emura T, Chen YH (2018). Analysis of Survival Data with Dependent Censoring, 
Copula-Based Approaches, JSS Research Series in Statistics, Springer, Singapore.

Emura T, Matsui S, Rondeau V (2019), Survival Analysis with Correlated Endpoints; Joint Frailty-Copula Models,
 JSS Research Series in Statistics, Springer

Emura T, Nakatochi M, Murotani K, Rondeau V (2017), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Stat Methods Med Res 26(6): 2649-66: Supplementary Material. 

Shih JH, Emura T (2020-), Penalized Cox regression with a five-parameter spline model, Commun Stat-Theor, in revision.
}
\author{Takeshi Emura}

\examples{
M.spline(time=c(1,2,3),xi1=1,xi3=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Spline }
\keyword{ M-spline }% __ONLY ONE__ keyword per line
