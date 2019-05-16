\name{splineCox.reg}
\alias{splineCox.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting the Cox model with penalized splines using univariate survival data}
\description{Fitting the Cox model when the baseline hazard function is approximated by splines.}
\usage{
splineCox.reg(t.event, event, Z, xi1 = min(t.event), xi3 = max(t.event), 
kappa = c(seq(10, 1e+17, length = 30)), LCV.plot = TRUE)
}
\arguments{
  \item{t.event}{a vector for time-to-event}
  \item{event}{a vector for event indicator (=1 event; =0 censoring)}
  \item{Z}{a matrix for covariates; nrow(Z)=sample size, ncol(Z)=the number of covariates}
  \item{xi1}{lower bound for the hazard function; the default is min(t.event)}
  \item{xi3}{upper bound for the hazard function; the default is max(t.event)}
  \item{kappa}{a vector for candidate smoothing parameters in LCV. Only positive values are allowed. Values too close to zero may yeild errors (see below).}
  \item{LCV.plot}{Plot the LCV curves if "TRUE"}
}
\details{The definition of LCV is given in Emura et al. (2019). This is a subroutine used to compute the optimal smoothing parameter (kappa) in the penalized likelihood of Emura et al. (2017).
 The error message "Error in nlm(l.func, p = rep(0, 5 + p), hessian = TRUE):non-finite value supplied by 'nlm'"
 may imply that some candidate parameters for kappa are too close to zero; please exclude such values from kappa.
}
\value{
  \item{beta }{Regression coefficient for Z}
  \item{h }{M-spline coefficients}
  \item{h_var }{Variance of M-spline coefficients}
  \item{kappa }{smoothing parameter at the optimal LCV}
  \item{DF }{degree of freedom at the optimal LCV}
  \item{LCV }{the optimal LCV(=logL-DF)}
}
\references{
Emura T, Nakatochi M, Murotani K, Rondeau V (2017), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Stat Methods Med Res 26(6): 2649-66

Emura T, Matsui S, Rondeau V (2019), Survival Analysis with Correlated Endpoint: Joint Frailty-copula models,
 JSS Research Series in Statistics, Springer
}
\author{Takeshi Emura}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(dataOvarian)
t.event=dataOvarian$t.event
event=dataOvarian$event
t.death=dataOvarian$t.death
death=dataOvarian$death
Z=dataOvarian$CXCL12
#splineCox.reg(t.event,event,Z,kappa=c(seq(10,1e+17,length=30)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Spline }
\keyword{ Penalized likelihood estimation }
\keyword{ Cox regression }
\keyword{ LCV }

