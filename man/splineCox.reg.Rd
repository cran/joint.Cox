\name{splineCox.reg}
\alias{splineCox.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Likelihood cross-validation (LCV) for the Cox model with penalized splines}
\description{Calculate likelihood cross-validation (LCV) values for the Cox model when the baseline hazard is
modeled via penalized splines.}
\usage{
splineCox.reg(t.event, event, Z, xi1 = min(t.event), xi3 = max(t.event), 
kappa_grid = c(seq(10, 1e+17, length = 30)), LCV_plot = TRUE)
}
\arguments{
  \item{t.event}{a vector object for time-to-event}
  \item{event}{a vector object for event indicator (=1 event; =0 censoring)}
  \item{Z}{a matrix object for covariates; nrow(Z)=sample size, ncol(Z)=the number of covariates}
  \item{xi1}{lower bound for the hazard function; the default is min(t.event)}
  \item{xi3}{upper bound for the hazard function; the default is max(t.event)}
  \item{kappa_grid}{a vector for candidate smoothing parameters in LCV. Only positive values are allowed. Values too close to zero may yeild errors (see below).}
  \item{LCV_plot}{Plot the LCV curves if "TRUE"}
}
\details{The definition of LCV is given in Emura et al. (2018). This is a subroutine used to compute the optimal smoothing parameter (kappa) in the penalized likelihood of Emura et al. (2015).
 The error message "Error in nlm(l.func, p = rep(0, 5 + p), hessian = TRUE):non-finite value supplied by 'nlm'"
 may imply that some candidate parameters for kappa_grid are too close to zero; please exclude such values from kappa_grid.
}
\value{
  \item{kappa }{smoothing parameter at the optimal LCV}
  \item{DF }{degree of freedom at the optimal LCV}
  \item{LCV }{the optimal LCV(=logL-DF)}
}
\references{
Emura T, Nakatochi M, Murotani K, Rondeau V (2015), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Statistical Methods in Medical Research, doi: 10.1177/0962280215604510

Emura T, Chen YH, Matsui S, and Rondeau V (2018), Survival Analysis with Dependent Censoring and Correlated Endpoints: Copula-Based Approaches,
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
Z=as.matrix(dataOvarian$CXCL12)
#splineCox.reg(t.event,event,Z,kappa_grid=c(seq(10,1e+17,length=30)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Spline }
\keyword{ Penalized likelihood estimation }
\keyword{ Cox regression }
\keyword{ LCV }

