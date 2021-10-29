\name{splineCox.reg}
\alias{splineCox.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting the Cox model for survival data using a penalized spline model}
\description{
Fitting the Cox proportional hazards model when the baseline hazard function is specified by a five-parameter spline model.
}
\usage{
splineCox.reg(t.event, event, Z, xi1 = min(t.event), xi3 = max(t.event), 
kappa = c(seq(10, 1e+17, length = 30)), LCV.plot = TRUE,p0=rep(0,5+p))
}
\arguments{
  \item{t.event}{a vector for time-to-event}
  \item{event}{a vector for event indicator (=1 event; =0 censoring)}
  \item{Z}{a matrix for covariates; nrow(Z)=sample size, ncol(Z)=the number of covariates}
  \item{xi1}{lower bound for the hazard function; the default is min(t.event)}
  \item{xi3}{upper bound for the hazard function; the default is max(t.event)}
  \item{kappa}{a vector for candidate smoothing parameters. Only positive values are allowed. Values too close to zero may yeild errors (see below).}
  \item{LCV.plot}{Plot the LCV curves if "TRUE". This plot is used to find the optimal value from the candidate smoothing parameters given by "kappa".}
  \item{p0}{Initial values to maximize the penalized likelihood (5+p parameters; five M-spline coefficients and p regression coefficients)}
}
\details{
One can perform Cox-type regression for censored survival data with covariates. The method is essentially the same as
 as Cox regression (Cox 1972) expect for the models of the baseline hazard function.
 Unlike the nonparametric model of Cox (1972), the method applies a five-parameter spline model as originally proposed by Emura et al. (2017).
 The method is detailed in Section 2.4 of Emura et al. (2019). See also Shih and Emura (2021) for more details.
 This method is also used as a subroutine for computing the optimal smoothing parameter (kappa1 and kappa2) for
 many advanced functions, such as "jointCox.reg", "cmprskCox.reg", and "condCox.reg".
 The definition of LCV is given in Section 3.7 of Emura et al. (2019). See also Shih and Emura (2021).
 The error message "Error in nlm(l.func, p = rep(0, 5 + p), hessian = TRUE):non-finite value supplied by 'nlm'"
 may imply that some candidate parameters for kappa are too close to zero; please exclude such values from kappa.
 The output values are usually similar to those given by "coxph(Surv(t.event,event)~Z)". 
 Unreasonable output values are usually caused by a wrong choice of "kappa" and occasionary caused by a wrong choice of p0.
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
Cox DR (1972), Regression models and life-tables, JRSS(B) 34(2):187-202.

Emura T, Matsui S, Rondeau V (2019), Survival Analysis with Correlated Endpoints; Joint Frailty-Copula Models,
 JSS Research Series in Statistics, Springer

Emura T, Nakatochi M, Murotani K, Rondeau V (2017), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Stat Methods Med Res 26(6): 2649-66: Supplementary Material. 

Shih JH, Emura T (2021) Penalized Cox regression with a five-parameter spline model, Commun Stat-Theor 50(16):3749-68
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

