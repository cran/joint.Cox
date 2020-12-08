\name{condCox.reg}
\alias{condCox.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Penalized Likelihood Estimation under the Joint Cox Models Between Tumour Progression 
 and Death for Meta-Analysis; A Conditional Copula Approach
}
\description{
Perform regression analyses under a copula-based joint Cox proportional hazards model between tumour progression and death for meta-analysis,
 which is proposed by Emura et al. (201x).
 The method extends the joint frailty copula model of Emura et al. (2017).
 The method is applicable for meta-analysis combining several studies or for cluster survival data.
}
\usage{
condCox.reg(t.event, event, t.death, death, Z1, Z2, Z12, group, alpha = 1,
 kappa1 = c(seq(10, 1e+17, length = 30)),kappa2 = c(seq(10, 1e+17, length = 30)),
 LCV.plot = TRUE, Randomize_num = 10, u.min = 0.001, u.max = 10,
 Adj = 500,convergence.par=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.event}{a vector for time-to-tumour progression (TTP)}
  \item{event}{a vector for progression indicator (=1 if progression; =0 if not progression)}
  \item{t.death}{a vector for overall survival (OS), i.e., time-to-death}
  \item{death}{a vector for death indicator(=1 if death; =0 if not death)}
  \item{Z1}{a matrix for covariates associated with TTP; ncol(Z1)=the number of covariates}
  \item{Z2}{a matrix for covariates associated with OS; ncol(Z2)=the number of covariates}
  \item{Z12}{a matrix for covariates associated with copula; ncol(Z12)=the number of covariates}
  \item{group}{a vector for group identification numbers, like 1,2,3....}
  \item{alpha}{A value related to the frailty (e.g., alpha=0 or =1); alpha=1 is default}
  \item{kappa1}{a vector for candidate smoothing parameters}
  \item{kappa2}{a vector for candidate smoothing parameters}
  \item{LCV.plot}{Plot the LCV curves if "TRUE"}
  \item{Randomize_num}{The number of randomizations for the ititial p0}
  \item{u.min}{the lower bound of the numerical integration for the frailty term}
  \item{u.max}{the upper bound of the numerical integration for the frailty term}
  \item{Adj}{Numerical adjustment to prevent overflow; Adj=500 is recommended}
  \item{convergence.par}{If TRUE, the converged estimate, gradient, and Hessian matrix are given (log-transformed)}
}

\details{
We employ "nlm" routine to maximize the penalized likelihood function with the initial value described in Emura et al. (2017).
If "nlm" does not converge, then we randomize the initial value by adding uniform random variables (Hu and Emura, 2015).
}
\value{
\item{count }{Count for event occurences}
\item{beta1 }{Regression coefficient for Z1}
\item{beta2 }{Regression coefficient for Z2}
\item{eta }{Frailty parameter (variance)}
\item{theta }{Baseline copula parameter under the Clayton copula}
\item{tau }{Kendall's tau corresponding to the baseline copula parameter}
\item{beta12 }{Regression coefficient for a copula parameter}
\item{LCV1}{Likelihood cross-validation for TTP}
\item{LCV2}{Likelihood cross-validation for OS}
\item{g }{M-spline coefficients for TTP}
\item{h }{M-spline coefficients for OS}
\item{g_var }{Variance of M-spline coefficients for TTP}
\item{h_var }{Variance of M-spline coefficients for OS}
\item{convergence }{convergence results for maximizing penalized likelihood}
\item{convergence.parameters }{converged estimate, gradient, and Hessian matrix (log-transformed)}
}

\references{
Emura T, Rondeau V (201x), Kendalls tau for individual-level surrogacy for failure time endpoints in meta-analysis (in preparation)

Emura T, Nakatochi M, Murotani K, Rondeau V (2017), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Stat Methods Med Res 26(6): 2649-66

Hu YH, Emura T (2015), Maximum likelihood estimation for a special exponential family under random double-truncation,
 Computational Stat 30 (4): 1199-1229
}
\author{Takeshi Emura}

\section{Error }{"Error in integrate(func1, 0.001, 10, stop.on.error = FALSE):non-finite function value",
 an error occurring when the penalized likelihood is maximizad by "nlm". The error may frequently
 occur during the iterations for maximizing the penalized likelihood, but is not crucial (can simply be ignored).}
 
\section{Warning }{"NA/Inf replaced by maximum positive value", an error occurring when the penalized likelihood
 is maximizad by "nlm". The error frequently occurs during the iterations for maximizing the penalized likelihood,
 but is not crucial (can simply be ignored).}

\examples{
##
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Penalized likelihood estimation }
\keyword{ Semi-competing risk }
\keyword{ Cox regression }
\keyword{ Clayton copula }
