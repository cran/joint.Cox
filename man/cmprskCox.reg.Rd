\name{cmprskCox.reg}
\alias{cmprskCox.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Competing Risks Version of Penalized Likelihood Estimation under the Joint Cox Models Between Tumour Progression 
 and Death for Meta-Analysis
}
\description{
Perform regression analyses under a copula-based joint Cox proportional hazards model between tumour progression and death for meta-analysis,
 which is proposed in Section 6.2 of Emura et al. (2017) and Section 5.1 of Emura et al. (2020-).
 This is the competing risks version of "jointCox.reg".
 To avoid the indentifiability problem, the copula parameter (theta) should be given by user, e.g., theta=2.
 The method is applicable for meta-analysis combining several studies or for cluster survival data.
}
\usage{
cmprskCox.reg(t.event, event1, event2, Z1, Z2, group, theta, alpha = 1,
 kappa1 = c(seq(10, 1e+17, length = 30)), kappa2 = c(seq(10, 1e+17, length = 30)),
 LCV.plot = TRUE, Randomize_num = 10, Adj = 500,convergence.par=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.event}{a vector for event tims}
  \item{event1}{a vector for event-type 1 indicators (=1 with event; =0 without event)}
  \item{event2}{a vector for event-type 2 indicators (=1 with event; =0 without event)}
  \item{Z1}{a matrix for covariates associated with event-type 1; ncol(Z1)=the number of covariates}
  \item{Z2}{a matrix for covariates associated with event-type 2; ncol(Z2)=the number of covariates}
  \item{group}{a vector for a group identification number, like 1,2,3....}
  \item{theta}{A copula parameter under the Clayton copula (theta > 0)}
  \item{alpha}{A value related to the frailty (e.g., alpha=0 or =1); alpha=1 is default}
  \item{kappa1}{a vector for candidate smoothing parameters}
  \item{kappa2}{a vector for candidate smoothing parameters}
  \item{LCV.plot}{Plot the LCV curves if "TRUE"}
  \item{Randomize_num}{The number of randomizations for the ititial p0}
  \item{Adj}{Numerical adjustment to prevent overflow; Adj=500 is recommended}
  \item{convergence.par}{If TRUE, the converged estimate, gradient, and Hessian matrix are given (log-transformed)}
}

\details{
We employ "nlm" routine to maximize the penalized likelihood function with the initial value described in Emura et al. (2015).
If "nlm" does not converge, then we randomize the initial value by adding uniform random variables (Hu and Emura, 2015).
}
\value{
\item{count }{Count for event occurences}
\item{beta1 }{Regression coefficient for Z1}
\item{beta2 }{Regression coefficient for Z2}
\item{eta }{Frailty parameter (variance)}
\item{theta }{Copula parameter under the Clayton copula (fixed by user)}
\item{tau }{Kendall's tau corresponding to the copula parameter}
\item{LCV1}{Likelihood cross-validation for event-type 1}
\item{LCV2}{Likelihood cross-validation for event-type 2}
\item{g }{M-spline coefficients for event-type 1}
\item{h }{M-spline coefficients for event-type 2}
\item{g_var }{Variance of M-spline coefficients for event-type 1}
\item{h_var }{Variance of M-spline coefficients for event-type 2}
\item{convergence }{convergence results for maximizing penalized likelihood}
\item{convergence.parameters }{converged estimate, gradient, and Hessian matrix (log-transformed)}
}

\references{
Emura T, Nakatochi M, Murotani K, Rondeau V (2017), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Stat Methods Med Res 26(6): 2649-66

Emura T, Shih JH, Ha ID, Wilke RA (2020-), Comparison of the Marginal Hazard Model and the Sub-distribution Hazard Model for Competing Risks under an Assumed Copula,
 Stat Methods Med Res, in press

Hu YH, Emura T (2015), Maximum likelihood estimation for a special exponential family under random double-truncation,
 Computational Stat 30 (4): 1199-1229
}
\author{Takeshi Emura, Shih JH}

\section{Error }{"Error in integrate(func1, 0.001, 10, stop.on.error = FALSE):non-finite function value",
 an error occurring when the penalized likelihood is maximizad by "nlm". The error may frequently
 occur during the iterations for maximizing the penalized likelihood, but is not crucial (can simply be ignored).}
 
\section{Warning }{"NA/Inf replaced by maximum positive value", an error occurring when the penalized likelihood
 is maximizad by "nlm". The error frequently occurs during the iterations for maximizing the penalized likelihood,
 but is not crucial (can simply be ignored).}

\examples{
data(dataOvarian)
t.event=dataOvarian$t.event
t.death=dataOvarian$t.death
event=dataOvarian$event
death=dataOvarian$death
non.event=which(event==1 & death==1 & t.event==t.death)
non.death=which(event==1 & death==1 & t.event<t.death)
event[non.event]=0 ## relapse before death
death[non.death]=0 ## death before relapse (tie is counted as death)
Z=as.matrix(dataOvarian$CXCL12)
group=dataOvarian$group
alpha_given=0
theta=2.35
kappa_grid=seq(10,1e+17,length = 30)

#set.seed(1)
#cmprskCox.reg(t.event=t.event,event1=event,event2=death,
#              Z1=Z,Z2=Z,group=group,theta=theta,alpha=alpha_given,
#              kappa1=kappa_grid,kappa2=kappa_grid,LCV.plot=TRUE,Adj=500)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Penalized likelihood estimation }
\keyword{ Cox regression }
\keyword{ Competing risk }
\keyword{ Clayton copula }
