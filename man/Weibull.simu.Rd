\name{Weibull.simu}
\alias{Weibull.simu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulating data from the Weibull joint frailty-copula model}
\description{
This function generate clustered (grouped) bivariate event times from the joint frailty-copula model
 with the Weibull baseline hazard functions.
 Simulating (X_ij,D_ij,C_ij), i=1,2,...,G, and j=1,2,...,N, where G is the number of studies (groups),
 and N is the number of individuals (patients) within each study.
 X_ij is time-to-event, D_ij is time-to-death, and C_ij is time-to-censoring.
 (X_ij, D_ij) and C_ij are independent.
 Dependence structure on (X_ij, D_ij) is modeled by a copula, which can be the Clayton (default), Frank, Gumbel, or BB1.
 Covariate effects are specified by the Cox models given a frailty term. 
}
\usage{
Weibull.simu(G,N,scale1,scale2,shape1,shape2,beta1,beta2,
 eta,copula="Clayton",theta,d=0,alpha,beta12=0,C.max,
 cmprsk=FALSE,tau=FALSE,Z.dist=runif,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{The number of studies or groups}
  \item{N}{The number of patients within each study}
  \item{scale1}{scale parameter related to the baseline hazard for progression}
  \item{scale2}{scale parameter related to the baseline hazard for death}
  \item{shape1}{shape parameter related to the baseline hazard for progression}
  \item{shape2}{shape parameter related to the baseline hazard for death}
  \item{beta1}{regression coefficients for progression}
  \item{beta2}{regression coefficients for death}
  \item{eta}{frailty variance}
  \item{copula}{copula function; "Clayton" (default), "Gumbel", "Frank", or "BB1"}
  \item{theta}{copula parameter}
  \item{d}{BB1 copula's departure parameter from the Clayton (d=0 is the default)}
  \item{alpha}{parameter related to frailty, e.g., alpha=1}
  \item{beta12}{regression coefficients for copula}
  \item{C.max}{the upper bound for the censoring distribution}
  \item{cmprsk}{if TRUE, simulated data follow the competing risks setting}
  \item{tau}{if TRUE, conditional Kendall's tau given Z is shown}
  \item{Z.dist}{the distribution of a covariate Z}
  \item{...}{parameters for Z.dist}
}

\details{
See Wu et al. (2020) for the algorithms for the Clayton copula.
 The method was later extended by including covariate effects
 on a copula (beta12) via the conditional copula model of Emura et al. (2021),
 The available copulas are the Frank, Gumbel, and BB1 copulas.
 For the BB1 copula, please see Supplementary Material:Additional simulation studies under
 the copula misspecification in Emura et al. (2021),
}
\value{
\item{X}{: time to event}
\item{D}{: time to death}
\item{C}{: time to independent censoring}
\item{t.event}{: time to event (=min(X,D,C))}
\item{event}{: event indicator (=I(X<=D,X<=C))}
\item{event1}{: indicator for Event 1 (=I(X<=D,X<=C))}
\item{t.death}{: time to death (=min(D,C))}
\item{death}{: death indicator (=I(D<=C))}
\item{event2}{: indicator for Event 2 (=I(D<X,D<=C))}
\item{group}{: study ID (=1,2,...,G)}
\item{Z}{: covariate}
\item{tau}{: Conditional Kendall's tau given Z}
}
\references{
Wu BH, Michimae H, Emura T (2020), Meta-analysis of individual patient data with semi-competing risks under
 the Weibull joint frailty-copula model. Comp Stat 35(4):1525-52

Emura T, Shih JH, Ha ID, Wilke RA (2020), Comparison of the marginal hazard model and the sub-distribution hazard model for competing risks under an assumed copula,
 Stat Methods Med Res 29(8):2307-27

Emura T, Sofeu C, Rondeau V (2021), Conditional copula models for correlated survival endpoints:
 individual patient data meta-analysis of randomized controlled trials,
 Stat Methods Med Res 30(12):2634-50

Supplementary Material:Additional simulation studies under the copula misspecification
 in "Emura T, Sofeu C, Rondeau V (2021), Conditional copula models for correlated survival endpoints:
 individual patient data meta-analysis of randomized controlled trials,
 Stat Methods Med Res 30(12):2634-50"
}
\author{Takeshi Emura}

\examples{
Weibull.simu(G=5,N=2,scale1=1,scale2=1,shape1=1,shape2=1,
        beta1=1,beta2=1,eta=0.5,theta=2,alpha=1,C.max=5)

Weibull.simu(G=5,N=2,scale1=1,scale2=1,shape1=1,shape2=1,
        beta1=1,beta2=1,eta=0.5,copula="Gumbel",theta=2,alpha=1,C.max=5)

Weibull.simu(G=5,N=2,scale1=1,scale2=1,shape1=1,shape2=1,
        beta1=1,beta2=1,eta=0.5,theta=2,alpha=1,C.max=5,Z.dist=rbinom,size=1,prob=0.5)

## simulated data follow the competing risks setting
Weibull.simu(G=5,N=2,scale1=1,scale2=1,shape1=1,shape2=1,
        beta1=1,beta2=1,eta=0.5,theta=2,alpha=1,C.max=5,cmprsk=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Simulation }
\keyword{ Semi-competing risk }
\keyword{ Competing risk }
\keyword{ Clayton copula }
\keyword{ Gumbel copula }
\keyword{ Frank copula }
\keyword{ BB1 copula }
\keyword{ Weibull }% __ONLY ONE__ keyword per line
