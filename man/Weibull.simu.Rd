\name{Weibull.simu}
\alias{Weibull.simu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulating data from the Weibull joint frailty-copula model}
\description{
Simulating data from the Weibull joint frailty-copula model.
 The form of copulas can be the Clayton (default), Frank, Gumbel, or BB1.
}
\usage{
Weibull.simu(G,N,scale1,scale2,shape1,shape2,beta1,beta2,
 eta,copula="Clayton",theta,d=0,alpha,beta12=0,C.max,tau=FALSE,Z.dist=runif,...)
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
  \item{copula}{copula function; either "Clayton" (default) or "Gumbel"}
  \item{theta}{copula parameter}
  \item{d}{BB1 copula's departure parameter from the Clayton (d=0 is default)}
  \item{alpha}{parameter related to frailty, e.g., alpha=1}
  \item{beta12}{regression coefficients for copula}
  \item{C.max}{the upper bound for the censoring distribution}
  \item{tau}{if TRUE, conditional Kendall's tau given Z is shown}
  \item{Z.dist}{the distribution of a covariate Z}
  \item{...}{parameters for Z.dist}
}

\details{
See Wu et al. (2020) for the algorithms for the Clayton copula. The model has been extended by including the covariate effect
 on a copula (beta12). See the conditional copula model of Emura et al. (2021).
 The Frank, Gumbel, and BBB1 copulas are also available.
}
\value{
\item{X}{: time to event}
\item{D}{: time to death}
\item{C}{: independent censoring time}
\item{t.event}{: time to event (censored)}
\item{event}{: event indicator (1=event, 0=no event)}
\item{t.death}{: time to death (censored)}
\item{death}{: death indicator (1=death, 0=alive)}
\item{group}{: study ID (1~G)}
\item{Z}{: covariate}
\item{tau}{: Kendall's tau}
}
\references{
Wu BH, Michimae H, Emura T (2020), Meta-analysis of individual patient data with semi-competing risks under
 the Weibull joint frailty-copula model. Comp Stat 35(4):1525-52

Emura T, Sofeu C, Rondeau V (2021), Conditional copula models for correlated survival endpoints in meta-analyses of randomized controlled trials, in revision
}
\author{Takeshi Emura}

\examples{
Weibull.simu(G=5,N=2,scale1=1,scale2=1,shape1=1,shape2=1,
        beta1=1,beta2=1,eta=0.5,theta=2,alpha=1,C.max=5)

Weibull.simu(G=5,N=2,scale1=1,scale2=1,shape1=1,shape2=1,
        beta1=1,beta2=1,eta=0.5,copula="Gumbel",theta=2,alpha=1,C.max=5)

Weibull.simu(G=5,N=2,scale1=1,scale2=1,shape1=1,shape2=1,
        beta1=1,beta2=1,eta=0.5,theta=2,alpha=1,C.max=5,Z.dist=rbinom,size=1,prob=0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Simulation }
\keyword{ Semi-competing risk }
\keyword{ Clayton copula }
\keyword{ Weibull }% __ONLY ONE__ keyword per line
