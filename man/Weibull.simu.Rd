\name{Weibull.simu}
\alias{Weibull.simu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulating data from the Weibull joint frailty-copula model}
\description{
Simulating data from the Weibull joint frailty-copula model.
}
\usage{
Weibull.simu(G,N,scale1,scale2,shape1,shape2,beta1,beta2,eta,theta,alpha,C.max)
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
  \item{theta}{copula parameter}
  \item{alpha}{parameter related to frailty; usually alpha=1}
  \item{C.max}{the upper bound for the censoring distribution}
}

\details{
To be discussed in the paper of Emura et al. (2019).
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
}
\references{
Emura T, Wu BH, Michimae H (2019), The Weibull joint frailty-copula model, in preparation
}
\author{Takeshi Emura}

\examples{
Weibull.simu(G=5,N=2,scale1=1,scale2=1,shape1=1,shape2=1,
        beta1=1,beta2=1,eta=0.5,theta=2,alpha=1,C.max=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Simulation }
\keyword{ Semi-competing risk }
\keyword{ Clayton copula }
\keyword{ Weibull }% __ONLY ONE__ keyword per line
