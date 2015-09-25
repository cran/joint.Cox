\name{F.window}
\alias{F.window}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction of death}
\description{
Dynamic prediction of death using under a joint model
}
\usage{
F.window(time, width, Z1, Z2, beta1, beta2, eta, theta, alpha, g, h, xi1, xi3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{time}
  \item{width}{window}
  \item{Z1}{a vector of covariates}
  \item{Z2}{a vector of covariates; usually Z1=Z2}
  \item{beta1}{regression coefficients for progression}
  \item{beta2}{regression coefficients for death}
  \item{eta}{frailty variance}
  \item{theta}{copula parameter}
  \item{alpha}{parameter related to frailty; usually alpha=1}
  \item{g}{parameters related to the baseline hazard for progression}
  \item{h}{parameters related to the baseline hazard for death}
  \item{xi1}{lower bound for time to event}
  \item{xi3}{upper bound for time to death}
}

\details{
Prediction probability of death is calculated given the event status (X<=t or X>t)
}
\value{
\item{time }{t}
\item{width }{w}
\item{F_event }{F(t+w|t,X<=t,Z1,Z2); usually Z1=Z2}
\item{F_noevent }{F(t+w|t,X>t,Z1,Z2); usually Z1=Z2}
}
\references{
Emura T*, Nakatochi M, Murotani K, Rondeau V (2015), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Statistical Methods in Medical Research, doi: 10.1177/0962280215604510
}
\author{Takeshi Emura}

\examples{
xi1=0.1
xi3=4
###### Both intra-subject and intra-study dependence #######
F.window(time=0.1,width=0.2,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.5,theta=3,alpha=1,
g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)
F.window(time=0.5,width=0.2,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.5,theta=3,alpha=1,
g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)
F.window(time=1.0,width=0.2,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.5,theta=3,alpha=1,
g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)

###### Only intra-study dependence #######
F.window(time=0.1,width=0.2,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.5,theta=0.1,alpha=1,
g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)
F.window(time=0.5,width=0.2,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.5,theta=0.1,alpha=1,
g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)
F.window(time=1.0,width=0.2,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.5,theta=0.1,alpha=1,
g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)

###### Prediction under nearly independence #######
F.window(time=0.1,width=0.2,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.1,theta=0.1,alpha=1,
g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)
F.window(time=0.5,width=0.2,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.1,theta=0.1,alpha=1,
g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)
F.window(time=1.0,width=0.2,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.1,theta=0.1,alpha=1,
g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prediction }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
