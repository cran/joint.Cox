\name{F.window}
\alias{F.window}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction of death}
\description{
Dynamic prediction of death using under a joint model
}
\usage{
F.window(time, width, X, Z1, Z2, beta1, beta2, eta, theta, alpha,
 g, h, xi1, xi3, Fplot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{prediction time (=t)}
  \item{width}{length of window (=w)}
  \item{X}{event time occurred < time}
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
  \item{Fplot}{if FALSE, the plot is not shown}
}

\details{
Prediction probability of death is calculated given the event status (X<=t or X>t)
}
\value{
\item{time }{t}
\item{width }{w}
\item{X }{X}
\item{F_event_at_X }{F(t+w|t,X=x,Z1,Z2); usually Z1=Z2; must be x<t}
\item{F_noevent }{F(t+w|t,X>t,Z1,Z2); usually Z1=Z2}
}
\references{
Emura T*, Nakatochi M, Murotani K, Rondeau V (2015), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Statistical Methods in Medical Research, doi: 10.1177/0962280215604510
}
\author{Takeshi Emura}

\examples{

xi1=0
xi3=3
time=1
w=1
### Both intra-subject and intra-cluster dependence ######
par(mfrow=c(1,2))
F.window(time=time,X=0.2,width=w,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.5,theta=8,
         alpha=1,g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)
F.window(time=time,X=0.8,width=w,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.5,theta=8,
         alpha=1,g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)

### Only intra-cluster dependence ######
F.window(time=time,X=0.2,width=w,Z1=1,Z2=1,beta1=1,beta2=1,eta=1,theta=0,
         alpha=1,g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)
F.window(time=time,X=0.8,width=w,Z1=1,Z2=1,beta1=1,beta2=1,eta=1,theta=0,
         alpha=1,g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)

### Only intra-subject dependence ######
F.window(time=time,X=0.2,width=w,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.05,theta=8,
         alpha=1,g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)
F.window(time=time,X=0.8,width=w,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.05,theta=8,
         alpha=1,g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)

### Almost no dependence ######
F.window(time=time,X=0.2,width=w,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.005,theta=0,
         alpha=1,g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)
F.window(time=time,X=0.8,width=w,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.005,theta=0,
         alpha=1,g=rep(1,5),h=rep(1,5),xi1=xi1,xi3=xi3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prediction }
\keyword{ dynamic prediction }% __ONLY ONE__ keyword per line
