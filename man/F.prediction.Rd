\name{F.prediction}
\alias{F.prediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dynamic prediction of death}
\description{
Dynamic prediction of death using a joint frailty-copula model. 
Probability of death between t and t+w is calculated given a tumour progression time X and covariates Z1 and Z2. 
If X<=t, the prediction probability is F(t,t+w|X=x, Z1, Z2). If X>t, the prediction probability is F(t,t+w|X>t, Z1, Z2). 
This function is a simpler version of F.windows. 
}
\usage{
F.prediction(time, widths, X, Z1, Z2, beta1, beta2, eta, theta, alpha,
 g, h, xi1, xi3, Fplot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{prediction time (=t)}
  \item{widths}{length of window (=w)}
  \item{X}{event time (=X) that occurrs before time t; if the event does not occur before time t, one can set arbitrary value X greater than t}
  \item{Z1}{a vector of covariates}
  \item{Z2}{a vector of covariates}
  \item{beta1}{regression coefficients for progression}
  \item{beta2}{regression coefficients for death}
  \item{eta}{frailty variance}
  \item{theta}{copula parameter}
  \item{alpha}{parameter related to frailty; usually alpha=1}
  \item{g}{parameters related to the baseline hazard for progression}
  \item{h}{parameters related to the baseline hazard for death}
  \item{xi1}{lower bound for time-to-event}
  \item{xi3}{upper bound for time-to-death}
  \item{Fplot}{if FALSE, the plot is not shown}
}

\details{
Predicted probability of death is calculated given the event status (X<=t or X>t) and covariates (Z1 and Z2).
}
\value{
\item{time }{t}
\item{widths }{w}
\item{X }{X}
\item{F }{F(t,t+w|X=x, Z1, Z2) or F(t,t+w|X>t, Z1, Z2)}
}
\references{
Emura T, Nakatochi M, Matsui S, Michimae H, Rondeau V (2017), Personalized dynamic prediction of death according to tumour progression and high-dimensional
 genetic factors: meta-analysis with a joint model Statistical Methods in Medical Research, doi:10.1177/0962280216688032
}
\author{Takeshi Emura}

\examples{
w=c(0,0.5,1,1.5,2)
par(mfrow=c(1,2))
F.prediction(time=1,X=0.8,widths=w,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.5,theta=8,
         alpha=1,g=rep(1,5),h=rep(1,5),xi1=0,xi3=3)
F.prediction(time=1,X=1.5,widths=w,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.5,theta=8,
         alpha=1,g=rep(1,5),h=rep(1,5),xi1=0,xi3=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Prediction }
\keyword{ Dynamic prediction }% __ONLY ONE__ keyword per line
