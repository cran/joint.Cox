\name{F.windows.Weibull}
\alias{F.windows.Weibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dynamic prediction of death under the joint frailty-copula model (the Weibull baseline hazard functions) }
\description{
Dynamic prediction of death using a joint frailty-copula model. 
Probability of death between t and t+w is calculated given a tumour progression time X and covariates Z1 and Z2. 
If X<=t, the prediction probability is F(t,t+w|X=x, Z1, Z2). If X>t, the prediction probability is F(t,t+w|X>t, Z1, Z2). 
This is a vector version of F.window.Weibull.
}
\usage{
F.windows.Weibull(time, widths, X, Z1, Z2, beta1, beta2, eta, theta, alpha, 
                  scale1, shape1,scale2,shape2, xi1, xi3, Fplot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{prediction time (=t)}
  \item{widths}{length of window (=w)}
  \item{X}{time of tumour progression < time}
  \item{Z1}{a vector of covariates for progression}
  \item{Z2}{a vector of covariates for death}
  \item{beta1}{a vector of regression coefficients for progression}
  \item{beta2}{a vector of regression coefficients for death}
  \item{eta}{frailty variance}
  \item{theta}{copula parameter}
  \item{alpha}{parameter related to frailty; usually alpha=1}
  \item{scale1}{scale parameter related to the baseline hazard for progression}
  \item{shape1}{shape parameter related to the baseline hazard for progression}
  \item{scale2}{scale parameter related to the baseline hazard for death}
  \item{shape2}{shape parameter related to the baseline hazard for death}
  \item{xi1}{lower bound for time to event}
  \item{xi3}{upper bound for time to death}
  \item{Fplot}{if FALSE, the plot is not shown}
}

\details{
Predicted probability of death is calculated given the event status (X<=t or X>t) and covariates (Z1 and Z2).
}
\value{
\item{time }{t}
\item{widths }{w}
\item{X }{X}
\item{F_event_at_X }{F(t,t+w|X=x, Z1, Z2)}
\item{F_noevent }{F(t,t+w|X>t, Z1, Z2)}
}
\references{
Emura T, Nakatochi M, Matsui S, Michimae H, Rondeau V (2018), Personalized dynamic prediction of death according to tumour progression and high-dimensional
 genetic factors: meta-analysis with a joint model, Stat Methods Med Res 27(9):2842-58
}
\author{Sayaka Shinohara, Takeshi Emura}

\examples{
w=c(0,0.5,1,1.5,2)
par(mfrow=c(1,2))
F.windows.Weibull(time=1,X=0.2,widths=w,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.5,theta=8,
                  alpha=1,scale1=1,shape1=1,scale2=1,shape2=1,xi1=0,xi3=3)
F.windows.Weibull(time=1,X=0.8,widths=w,Z1=1,Z2=1,beta1=1,beta2=1,eta=0.5,theta=8,
                  alpha=1,scale1=1,shape1=1,scale2=1,shape2=1,xi1=0,xi3=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Prediction }
\keyword{ Dynamic prediction }% __ONLY ONE__ keyword per line
