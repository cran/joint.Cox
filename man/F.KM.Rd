\name{F.KM}
\alias{F.KM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction of death using the Kaplan-Meier estimator}
\description{
Dynamic prediction of death using using the Kaplan-Meier estimator. 
Probability of death between t and t+w is calculated. The prediction probability is F(t,t+w)=1-S(t+w)/S(t), where S is the Kaplan-Meier estimator.
}
\usage{
F.KM(time, widths, t.death, death)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{prediction time (=t)}
  \item{widths}{length of window (=w)}
  \item{t.death}{a vector object for overall survival (OS), i.e., time-to-death}
  \item{death}{a vector object for death indicator(=1 if death; =0 if not death)}
}

\details{
Prediction probability of death is calculated without covariates.
}
\value{
\item{time }{t}
\item{widths }{w}
\item{F }{F(t,t+w)}
}
\references{
Emura T, Nakatochi M, Matsui S, Michimae H, Rondeau V (2017), Personalized dynamic prediction of death according to tumour progression and high-dimensional
 genetic factors: meta-analysis with a joint model Statistical Methods in Medical Research, doi:10.1177/0962280216688032
}
\author{Takeshi Emura}

\examples{

time=1
widths=c(0,0.5,1,1.5,2)
t.death=c(0.5,1,1.5,2,2.5,3)
death=c(1,1,1,1,1,1)
F.KM(time=time,width=widths,t.death=t.death,death=death)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Prediction }
\keyword{ Dynamic prediction }% __ONLY ONE__ keyword per line
