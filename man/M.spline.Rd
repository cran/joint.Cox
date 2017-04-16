\name{M.spline}
\alias{M.spline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ M-Spline function }
\description{ Calculate M-Spline bases (5 bases) suggested in Emura et al. (2015). }
\usage{
M.spline(time, xi1, xi3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{ a vector of times }
  \item{xi1}{ lower bound of times }
  \item{xi3}{ upper bound of times }
}

\details{
The "time" argument is a vector satisfying the contraints xi1<=time<=xi3. Otherwise, error messages will be produced.
}
\value{
\item{NULL }{M-Spline bases (5 bases) evaluated at "time".}
}
\references{
Supplementary Material to: Emura T, Nakatochi M, Murotani K, Rondeau V (2015), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Statistical Methods in Medical Research, doi: 10.1177/0962280215604510
}
\author{Takeshi Emura}

\examples{
M.spline(c(1,1.5,2,2.5,3),xi1=1,xi3=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Spline }
\keyword{ M-spline }% __ONLY ONE__ keyword per line
