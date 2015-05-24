\name{dataOvarian}
\alias{dataOvarian}
\docType{data}
\title{
Meta-analysis data of ovarian cancer patients
}
\description{
Meta-analysis data for studying the CXCL12 gene expression as a predictive
 biomarker of survival in ovarian cancer.
 This is a subset of the curated ovarian data of Ganzfried et al (2013). 
}
\usage{data("dataOvarian")}
\format{
  A data frame with 1003 observations on the following 6 variables.
  \describe{
    \item{\code{t.event}}{: time to event}
    \item{\code{event}}{: event indicator}
    \item{\code{t.death}}{: time to death}
    \item{\code{death}}{: death indicator}
    \item{\code{group}}{: study ID; groud=4, 8, 11, or 14}
    \item{\code{CXCL12}}{: CXCL12 expression}
  }
}
\details{
4 studies are combined (groud=4, 8, 11, and 14).
}
\source{
Ganzfried BF et al. (2013), Curated ovarian data: clinically annotated data for the ovarian cancer transcriptome, Database, Article ID bat013, doi:10.1093/database/bat013.
}
\references{
Ganzfried BF et al. (2013), Curated ovarian data: clinically annotated data for the ovarian cancer transcriptome, Database, Article ID bat013, doi:10.1093/database/bat013.
}
\examples{
data(dataOvarian)
GSE17260=dataOvarian[dataOvarian$group==4,] # extract the study GSE17260 
GSE17260
}
\keyword{datasets}
