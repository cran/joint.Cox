\name{dataOvarian}
\alias{dataOvarian}
\docType{data}
\title{
Meta-analytic data of ovarian cancer patients combining 4 independent studies.
}
\description{
Meta-analytic data for studying the CXCL12 gene expression as a predictive
 biomarker of survival in ovarian cancer.
 The dataset is a subset of the curated ovarian data of Ganzfried et al (2013).
 We prepared the dataset by using "patientselection.config" in "Curated ovarian data" around May 2015. 
}
\usage{data("dataOvarian")}
\format{
  A data frame with 1003 observations on the following 6 variables.
  \describe{
    \item{\code{t.event}}{: time to event in days}
    \item{\code{event}}{: event indicator (1=recurrence, 0=no recurrence)}
    \item{\code{t.death}}{: time to death in days}
    \item{\code{death}}{: death indicator (1=death, 0=alive)}
    \item{\code{group}}{: study ID; group=4, 8, 11, or 14}
    \item{\code{CXCL12}}{: CXCL12 expression}
  }
}
\details{
4 studies are combined (group=4, 8, 11, and 14). The numbers 4, 8, 11 and 14 corresponds to the IDs 
from the original data of Ganzfried et al. (2013).}
\source{
Ganzfried BF et al. (2013), Curated ovarian data: clinically annotated data for the ovarian cancer transcriptome, Database, Article ID bat013, doi:10.1093/database/bat013.
}
\references{
Ganzfried BF et al. (2013), Curated ovarian data: clinically annotated data for the ovarian cancer transcriptome, Database, Article ID bat013, doi:10.1093/database/bat013.
}
\examples{
data(dataOvarian)
study4=dataOvarian[dataOvarian$group==4,] # extract one study
study4
}
\keyword{Datasets}
\keyword{Ovarian cancer}