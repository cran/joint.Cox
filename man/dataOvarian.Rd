\name{dataOvarian}
\alias{dataOvarian}
\docType{data}
\title{
Survival data of 1003 ovarian cancer patients from 4 independent studies.
}
\description{
The survival data was used to study if the CXCL12 gene expression is a predictive biomarker of survival endpoints in ovarian cancer.
 The dataset was used in Emura et al. (2017), which is a subset of the curated ovarian data of Ganzfried et al (2013).
 We prepared the dataset by using "patientselection.config" in "Curated ovarian data" around May 2015 in the analysis of Emura et al. (2017).
}
\usage{data("dataOvarian")}
\format{
  A data frame with 1003 observations on the following 6 variables.
  \describe{
    \item{\code{t.event}}{: time to event in days}
    \item{\code{event}}{: event indicator (1=recurrence, 0=no recurrence)}
    \item{\code{t.death}}{: time to death in days}
    \item{\code{death}}{: death indicator (1=death, 0=alive)}
    \item{\code{group}}{: study ID; group=4, 8, 11, or 14; see the details below}
    \item{\code{CXCL12}}{: CXCL12 gene expression}
  }
}
\details{
The data include individual-patient information on 1003 patients from 4 studies (group=4, 8, 11, and 14). The numbers 4, 8, 11 and 14 corresponds to the study IDs from the original data of Ganzfried et al. (2013).
"group=4" corresponds to 110 Japanese patients from the study of Yoshihara et al. (2010) (GEO accession number: GSE17260). Other groupds are the studies of GSE30161 (58 patients), GSE9891 (278 patients), and TCGA (557 patients).
}
\source{
Ganzfried BF et al. (2013), Curated ovarian data: clinically annotated data for the ovarian cancer transcriptome, Database, Article ID bat013, doi:10.1093/database/bat013.
}
\references{
Emura T, Nakatochi M, Murotani K, Rondeau V (2017), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Stat Methods Med Res 26(6): 2649-66

Ganzfried BF et al. (2013), Curated ovarian data: clinically annotated data for the ovarian cancer transcriptome, Database, Article ID bat013, doi:10.1093/database/bat013.

Yoshihara K et al. (2010) Gene expression profile for predicting survival in advanced-stage serous ovarian cancer across two independent datasets. PLoS One 12;5(3):e9615
}
\examples{
data(dataOvarian)
study4=dataOvarian[dataOvarian$group==4,] # extract one study
study4
}
\keyword{Datasets}
\keyword{Ovarian cancer}
