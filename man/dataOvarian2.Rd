\name{dataOvarian2}
\alias{dataOvarian2}
\docType{data}
\title{
Data on time-to-death and 128 gene expressions for 912 ovarian cancer patients from 4 independent studies.
}
\description{
Meta-analytic data containing 128 gene expressions and time-to-death information for ovarian cancer patients.
 The data include time-to-death, residual tumour size (>=1cm> vs. <1cm), and associated 128 gene expressions.
 The dataset is a subset of the curated ovarian data of Ganzfried et al (2013).
 We prepared the dataset by using "patientselection.config" in "Curated ovarian data" around October 2016. 
}
\usage{data("dataOvarian2")}
\format{
  A data frame with 912 observations on the following 132 variables.
  \describe{
    \item{\code{t.death}}{: time to death in days}
    \item{\code{death}}{: death indicator (1=death, 0=alive)}
    \item{\code{group}}{: study ID; group=4, 9, 12, or 16}
    \item{\code{debulk}}{: residual tumour size (>=1cm> vs. <1cm)}
    \item{\code{ANKRD27}}{a numeric vector}
    \item{\code{AP3S1}}{a numeric vector}
    \item{\code{APMAP}}{a numeric vector}
    \item{\code{ARHGAP28}}{a numeric vector}
    \item{\code{ASAP1}}{a numeric vector}
    \item{\code{ASAP3}}{a numeric vector}
    \item{\code{ASB7}}{a numeric vector}
    \item{\code{B4GALT5}}{a numeric vector}
    \item{\code{BYSL}}{a numeric vector}
    \item{\code{C1QTNF3}}{a numeric vector}
    \item{\code{CASP8}}{a numeric vector}
    \item{\code{CCL18}}{a numeric vector}
    \item{\code{CD79A}}{a numeric vector}
    \item{\code{CDK19}}{a numeric vector}
    \item{\code{CLIC4}}{a numeric vector}
    \item{\code{COL11A1}}{a numeric vector}
    \item{\code{COL16A1}}{a numeric vector}
    \item{\code{COL3A1}}{a numeric vector}
    \item{\code{COL5A1}}{a numeric vector}
    \item{\code{COL5A2}}{a numeric vector}
    \item{\code{COMP}}{a numeric vector}
    \item{\code{COX7A2P2}}{a numeric vector}
    \item{\code{CPNE1}}{a numeric vector}
    \item{\code{CRISPLD2}}{a numeric vector}
    \item{\code{CRYAB}}{a numeric vector}
    \item{\code{CTNNBL1}}{a numeric vector}
    \item{\code{CXCL12}}{a numeric vector}
    \item{\code{CXCL9}}{a numeric vector}
    \item{\code{CYBRD1}}{a numeric vector}
    \item{\code{CYR61}}{a numeric vector}
    \item{\code{CYTH3}}{a numeric vector}
    \item{\code{DDX27}}{a numeric vector}
    \item{\code{DLGAP4}}{a numeric vector}
    \item{\code{DNAJC13}}{a numeric vector}
    \item{\code{DYNLRB1}}{a numeric vector}
    \item{\code{EFNB2}}{a numeric vector}
    \item{\code{EIF3K}}{a numeric vector}
    \item{\code{ELN}}{a numeric vector}
    \item{\code{EMP1}}{a numeric vector}
    \item{\code{ENPP1}}{a numeric vector}
    \item{\code{FABP4}}{a numeric vector}
    \item{\code{FAP}}{a numeric vector}
    \item{\code{FBL}}{a numeric vector}
    \item{\code{FGF1}}{a numeric vector}
    \item{\code{FOXN3}}{a numeric vector}
    \item{\code{FSTL1}}{a numeric vector}
    \item{\code{GABRG3}}{a numeric vector}
    \item{\code{GAS1}}{a numeric vector}
    \item{\code{GFRA1}}{a numeric vector}
    \item{\code{GJC1}}{a numeric vector}
    \item{\code{GPATCH1}}{a numeric vector}
    \item{\code{GZMB}}{a numeric vector}
    \item{\code{HLA.DOB}}{a numeric vector}
    \item{\code{HOXA5}}{a numeric vector}
    \item{\code{HP1BP3}}{a numeric vector}
    \item{\code{HSD17B6}}{a numeric vector}
    \item{\code{IL2RG}}{a numeric vector}
    \item{\code{INHBA}}{a numeric vector}
    \item{\code{ITGB1}}{a numeric vector}
    \item{\code{ITPKC}}{a numeric vector}
    \item{\code{JAM2}}{a numeric vector}
    \item{\code{JUN}}{a numeric vector}
    \item{\code{KCNH4}}{a numeric vector}
    \item{\code{KDELC1}}{a numeric vector}
    \item{\code{KIAA0355}}{a numeric vector}
    \item{\code{KIN}}{a numeric vector}
    \item{\code{LEP}}{a numeric vector}
    \item{\code{LOX}}{a numeric vector}
    \item{\code{LPL}}{a numeric vector}
    \item{\code{LSM14A}}{a numeric vector}
    \item{\code{LUM}}{a numeric vector}
    \item{\code{LUZP1}}{a numeric vector}
    \item{\code{MAPRE1}}{a numeric vector}
    \item{\code{MCL1}}{a numeric vector}
    \item{\code{MEOX2}}{a numeric vector}
    \item{\code{MMP12}}{a numeric vector}
    \item{\code{N4BP2L2}}{a numeric vector}
    \item{\code{NCOA3}}{a numeric vector}
    \item{\code{NCOA6}}{a numeric vector}
    \item{\code{NOTCH2NL}}{a numeric vector}
    \item{\code{NR1H3}}{a numeric vector}
    \item{\code{NUAK1}}{a numeric vector}
    \item{\code{OAT}}{a numeric vector}
    \item{\code{OMD}}{a numeric vector}
    \item{\code{PAK4}}{a numeric vector}
    \item{\code{PCDH9}}{a numeric vector}
    \item{\code{PDP1}}{a numeric vector}
    \item{\code{PDPN}}{a numeric vector}
    \item{\code{PHF20}}{a numeric vector}
    \item{\code{PLXNA1}}{a numeric vector}
    \item{\code{PSMC4}}{a numeric vector}
    \item{\code{PSMD8}}{a numeric vector}
    \item{\code{RAB13}}{a numeric vector}
    \item{\code{RAI14}}{a numeric vector}
    \item{\code{RARRES1}}{a numeric vector}
    \item{\code{RBM39}}{a numeric vector}
    \item{\code{RECQL}}{a numeric vector}
    \item{\code{RIN2}}{a numeric vector}
    \item{\code{RND3}}{a numeric vector}
    \item{\code{RPS16}}{a numeric vector}
    \item{\code{SACS}}{a numeric vector}
    \item{\code{SH3PXD2A}}{a numeric vector}
    \item{\code{SKI}}{a numeric vector}
    \item{\code{SLAMF7}}{a numeric vector}
    \item{\code{SLC37A4}}{a numeric vector}
    \item{\code{SMG5}}{a numeric vector}
    \item{\code{SOCS5}}{a numeric vector}
    \item{\code{SPARC}}{a numeric vector}
    \item{\code{SSR4}}{a numeric vector}
    \item{\code{STAU1}}{a numeric vector}
    \item{\code{SUPT5H}}{a numeric vector}
    \item{\code{TBCB}}{a numeric vector}
    \item{\code{TBCC}}{a numeric vector}
    \item{\code{TEAD1}}{a numeric vector}
    \item{\code{TESK1}}{a numeric vector}
    \item{\code{TIMP3}}{a numeric vector}
    \item{\code{TJP1}}{a numeric vector}
    \item{\code{TP53BP2}}{a numeric vector}
    \item{\code{TSPAN9}}{a numeric vector}
    \item{\code{TTI1}}{a numeric vector}
    \item{\code{TUBB2A}}{a numeric vector}
    \item{\code{TUBB6}}{a numeric vector}
    \item{\code{URI1}}{a numeric vector}
    \item{\code{USP48}}{a numeric vector}
    \item{\code{YWHAB}}{a numeric vector}
    \item{\code{ZFP36}}{a numeric vector}
    \item{\code{ZFP36L2}}{a numeric vector}
    \item{\code{ZNF148}}{a numeric vector}
  }
}
\details{
4 studies are combined (group=4, 9, 12, and 16). The numbers 4, 9, 12 and 16 corresponds to the IDs 
from the original data of Ganzfried et al. (2013).}
\source{
Ganzfried BF et al. (2013), Curated ovarian data: clinically annotated data for the ovarian cancer transcriptome, Database, Article ID bat013, doi:10.1093/database/bat013.
}
\references{
Ganzfried BF et al. (2013), Curated ovarian data: clinically annotated data for the ovarian cancer transcriptome, Database, Article ID bat013, doi:10.1093/database/bat013.
}
\examples{
data(dataOvarian2)
######## univariate Cox ##########
t.death=dataOvarian2$t.death
death=dataOvarian2$death
X.mat=dataOvarian2[,-c(1,2,3,4)] ## gene expression
Symbol=colnames(dataOvarian2)[-c(1,2,3,4)] ## gene symbol

p=ncol(X.mat)
P_value=coef=NULL
for(j in 1:p){
  res=summary(coxph(Surv(t.death,death)~X.mat[,j]))$coefficients
  P_value=c(P_value,res[5])
  coef=c(coef,res[1])
}
data.frame( gene=Symbol[order(P_value)], P=P_value[order(P_value)], 
coef=round(coef[order(P_value)],3) )
}
\keyword{Datasets}
\keyword{Ovarian cancer}
