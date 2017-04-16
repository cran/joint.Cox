\name{dataOvarian1}
\alias{dataOvarian1}
\docType{data}
\title{
Data on time-to-recurrence and 158 gene expressions for 912 ovarian cancer patients from 4 independent studies. 
}
\description{
Meta-analytic data containing 158 gene expressions and time-to-relapse information for ovarian cancer patients.
 The data include time-to-recurrence, residual tumour size (>=1cm> vs. <1cm), and associated 158 gene expressions.
 The dataset is a subset of the curated ovarian data of Ganzfried et al (2013).
 We prepared the dataset by using "patientselection.config" in "Curated ovarian data" around October 2016. 
}
\usage{data("dataOvarian1")}
\format{
  A data frame with 912 observations on the following 162 variables.
  \describe{
    \item{\code{t.event}}{: time-to-recurrence in days}
    \item{\code{event}}{: event indicator (1=recurrence, 0=no recurrence)}
    \item{\code{group}}{: study ID; group=4, 9, 12, or 16}
    \item{\code{debulk}}{: residual tumour size (>=1cm> vs. <1cm)}
    \item{\code{ABI3BP}}{a numeric vector}
    \item{\code{ADAM12}}{a numeric vector}
    \item{\code{ADORA3}}{a numeric vector}
    \item{\code{ANKRD27}}{a numeric vector}
    \item{\code{AP2M1}}{a numeric vector}
    \item{\code{AP3S1}}{a numeric vector}
    \item{\code{ARHGAP28}}{a numeric vector}
    \item{\code{ARHGAP29}}{a numeric vector}
    \item{\code{ARTN}}{a numeric vector}
    \item{\code{ASAP3}}{a numeric vector}
    \item{\code{B4GALT5}}{a numeric vector}
    \item{\code{BCAP31}}{a numeric vector}
    \item{\code{BRD4}}{a numeric vector}
    \item{\code{C1QTNF3}}{a numeric vector}
    \item{\code{CALD1}}{a numeric vector}
    \item{\code{CCNE1}}{a numeric vector}
    \item{\code{CCNL1}}{a numeric vector}
    \item{\code{CDC42}}{a numeric vector}
    \item{\code{CDV3}}{a numeric vector}
    \item{\code{CEBPB}}{a numeric vector}
    \item{\code{CLIC4}}{a numeric vector}
    \item{\code{COL10A1}}{a numeric vector}
    \item{\code{COL11A1}}{a numeric vector}
    \item{\code{COL16A1}}{a numeric vector}
    \item{\code{COL3A1}}{a numeric vector}
    \item{\code{COL5A1}}{a numeric vector}
    \item{\code{COL5A2}}{a numeric vector}
    \item{\code{COMP}}{a numeric vector}
    \item{\code{CRISPLD2}}{a numeric vector}
    \item{\code{CRYAB}}{a numeric vector}
    \item{\code{CSE1L}}{a numeric vector}
    \item{\code{CTSK}}{a numeric vector}
    \item{\code{CXCL12}}{a numeric vector}
    \item{\code{CYR61}}{a numeric vector}
    \item{\code{DCUN1D1}}{a numeric vector}
    \item{\code{DDX27}}{a numeric vector}
    \item{\code{DIAPH3}}{a numeric vector}
    \item{\code{DNAJB4}}{a numeric vector}
    \item{\code{DNAJC13}}{a numeric vector}
    \item{\code{DNAJC8}}{a numeric vector}
    \item{\code{DPYSL3}}{a numeric vector}
    \item{\code{DVL3}}{a numeric vector}
    \item{\code{EFNB2}}{a numeric vector}
    \item{\code{EIF3K}}{a numeric vector}
    \item{\code{ELK1}}{a numeric vector}
    \item{\code{ENPP1}}{a numeric vector}
    \item{\code{EPYC}}{a numeric vector}
    \item{\code{FABP4}}{a numeric vector}
    \item{\code{FAM69A}}{a numeric vector}
    \item{\code{FAP}}{a numeric vector}
    \item{\code{FERMT2}}{a numeric vector}
    \item{\code{FGF1}}{a numeric vector}
    \item{\code{FN1}}{a numeric vector}
    \item{\code{FOSL2}}{a numeric vector}
    \item{\code{FSTL1}}{a numeric vector}
    \item{\code{GABRG3}}{a numeric vector}
    \item{\code{GAS1}}{a numeric vector}
    \item{\code{GFRA1}}{a numeric vector}
    \item{\code{GFRA3}}{a numeric vector}
    \item{\code{GJC1}}{a numeric vector}
    \item{\code{GLIPR1}}{a numeric vector}
    \item{\code{GPATCH1}}{a numeric vector}
    \item{\code{HLTF}}{a numeric vector}
    \item{\code{HP1BP3}}{a numeric vector}
    \item{\code{HSD17B6}}{a numeric vector}
    \item{\code{INHBA}}{a numeric vector}
    \item{\code{ITGB1}}{a numeric vector}
    \item{\code{JUN}}{a numeric vector}
    \item{\code{KIAA0226}}{a numeric vector}
    \item{\code{KIAA0355}}{a numeric vector}
    \item{\code{KIAA1598}}{a numeric vector}
    \item{\code{KIN}}{a numeric vector}
    \item{\code{KLHL25}}{a numeric vector}
    \item{\code{KPNA6}}{a numeric vector}
    \item{\code{KRT7}}{a numeric vector}
    \item{\code{KRTAP5.8}}{a numeric vector}
    \item{\code{L2HGDH}}{a numeric vector}
    \item{\code{LGALS1}}{a numeric vector}
    \item{\code{LOX}}{a numeric vector}
    \item{\code{LPP}}{a numeric vector}
    \item{\code{LUM}}{a numeric vector}
    \item{\code{LUZP1}}{a numeric vector}
    \item{\code{MAP7D1}}{a numeric vector}
    \item{\code{MAPRE1}}{a numeric vector}
    \item{\code{MCL1}}{a numeric vector}
    \item{\code{MEOX2}}{a numeric vector}
    \item{\code{METTL9}}{a numeric vector}
    \item{\code{MFN1}}{a numeric vector}
    \item{\code{MICAL2}}{a numeric vector}
    \item{\code{MMP12}}{a numeric vector}
    \item{\code{MRPS22}}{a numeric vector}
    \item{\code{MXD1}}{a numeric vector}
    \item{\code{MXRA8}}{a numeric vector}
    \item{\code{N4BP2L2}}{a numeric vector}
    \item{\code{NCOA3}}{a numeric vector}
    \item{\code{NDRG3}}{a numeric vector}
    \item{\code{NINJ1}}{a numeric vector}
    \item{\code{NNMT}}{a numeric vector}
    \item{\code{NOTCH2}}{a numeric vector}
    \item{\code{NPY}}{a numeric vector}
    \item{\code{NTM}}{a numeric vector}
    \item{\code{NUAK1}}{a numeric vector}
    \item{\code{OAT}}{a numeric vector}
    \item{\code{OLFML2B}}{a numeric vector}
    \item{\code{PARD3}}{a numeric vector}
    \item{\code{PCYT1A}}{a numeric vector}
    \item{\code{PDE1A}}{a numeric vector}
    \item{\code{PDGFD}}{a numeric vector}
    \item{\code{PDPN}}{a numeric vector}
    \item{\code{PGRMC1}}{a numeric vector}
    \item{\code{PLAU}}{a numeric vector}
    \item{\code{PLOD2}}{a numeric vector}
    \item{\code{PLSCR4}}{a numeric vector}
    \item{\code{POSTN}}{a numeric vector}
    \item{\code{PPIC}}{a numeric vector}
    \item{\code{PRDM2}}{a numeric vector}
    \item{\code{PSMC4}}{a numeric vector}
    \item{\code{RAB22A}}{a numeric vector}
    \item{\code{RAB31}}{a numeric vector}
    \item{\code{RAB32}}{a numeric vector}
    \item{\code{RARRES1}}{a numeric vector}
    \item{\code{RPS16}}{a numeric vector}
    \item{\code{SERPINE1}}{a numeric vector}
    \item{\code{SGK1}}{a numeric vector}
    \item{\code{SH3PXD2A}}{a numeric vector}
    \item{\code{SKIL}}{a numeric vector}
    \item{\code{SLC12A8}}{a numeric vector}
    \item{\code{SPARC}}{a numeric vector}
    \item{\code{SPHK1}}{a numeric vector}
    \item{\code{STAU1}}{a numeric vector}
    \item{\code{SULF1}}{a numeric vector}
    \item{\code{SUPT5H}}{a numeric vector}
    \item{\code{TAGLN}}{a numeric vector}
    \item{\code{TBCB}}{a numeric vector}
    \item{\code{TEAD1}}{a numeric vector}
    \item{\code{TESK1}}{a numeric vector}
    \item{\code{TGM5}}{a numeric vector}
    \item{\code{THEMIS2}}{a numeric vector}
    \item{\code{TIMP2}}{a numeric vector}
    \item{\code{TIMP3}}{a numeric vector}
    \item{\code{TJP1}}{a numeric vector}
    \item{\code{TP73.AS1}}{a numeric vector}
    \item{\code{TPM2}}{a numeric vector}
    \item{\code{TPM4}}{a numeric vector}
    \item{\code{TSC22D2}}{a numeric vector}
    \item{\code{TUBB2A}}{a numeric vector}
    \item{\code{TUBB6}}{a numeric vector}
    \item{\code{TUFT1}}{a numeric vector}
    \item{\code{URI1}}{a numeric vector}
    \item{\code{USP48}}{a numeric vector}
    \item{\code{VCAN}}{a numeric vector}
    \item{\code{VSIG4}}{a numeric vector}
    \item{\code{YWHAB}}{a numeric vector}
    \item{\code{ZFP36}}{a numeric vector}
    \item{\code{ZFP36L2}}{a numeric vector}
    \item{\code{ZMYM1}}{a numeric vector}
    \item{\code{ZNF148}}{a numeric vector}
    \item{\code{ZNF79}}{a numeric vector}
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
data(dataOvarian1)
######## univariate Cox ##########
t.event=dataOvarian1$t.event
event=dataOvarian1$event
X.mat=dataOvarian1[,-c(1,2,3,4)] ## gene expression
Symbol=colnames(dataOvarian1)[-c(1,2,3,4)] ## gene symbol

p=ncol(X.mat)
P_value=coef=NULL
for(j in 1:p){
  res=summary(coxph(Surv(t.event,event)~X.mat[,j]))$coefficients
  P_value=c(P_value,res[5])
  coef=c(coef,res[1])
}
data.frame( gene=Symbol[order(P_value)], P=P_value[order(P_value)], 
coef=round(coef[order(P_value)],3) )
}
\keyword{Datasets}
\keyword{Ovarian cancer}