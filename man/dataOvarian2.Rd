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
    \item{\code{CXCL12}}{a numeric vector of gene expressions.
 The CXCL12 gene expression is a predictive biomarker of survival in ovarian cancer (Popple et al. 2012).
 It has been known that CXCL12 promotes tumour growth, participates in tumour metastasis,
 and suppresses tumour immunity (Kryczek et al. 2007).
 The statistical significance of the CXCL12 expression on survival is first examined by Popple et al. (2012), and
 is further confirmed by Ganzfried et al. (2013) based on the meta-analysis of 14 independent studies.
 A meta-analysis using a joint model further confirmed that the expression of CXCL12 gene is predictive of
 both cancer relapse and death (Emura et al. 2017; 2018)}
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
    \item{\code{NCOA3}}{a numeric vector of gene expressions.
 The NCOA3 gene encodes a nuclear receptor coactivator, and amplification of the gene
 occurs in breast and ovarian cancers (Anzick et al. 1997).
 The overexpression of NCOA3 is associated with tumor size (Spears et al. 2012) and
 tamoxifen resistance (Osborne et al. 2003), which are involved in the progression.
 Yoshida et al. (2005) reported that NCOA3 could contribute to ovarian cancer progression by promoting cell migration.
 In Emura et al. (2018), the overexpression of the gene was highly associated with time-to-relapse
 (Coefficient=0.194, P-value<0.00001) and time-to-death (Coefficient=0.237, P-value<0.00001).
 This result is consistent with the function of these reports.}
    \item{\code{NCOA6}}{a numeric vector of gene expressions}
    \item{\code{NOTCH2NL}}{a numeric vector}
    \item{\code{NR1H3}}{a numeric vector}
    \item{\code{NUAK1}}{a numeric vector}
    \item{\code{OAT}}{a numeric vector}
    \item{\code{OMD}}{a numeric vector}
    \item{\code{PAK4}}{a numeric vector}
    \item{\code{PCDH9}}{a numeric vector}
    \item{\code{PDP1}}{a numeric vector}
    \item{\code{PDPN}}{a numeric vector of gene expressions. The PDPN gene encodes the podoplanin protein.
 It is reported that cancer cells with higher PDPN expression have higher malignant potential due to enhanced platelet aggregation,
 which promotes alteration of metastasis, cell motility, and epithelial-mesenchymal transition (Shindo et al. 2013).
 Zhang et al. (2011) reported that overexpression of PDPN in fibroblasts is significantly associated with a poor prognosis in ovarian carcinoma.
 In Emura et al. (2018), the overexpression of the gene was highly associated with time-to-relapse
 (Coefficient=0.222, P-value<0.00001) and time-to-death (Coefficient=0.161, P-value<0.0001).}
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
    \item{\code{TEAD1}}{a numeric vector of gene expressions.
 TEAD1 encodes a ubiquitous transcriptional enhancer factor that is a member of the TEA/ATTS domain family.
 It is reported that the protein level of TEAD1 was associated with poor prognosis in prostate cancer patients (Knight et al. 2008).
 In Emura et al. (2018), the overexpression of the gene was highly associated with
 time-to-relapse (Coefficient=0.195, P-value<0.00001) and time-to-death (Coefficient=0.223, P-value<0.00001).}
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
    \item{\code{YWHAB}}{a numeric vector of gene expressions.
 YWHAB encodes a protein belonging to the 14-3-3 family of proteins,
 members of which mediate signal transduction by binding to phosphoserine-containing proteins.
 It is reported that the protein of YWHAB can regulate cell survival, proliferation, and motility (Tzivion 2006).
 Actually, it is reported that overexpression of this gene promotes tumor progression and was associated with
 extrahepatic metastasis and worse survival in hepatocellular carcinoma (Liu et al. 2011).
 In Emura et al. (2018), the overexpression of the gene was highly associated with
 time-to-relapse (Coefficient=0.169, P-value<0.0001) and time-to-death (Coefficient=0.263, P-value<0.00001). }
    \item{\code{ZFP36}}{a numeric vector}
    \item{\code{ZFP36L2}}{a numeric vector}
    \item{\code{ZNF148}}{a numeric vector}
  }
}
\details{
4 studies are combined (group=4, 9, 12, and 16). The numbers 4, 9, 12 and 16 corresponds to the IDs 
from the original data of Ganzfried et al. (2013).}
\source{
Ganzfried BF et al. (2013), Curated ovarian data: clinically annotated data for the ovarian cancer transcriptome, Database, Article ID bat013
}
\references{
Emura T, Nakatochi M, Murotani K, Rondeau V (2017), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Stat Methods Med Res 26(6):2649-66

Emura T, Nakatochi M, Matsui S, Michimae H, Rondeau V (2018), Personalized dynamic prediction of death according to tumour progression and high-dimensional
 genetic factors: meta-analysis with a joint model, Stat Methods Med Res 27(9):2842-58

Ganzfried BF, et al. (2013), Curated ovarian data: clinically annotated data for the ovarian cancer transcriptome, Database, Article ID bat013.

Knight JF, et al. (2008), TEAD1 and c-Cbl are novel prostate basal cell markers that correlate with poor clinical outcome in prostate cancer. Br J Cancer 99:1849-58

Kryczek I, et al. (2007), Stroma-derived factor (SDF-1/CXCL12) and human tumor pathogenesis. Am J Physiol 292:987-95

Liu TA, et al. (2011), Increased expression of 14-3-3beta promotes tumor progression and predicts extrahepatic metastasis and worse survival
 in hepatocellular carcinoma. Am J Pathol 179:2698-708

Osborne CK, et al. (2003), Role of the estrogen receptor coactivator AIB1 (SRC-3) and HER-2/neu in tamoxifen resistance in breast cancer. J Natl Cancer Inst 95:353-61

Popple A, et al. (2012), The chemokine, CXCL12, is an independent predictor of poor survival in ovarian cancer. Br J Cancer 106:1306-13

Shindo K, et al. (2013), Podoplanin expression in cancer-associated fibroblasts enhances tumor progression of invasive ductal carcinoma of the pancreas. Mol Cancer 12:168

Tzivion G, et al. (2006), 14-3-3 proteins as potential oncogenes. Semin Cancer Biol 16:203-13

Yoshida H, et al. (2005), Steroid receptor coactivator-3, a homolog of Taiman that controls cell migration in the Drosophila ovary,
 regulates migration of human ovarian cancer cells. Mol Cell Endocrinol 245:77-85

Zhang Y, et al. (2011), Ovarian cancer-associated fibroblasts contribute to epithelial ovarian carcinoma metastasis by promoting angiogenesis,
 lymphangiogenesis and tumor cell invasion. Cancer Lett 303:47-55
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
