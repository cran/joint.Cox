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
    \item{\code{CXCL12}}{a numeric vector of gene expressions.
 The CXCL12 gene expression is a predictive biomarker of survival in ovarian cancer (Popple et al. 2012).
 It has been known that CXCL12 promotes tumour growth, participates in tumour metastasis,
 and suppresses tumour immunity (Kryczek et al. 2007).
 The statistical significance of the CXCL12 expression on survival is first examined by Popple et al. (2012), and
 is further confirmed by Ganzfried et al. (2013) based on the meta-analysis of 14 independent studies.
 A meta-analysis using a joint model further confirmed that the expression of CXCL12 gene is predictive of
 both cancer relapse and death (Emura et al. 2017; 2018).}
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
    \item{\code{NCOA3}}{a numeric vector of gene expressions.
 The NCOA3 gene encodes a nuclear receptor coactivator, and amplification of the gene
 occurs in breast and ovarian cancers (Anzick et al. 1997).
 The overexpression of NCOA3 is associated with tumor size (Spears et al. 2012) and
 tamoxifen resistance (Osborne et al. 2003), which are involved in the progression.
 Yoshida et al. (2005) reported that NCOA3 could contribute to ovarian cancer progression by promoting cell migration.
 In Emura et al. (2018), the overexpression of the gene was highly associated with time-to-relapse
 (Coefficient=0.194, P-value<0.00001) and time-to-death (Coefficient=0.237, P-value<0.00001).
 This result is consistent with the function of these reports.}
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
    \item{\code{PDPN}}{a numeric vector of gene expressions. The PDPN gene encodes the podoplanin protein.
 It is reported that cancer cells with higher PDPN expression have higher malignant potential due to enhanced platelet aggregation,
 which promotes alteration of metastasis, cell motility, and epithelial-mesenchymal transition (Shindo et al. 2013).
 Zhang et al. (2011) reported that overexpression of PDPN in fibroblasts is significantly associated with a poor prognosis in ovarian carcinoma.
 In Emura et al. (2018), the overexpression of the gene was highly associated with time-to-relapse
 (Coefficient=0.222, P-value<0.00001) and time-to-death (Coefficient=0.161, P-value<0.0001).}
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
    \item{\code{TEAD1}}{a numeric vector of gene expressions.
 TEAD1 encodes a ubiquitous transcriptional enhancer factor that is a member of the TEA/ATTS domain family.
 It is reported that the protein level of TEAD1 was associated with poor prognosis in prostate cancer patients (Knight et al. 2008).
 In Emura et al. (2018), the overexpression of the gene was highly associated with
 time-to-relapse (Coefficient=0.195, P-value<0.00001) and time-to-death (Coefficient=0.223, P-value<0.00001).}
    \item{\code{TESK1}}{a numeric vector}
    \item{\code{TGM5}}{a numeric vector}
    \item{\code{THEMIS2}}{a numeric vector}
    \item{\code{TIMP2}}{a numeric vector of gene expressions. TIMP2 is a member of the TIMP gene family.
 The proteins encoded by this gene family are natural inhibitors of the matrix metalloproteinases (MMPs).
 MMPs and their inhibitors (TIMP gene family) play an important regulatory role in the homeostasis of
 the extracellular matrix (Halon et al. 2012). In addition to inhibitors of MMPs,
 TIMP2 has additional functions that are associated with cell proliferation and survival (Bourboulia et al., 2011).
 In Emura et al. (2018), the overexpression of the gene was highly associated with time-to-relapse (Coefficient=0.235, P-value<0.00001).}
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
    \item{\code{YWHAB}}{a numeric vector of gene expressions.
 YWHAB encodes a protein belonging to the 14-3-3 family of proteins,
 members of which mediate signal transduction by binding to phosphoserine-containing proteins.
 It is reported that the protein of YWHAB can regulate cell survival, proliferation, and motility (Tzivion 2006).
 Actually, it is reported that overexpression of this gene promotes tumor progression and was associated with
 extrahepatic metastasis and worse survival in hepatocellular carcinoma (Liu et al. 2011).
 In Emura et al. (2018), the overexpression of the gene was highly associated with
 time-to-relapse (Coefficient=0.169, P-value<0.0001) and time-to-death (Coefficient=0.263, P-value<0.00001)}
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
Ganzfried BF et al. (2013), Curated ovarian data: clinically annotated data for the ovarian cancer transcriptome, Database, Article ID bat013.
}
\references{
Bourboulia D, et al. (2011), Endogenous angiogenesis inhibitor blocks tumor growth via direct and indirect effects on tumor microenvironment. Am J Pathol 179:2589-600

Emura T, Nakatochi M, Murotani K, Rondeau V (2017), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Stat Methods Med Res 26(6):2649-66

Emura T, Nakatochi M, Matsui S, Michimae H, Rondeau V (2018), Personalized dynamic prediction of death according to tumour progression and high-dimensional
 genetic factors: meta-analysis with a joint model, Stat Methods Med Res 27(9):2842-58

Ganzfried BF et al. (2013), Curated ovarian data: clinically annotated data for the ovarian cancer transcriptome, Database, Article ID bat013.

Halon A, et al. (2012), Enhanced immunoreactivity of TIMP-2 in the stromal compartment of tumor as a marker of favorable prognosis in ovarian cancer patients.
 J Histochem Cytochem 60:491-501

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