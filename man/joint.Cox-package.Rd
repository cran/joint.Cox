\name{joint.Cox-package}
\alias{joint.Cox-package}
\alias{joint.Cox}
\docType{package}
\title{
Joint Frailty-Copula Models for Tumour Progression and Death in Meta-Analysis
}
\description{
Fit survival data and perform dynamic prediction under joint frailty-copula models for tumour progression and death.
 Likelihood-based methods are employed for estimating model parameters, where the baseline hazard functions are modeled by the cubic M-spline or the Weibull model.
 The methods are applicable for meta-analytic data containing individual-patient information from several studies.
 Survival outcomes need information on both terminal event time (e.g., time-to-death) and non-terminal event time (e.g., time-to-tumour progression).
 Methodologies were published in Emura et al. (2017), Emura et al. (2018), Emura et al. (2020), Wu et al. (2020), and Shinohara et al. (2020).
 See also the book of Emura et al. (2019).
 Survival data from ovarian cancer patients are also available.
}
\details{
\tabular{ll}{
Package: \tab joint.Cox\cr
Type: \tab Package\cr
Version: \tab 3.10\cr
Date: \tab 2020-12-8\cr
License: \tab GPL-2\cr
}
}
\author{
Takeshi Emura
Maintainer: Takeshi Emura <takeshiemura@gmail.com>
}
\references{
Emura T, Nakatochi M, Murotani K, Rondeau V (2017), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Stat Methods Med Res 26(6):2649-66

Emura T, Nakatochi M, Matsui S, Michimae H, Rondeau V (2018), Personalized dynamic prediction of death according to tumour progression and high-dimensional
 genetic factors: meta-analysis with a joint model, Stat Methods Med Res 27(9):2842-58

Emura T, Matsui S, Rondeau V (2019), Survival Analysis with Correlated Endpoints, Joint Frailty-Copula Models, JSS Research Series in Statistics, Springer

Emura T, Shih JH, Ha ID, Wilke RA (2020), Comparison of the marginal hazard model and the sub-distribution hazard model for competing risks under an assumed Copula,
 Stat Methods Med Res 29(8):2307-27

Shinohara S, Lin YH, Michimae H, Emura T (2020), Dynamic lifetime prediction using a Weibull-based bivariate failure time model:
 a meta-analysis of individual-patient data, Comm Stat-Simul, DOI:10.1080/03610918.2020.1855449

Wu BH, Michimae H, Emura T (2020), Meta-analysis of individual patient data with semi-competing risks under
 the Weibull joint frailty-copula model, Comp Stat 35(4):1525-52
}
\keyword{ package }

