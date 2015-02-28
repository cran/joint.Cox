\name{jointCox.reg}
\alias{jointCox.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Penalized Likelihood Estimation under the Joint Cox Models Between TTP 
 and OS for Meta-Analysis
}
\description{
Perform regression analyses under a copula-based joint Cox proportional hazards model between TTP and OS for meta-analysis,
 which is proposed by Emura, Murotani & Rondeau (2015).
 The method extends the joint frailty model of Rondeau et al. (2011 Stat Methods Med Res) such that intra-subject 
dependence is modeled via the Clayton copula.
 The method is applicable for meta-analysis combining several studies or for cluster survival data.
Below, the symbol G represents the number of studies or clusters.   
}
\usage{
jointCox.reg(G, T1, d1, T2, d2, Z1, Z2, alpha = 1, p0 = rep(0, 14),
 K_low = 50, K_up = 1000, LCV_plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{The number of studies, which must be the number of objects in T1}
  \item{T1}{A list object for time-to-progression (TTP)}
  \item{d1}{A list object for the censoring indicator of TTP}
  \item{T2}{A list object for time-to-death (OS)}
  \item{d2}{A list object for the censoring indicator of OS}
  \item{Z1}{A list object for a covariate associated with TTP}
  \item{Z2}{A list object for a covariate associated with OS; often Z1=Z2}
  \item{alpha}{A nonnegative value (e.g., alpha=0 or =1); alpha=1 is default}
  \item{p0}{A vector of 14 initial values used for maximizing the log-likelihood; detault is rep(0,14)}
  \item{K_low}{The lower limit of tuning parameters in LCV; default=50}
  \item{K_up}{The upper limit of tuning parameters in LCV; default=1000}
  \item{LCV_plot}{Plot the LCV curves if "TRUE"}
}
\details{
Currently, only one covariate is allowed (p=1).
}
\value{
\item{beta1 }{Regression coefficient for Z1}
\item{beta2 }{Regression coefficient for Z2}
\item{eta }{Frailty parameter}
\item{theta }{Copula parameter under the Clayton copula}
\item{tau }{Kendall's tau corresponding to the copula parameter}
\item{K1}{Smoothing parameter for TTP}
\item{K2 }{Smoothing parameter for OS}
\item{g }{M spline coefficients for TTP}
\item{h }{M spline coefficients for OS}
\item{g_var }{Variance of M spline coefficients for TTP}
\item{h_var }{Variance of M spline coefficients for OS}
}
\references{
Emura T*, Murotani K, Rondeau V (2015), A joint frailty-copula model between disease progression and death for meta-analysis (to be submitted)
}
\author{Takeshi Emura}

\examples{

G=5
N=50
theta_true=2 ## copula parameter ##

beta1_true=1
beta2_true=1
r1_true=1
r2_true=1

eta_true=0.5
alpha_true=1

T1=T2=d1=d2=Z1=rep(list(rep(0,N)),G)

set.seed(2)
for(i in 1:G){
  
  u_i=rgamma(1,shape=1/eta_true,scale=eta_true)
  for(j in 1:N){
    
    Z1[[i]][j]=runif(1)
    r1_ij=r1_true*u_i*exp(beta1_true*Z1[[i]][j])
    r2_ij=r2_true*(u_i^alpha_true)*exp(beta2_true*Z1[[i]][j])
    V1=runif(1)
    V2=runif(1)
    X_ij=-1/r1_ij*log(1-V1);W=(1-V1)^(-theta_true)
    D_ij=1/theta_true/r2_ij*log(  1-W+W*(1-V2)^(-theta_true/(theta_true+1))  )
    C_ij=runif(1,min=0,max=5)
    T1[[i]][j]=min(X_ij,D_ij,C_ij)
    T2[[i]][j]=min(D_ij,C_ij)
    d1[[i]][j]=as.numeric( T1[[i]][j]==X_ij )
    d2[[i]][j]=as.numeric( T2[[i]][j]==D_ij )
  }
  
}

T1
T2
d1
d2
Z1
#jointCox.reg(G,T1,d1,T2,d2,Z1,Z1,alpha=alpha_true)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
