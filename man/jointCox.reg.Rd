\name{jointCox.reg}
\alias{jointCox.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Penalized Likelihood Estimation under the Joint Cox Models Between Tumour Progression 
 and Death for Meta-Analysis
}
\description{
Perform regression analyses under a copula-based joint Cox proportional hazards model between tumour progression and death for meta-analysis,
 which is proposed by Emura et al. (2015).
 The method extends the joint frailty model of Rondeau et al. (2011 Stat Methods Med Res) such that intra-subject 
dependence is modeled via the Clayton copula.
 The method is applicable for meta-analysis combining several studies or for cluster survival data.
}
\usage{
jointCox.reg(t.event, event, t.death, death, Z1, Z2, group, alpha = 1,
 kappa_grid = c(seq(10, 1e+17, length = 30)), LCV_plot = TRUE, RNR_num = 10, Adj = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.event}{a vector object for time-to-tumour progression (TTP)}
  \item{event}{a vector object for progression indicator (=1 if progression; =0 if not progression)}
  \item{t.death}{a vector object for overall survival (OS), i.e., time-to-death}
  \item{death}{a vector object for death indicator(=1 if death; =0 if not death)}
  \item{Z1}{a matrix object for covariates associated with TTP; ncol(Z1)=the number of covariates}
  \item{Z2}{a matrix object for covariates associated with OS; ncol(Z2)=the number of covariates}
  \item{group}{a vector object for a group identification number, like 1,2,3....}
  \item{alpha}{A value related to the frailty (e.g., alpha=0 or =1); alpha=1 is default}
  \item{kappa_grid}{a vector for candidate smoothing parameters in likelihood cross-validation (LCV)}
  \item{LCV_plot}{Plot the LCV curves if "TRUE"}
  \item{RNR_num}{The number of randomizations for the ititial p0}
  \item{Adj}{Numerical adjustment to prevent overflow; Adj=500 is recommended}
}

\details{
Usually Z1=Z2 as described in Emura et al. (2015).
}
\value{
\item{count }{Count for event occurences}
\item{beta1 }{Regression coefficient for Z1}
\item{beta2 }{Regression coefficient for Z2}
\item{eta }{Frailty parameter}
\item{theta }{Copula parameter under the Clayton copula}
\item{tau }{Kendall's tau corresponding to the copula parameter}
\item{K1}{Smoothing parameter for TTP}
\item{K2}{Smoothing parameter for OS}
\item{g }{M spline coefficients for TTP}
\item{h }{M spline coefficients for OS}
\item{g_var }{Variance of M spline coefficients for TTP}
\item{h_var }{Variance of M spline coefficients for OS}
\item{convergence_code }{convergence code for "nlm" for maximizing the penalized likelihood}
\item{iterations }{the number of iterations until convergence}
\item{Randomiz_num }{the number of randomizing p0}
\item{gradient }{gradient vector; should be close to 0}
\item{ML}{Maximized penalized log-likelihood}
}

\references{
Emura T*, Nakatochi M, Murotani K, Rondeau V (2015), A joint frailty-copula model between tumour progression and death for meta-analysis,
 Statistical Methods in Medical Research, doi: 10.1177/0962280215604510
}
\author{Takeshi Emura}



\examples{
G=5
N=100
theta_true=2 ## copula parameter ##

beta1_true=1
beta2_true=1
r1_true=1
r2_true=1
eta_true=0.5
alpha_true=1

t.event=t.death=event=death=Z1=group=NULL

ij=0

set.seed(1)

for(i in 1:G){
  
  u_i=rgamma(1,shape=1/eta_true,scale=eta_true)
  for(j in 1:N){
    
    ij=ij+1
    
    group[ij]=i
    Z1[ij]=runif(1)
    r1_ij=r1_true*u_i*exp(beta1_true*Z1[ij])
    r2_ij=r2_true*(u_i^alpha_true)*exp(beta2_true*Z1[ij])
    V1=runif(1)
    V2=runif(1)
    X_ij=-1/r1_ij*log(1-V1);W=(1-V1)^(-theta_true)
    D_ij=1/theta_true/r2_ij*log(  1-W+W*(1-V2)^(-theta_true/(theta_true+1))  )
    C_ij=runif(1,min=0,max=5)
    t.event[ij]=min(X_ij,D_ij,C_ij)
    t.death[ij]=min(D_ij,C_ij)
    event[ij]=as.numeric( t.event[ij]==X_ij )
    death[ij]=as.numeric( t.death[ij]==D_ij )
  }
  
}

Z1=as.matrix(Z1)
#jointCox.reg(t.event=t.event,event=event,t.death=t.death,
#                          death=death,Z1=Z1,Z2=Z1,group=group,alpha=alpha_true,
#                          kappa_grid=seq(10,1000,length=10),LCV_plot=TRUE)

#jointCox.indep.reg(t.event=t.event,event=event,t.death=t.death,
#             death=death,Z1=Z1,Z2=Z1,group=group,alpha=alpha_true,
#             kappa_grid=seq(10,1000,length=10),LCV_plot=TRUE)

### takes a few minutes to finish ###
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Estimation }
\keyword{ Copula }% __ONLY ONE__ keyword per line
