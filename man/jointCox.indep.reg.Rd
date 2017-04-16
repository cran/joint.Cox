\name{jointCox.indep.reg}
\alias{jointCox.indep.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Penalized Likelihood Estimation under the Joint Cox Models Between Tumour Progression 
 and Death for Meta-Analysis
}
\description{
Perform regression analyses under a joint Cox proportional hazards model between tumour progression and death for meta-analysis,
 which is proposed by Rondeau et al. (2015). 
 The method is applicable for meta-analysis combining several studies or for cluster survival data.
}
\usage{
jointCox.indep.reg(t.event, event, t.death, death, Z1, Z2, group, alpha = 1,
 kappa_grid = c(seq(10, 1e+17, length = 30)), LCV_plot = TRUE,
 Randomize_num = 10, Adj = 500,convergence.par=FALSE)
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
  \item{Randomize_num}{The number of randomizations for the ititial p0}
  \item{Adj}{Numerical adjustment to prevent overflow; Adj=500 is recommended}
  \item{convergence.par}{If TRUE, the converged estimate, gradient, and Hessian matrix are given (log-transformed)}
}
\details{
We employ "nlm" routine to maximize the penalized likelihood function with the initial value described in Emura et al. (2015).
If "nlm" does not converge, then we randomize the initial value by adding uniform random variables (Hu and Emura, 2015).
}
\value{
\item{count }{Count for event occurences}
\item{beta1 }{Regression coefficient for Z1}
\item{beta2 }{Regression coefficient for Z2}
\item{eta }{Frailty parameter (variance)}
\item{LCV1}{Likelihood cross-validation for TTP}
\item{LCV2}{Likelihood cross-validation for OS}
\item{g }{M-spline coefficients for TTP}
\item{h }{M-spline coefficients for OS}
\item{g_var }{Variance of M-spline coefficients for TTP}
\item{h_var }{Variance of M-spline coefficients for OS}
\item{convergence }{convergence results for maximizing penalized likelihood}
\item{convergence.parameters }{converged estimate, gradient, and Hessian matrix (log-transformed)}
}

\references{
Rondeau V, Pignon JP, Michiels S (2015). A joint model for dependence between clustered times to tumour progression and deaths: A meta-analysis of chemotherapy in head and neck cancer.
 Statistical Method in Medical Research 24(6):711-729.

Hu YH, Emura T (2015), Maximum likelihood estimation for a special exponential family under random double-truncation,
 Computational Statistics, 30 (No. 4): 1199-1229
}
\author{Takeshi Emura}

\section{Error }{"Error in integrate(func1, 0.001, 10, stop.on.error = FALSE):non-finite function value",
 an error occurring when the penalized likelihood is maximizad by "nlm". The error frequently
 occurs during the iterations for maximizing the penalized likelihood. The error is not crucial and can be ignored.}
 
\section{Warning }{"NA/Inf replaced by maximum positive value", an error occurring when the penalized likelihood
 is maximizad by "nlm". The error frequently occurs during the iterations for maximizing the penalized likelihood.
 The error is not crucial and can be ignored.}

\examples{
G=5
N=20
theta_true=2 ## copula parameter ##

beta1_true=2
beta2_true=2
r1_true=1
r2_true=1
eta_true=1
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
    C_ij=runif(1,min=0,max=10)
    t.event[ij]=min(X_ij,D_ij,C_ij)
    t.death[ij]=min(D_ij,C_ij)
    event[ij]=as.numeric( t.event[ij]==X_ij )
    death[ij]=as.numeric( t.death[ij]==D_ij )
  }
  
}

Z1=as.matrix(Z1)

jointCox.indep.reg(t.event=t.event,event=event,t.death=t.death,
             death=death,Z1=Z1,Z2=Z1,group=group,alpha=alpha_true,
             kappa_grid=seq(10,5000,length=10),LCV_plot=TRUE)

############## Reproduce the results of Emura et al. (2015) #############
data(dataOvarian)
t.event=dataOvarian$t.event
event=dataOvarian$event
t.death=dataOvarian$t.death
death=dataOvarian$death
Z1=as.matrix(dataOvarian$CXCL12)
group=dataOvarian$group
alpha_given=0
kappa_grid=seq(10,1e+17,length=30)
set.seed(1)
#jointCox.indep.reg(t.event=t.event,event=event,t.death=t.death,death=death,
#                 Z1=Z1,Z2=Z1,group=group,alpha=alpha_given,
#                 kappa_grid=kappa_grid,LCV_plot=TRUE,Adj=500)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Estimation }
\keyword{ Frailty }% __ONLY ONE__ keyword per line
