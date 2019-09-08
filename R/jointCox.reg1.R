jointCox.reg1 <-
function(t.event,event,t.death,death,Z1,Z2,Z12,group,alpha=1,
         kappa1=c(seq(10,1e+17,length=30)),
         kappa2=c(seq(10,1e+17,length=30)),
         LCV.plot=TRUE,Randomize_num=10,Adj=500,convergence.par=FALSE){

T1=t.event
T2=t.death
d1=event
d2=death
Z1=as.matrix(Z1)
Z2=as.matrix(Z2)
p1=ncol(Z1)
p2=ncol(Z2)
p12=ncol(Z12)

G_id=as.numeric((levels(factor(group))))
G=length(G_id)

########### Summary ###########
n.event=xtabs(d1~group)
n.death=xtabs(d2~group)
n.censor=xtabs(1-d2~group)
count=cbind(table(group),n.event,n.death,n.censor)  
colnames(count)=c("No.of samples","No.of events","No.of deaths","No.of censors")

xi1=min( T1 )
xi3=max( T2 )

######## Choose smoothing parameters K1 and K2 ##########
res1=splineCox.reg(t.event,event,Z1,xi1=xi1,xi3=xi3,
              kappa=kappa1,LCV.plot=LCV.plot)
res2=splineCox.reg(t.death,death,Z2,xi1=xi1,xi3=xi3,
              kappa=kappa2,LCV.plot=LCV.plot)

K1_est=res1$kappa
LCV1_res=c(K1=K1_est,LCV1=res1$LCV)
K2_est=res2$kappa
LCV2_res=c(K2=K2_est,LCV2=res2$LCV)

Omega=c(192,-132,24,12,0,
        -132,96,-24,-12,12,
        24,-24,24,-24,24,
        12,-12,-24,96,-132,
        0,12,24,-132,192)
Omega=matrix(Omega,5,5)/( (xi3-xi1)/2 )^5

############ Likelihood function ###############
l.func=function(phi){

  g1=exp( pmax( pmin(phi[1:5],500), -500)  )
  g2=exp( pmax( pmin(phi[6:10],500),-500)  )
  eta=exp(phi[11])
  gamma0=phi[12]
  beta1=phi[(12+1):(12+p1)]
  beta2=phi[(12+p1+1):(12+p1+p2)]
  beta12=phi[(12+p1+p2+1):(12+p1+p2+p12)]
  theta=pmin( exp(gamma0+beta12*Z12),exp(5) )
  
  l=-K1_est*t(g1)%*%Omega%*%g1-K2_est*t(g2)%*%Omega%*%g2
  bZ1=as.vector( Z1%*%beta1 )
  bZ2=as.vector( Z2%*%beta2 )
  r1=as.vector( M.spline(T1,xi1=xi1,xi3=xi3)%*%g1 )
  r2=as.vector( M.spline(T2,xi1=xi1,xi3=xi3)%*%g2 )
  R1=as.vector( I.spline(T1,xi1=xi1,xi3=xi3)%*%g1 )
  R2=as.vector( I.spline(T2,xi1=xi1,xi3=xi3)%*%g2 )
  l=l+sum( d1*(log(r1)+bZ1) )+sum( d2*(log(r2)+bZ2) )

  for(i in G_id){
    
    Gi=c(group==i)
    m1=sum(d1[Gi])
    m2=sum(d2[Gi])
    #m12=sum(d1[Gi]*d2[Gi])
    EZ1=exp(bZ1[Gi])*R1[Gi]
    EZ2=exp(bZ2[Gi])*R2[Gi]
    D1=as.logical(d1[Gi])
    D2=as.logical(d2[Gi])
    
    func1=function(u){
      S1=pmin( exp( u%*%t( EZ1*theta[Gi] ) ), exp(500) )
      S2=pmin( exp( u^alpha%*%t( EZ2*theta[Gi] ) ), exp(500) )
      A=(S1+S2-1)
      E1=apply((S1/A)[,D1,drop=FALSE],MARGIN=1,FUN=prod)
      E2=apply((S2/A)[,D2,drop=FALSE],MARGIN=1,FUN=prod)
      Psi=log(A)%*%(1/theta[Gi])
      #Psi=rowSums( (1/t(theta[Gi]))*log(A) )
      D12=exp(-Psi+Adj)  ### Adjustment to avoid too small D12 ###
      u^(m1+alpha*m2)*E1*E2*D12*prod( (1+theta[Gi])^(d1[Gi]*d2[Gi]) )*dgamma(u,shape=1/eta,scale=eta)
    }
    
    Int=try( integrate(func1,0.001,10,stop.on.error = FALSE) ) 
    if( class(Int)=="try-error" ){l=l-500000}else{
      if(Int$value==0){l=l-500000}else{
        l=l+log(Int$value)-Adj ### Re-adjustment to avoid too small D12 ### 
      }
    }
    
  }
  
  -l  
}

p0=rep(0,12+p1+p2+p12)
res=nlm(l.func,p=p0,hessian=TRUE)
MPL=-res$minimum

R_num=0
repeat{
  if( (min( eigen(res$hessian)$values )>0)&(res$code==1) ){break}
  if(R_num>=Randomize_num){break}
  R_num=R_num+1
  p0_Rand=runif(12+p1+p2+p12,-1,1)
  res_Rand=nlm(l.func,p=p0_Rand,hessian=TRUE)
  MPL_Rand=-res_Rand$minimum
  if(MPL_Rand>MPL){
    res=res_Rand
    MPL=-res$minimum
  }
}
H_PL=-res$hessian
DF_upper=18+p1+p2+p12

temp=(det(H_PL)==0)|is.na(det(H_PL))
if(temp){V=solve( -H_PL+diag(rep(0.0001,12+p1+p2+p12)) ,tol=10^(-50))}else{  V=solve(-H_PL,tol=10^(-50))}

D_PL=diag( c(1/exp(res$estimate[1:12]),rep(1,p1+p2+p12)) )
H_PL=D_PL%*%H_PL%*%D_PL
H=H_PL
H[1:5,1:5]=H[1:5,1:5]+2*K1_est*Omega
H[6:10,6:10]=H[6:10,6:10]+2*K2_est*Omega
if( is.na(det(H_PL))|det(H_PL)==0 ){DF=DF_upper}else{
  DF=min( max( sum( diag(solve(H_PL,tol=10^(-50))%*%H) ), p1+p2+p12+2), DF_upper)
}
K1_est=K2_est=0
LCV=-l.func(res$estimate)-DF 

convergence_res=c(MPL=MPL,DF=DF,LCV=LCV,code=res$code,
          No.of.iterations=res$iterations,No.of.randomizations=R_num)


g_est=exp(res$est[1:5])
h_est=exp(res$est[6:10])
eta_est=exp(res$est[11])
theta0_est=exp(res$est[12])
gamma0_est=res$est[12]
tau0_est=theta0_est/(theta0_est+2)
beta1_est=res$est[(12+1):(12+p1)]
beta2_est=res$est[(12+p1+1):(12+p1+p2)]
beta12_est=res$est[(12+p1+p2+1):(12+p1+p2+p12)]

eta_se=eta_est*sqrt(diag(V)[11])
theta0_se=theta0_est*sqrt(diag(V)[12])
tau0_se=2/((theta0_est+2)^2)*theta0_se
beta1_se=sqrt(diag(V)[(12+1):(12+p1)])
beta2_se=sqrt(diag(V)[(12+p1+1):(12+p1+p2)])
beta12_se=sqrt(diag(V)[(12+p1+p2+1):(12+p1+p2+p12)])

theta0_Lower=theta0_est*exp(-1.96*sqrt(diag(V)[12]))
theta0_Upper=theta0_est*exp(1.96*sqrt(diag(V)[12]))

g_var=diag(g_est)%*%V[1:5,1:5]%*%diag(g_est)
h_var=diag(h_est)%*%V[6:10,6:10]%*%diag(h_est)

beta1_res=c(estimate=beta1_est,SE=beta1_se,
  Lower=beta1_est-1.96*beta1_se,Upper=beta1_est+1.96*beta1_se)
beta2_res=c(estimate=beta2_est,SE=beta2_se,
            Lower=beta2_est-1.96*beta2_se,Upper=beta2_est+1.96*beta2_se)
eta_res=c(estimate=eta_est,SE=eta_se,
          Lower=eta_est*exp(-1.96*sqrt(diag(V)[11])),
          Upper=eta_est*exp(1.96*sqrt(diag(V)[11])))
theta0_res=c(estimate=theta0_est,SE=theta0_se,Lower=theta0_Lower,Upper=theta0_Upper)
tau0_res=c(estimate=tau0_est,tau_se=tau0_se,
          Lower=theta0_Lower/(theta0_Lower+2),Upper=theta0_Upper/(theta0_Upper+2))
beta12_res=c(estimate=beta12_est,SE=beta12_se,
            Lower=beta12_est-1.96*beta12_se,Upper=beta12_est+1.96*beta12_se)

if(convergence.par==FALSE){convergence.parameters=NULL}else{
   convergence.parameters=list(log_estimate=res$est,gradient=-res$gradient,log_var=V)
}

list(count=count,
     beta1=beta1_res,beta2=beta2_res,eta=eta_res,theta=theta0_res,
     tau=tau0_res,beta12=beta12_res,
     LCV1=LCV1_res,LCV2=LCV2_res,
     g=g_est,h=h_est,g_var=g_var,h_var=h_var,
     convergence=convergence_res,convergence.parameters=convergence.parameters
     )
}
