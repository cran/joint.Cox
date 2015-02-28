jointCox.reg <-
function(G,T1,d1,T2,d2,Z1,Z2,alpha=1,p0=rep(0,14),K_low=50,K_up=1000,LCV_plot=TRUE){
  
########### M-spline matrix #############
xi1=min( unlist(T1) )
xi3=max( unlist(T2) )
D=(xi3-xi1)/2
xi2=xi1+D

M_func=function(t){
  z1=(t-xi1)/D;z2=(t-xi2)/D;z3=(t-xi3)/D
  
  M1=-(4*z2^3/D)*(t<xi2)+0*(t>=xi2)
  M2=(7*z1^3-18*z1^2+12*z1)/2/D*(t<xi2)-z3^3/2/D*(t>=xi2)
  M3=(-2*z1^3+3*z1^2)/D*(t<xi2)+(2*z2^3-3*z2^2+1)/D*(t>=xi2)
  M4=z1^3/2/D*(t<xi2)+(-7*z2^3+3*z2^2+3*z2+1)/2/D*(t>=xi2)
  M5=4*z2^3/D*(t>=xi2)
  
  cbind(M1,M2,M3,M4,M5)
}

########## I-spline matrix ##############
I_func=function(t){
  z1=(t-xi1)/D;z2=(t-xi2)/D;z3=(t-xi3)/D
  
  I1=(1-z2^4)*(t<xi2)+1*(t>=xi2)
  I2=( 7/8*z1^4-3*z1^3+3*z1^2 )*(t<xi2)+( 1-z3^4/8 )*(t>=xi2)
  I3=( -z1^4/2+z1^3 )*(t<xi2)+( 1/2+z2^4/2-z2^3+z2 )*(t>=xi2)
  I4=( z1^4/8 )*(t<xi2)+( 1/8-7/8*z2^4+1/2*z2^3+3/4*z2^2+1/2*z2 )*(t>=xi2)
  I5=z2^4*(t>=xi2)
  
  cbind(I1,I2,I3,I4,I5)
}

######## Penalization term #########
Omega=c(192,-132,24,12,0,
        -132,96,-24,-12,12,
        24,-24,24,-24,24,
        12,-12,-24,96,-132,
        0,12,24,-132,192)
Omega=matrix(Omega,5,5)/D^5

############ LCV for K1 ###############
T1_vec=unlist(T1)
Z1_vec=unlist(Z1)
d1_vec=unlist(d1)

l1.func=function(phi){
  beta1=phi[1]
  g1=exp(  pmin(phi[c(2,3,4,5,6)],500)  ) ## M-spline coefficients ##
  l=-K1*t(g1)%*%Omega%*%g1
  r1=as.vector( M_func(T1_vec)%*%g1 )
  R1=as.vector( I_func(T1_vec)%*%g1 )
  l=l+sum( d1_vec*(log(r1)+beta1*Z1_vec) )
  l=l-sum(  pmin( exp(beta1*Z1_vec)*R1, exp(500) )  )
  -l  
}

p1=rep(0,6)
K1=K_low
repeat{
  res1=nlm(l1.func,p=p1,hessian=TRUE)
  if( res1$code!=2 ){break}
  p1=runif(6,-1,1)
}
H1=res1$hessian

K_grid=seq(K_low,K_up,length=10)
L1=DF1=NULL

for(k in 1:length(K_grid)){
  K1=K_grid[k]
  res1=nlm(l1.func,p=rep(0,6),hessian=TRUE)
  H1_PL=res1$hessian
  K1=K_low
  L1[k]=-l1.func(res1$estimate) 
  if( is.na(det(H1_PL))|det(H1_PL)==0 ){DF1[k]=10}else{
    if(det(H1_PL)<=0.000000000000001){H1_PL=H1_PL+diag(rep(0.0001,6))}  
    DF1[k]=max( sum( diag(solve(H1_PL)%*%H1) ), 3)
  }
}

K1_est=K_grid[L1-DF1==max(L1-DF1)][1]

############ LCV for K2 ###############
T2_vec=unlist(T2)
Z2_vec=unlist(Z2)
d2_vec=unlist(d2)

l2.func=function(phi){
  beta2=phi[1]
  g2=exp(  pmin(phi[c(2,3,4,5,6)],500)  ) ## M-spline coefficients ##
  l=-K2*t(g2)%*%Omega%*%g2
  r2=as.vector( M_func(T2_vec)%*%g2 )
  R2=as.vector( I_func(T2_vec)%*%g2 )
  l=l+sum( d2_vec*(log(r2)+beta2*Z2_vec) )
  l=l-sum(  pmin( exp(beta2*Z2_vec)*R2, exp(500) )  )
  -l  
}

p2=rep(0,6)
K2=K_low
repeat{
  res2=nlm(l2.func,p=p2,hessian=TRUE)
  if( res2$code!=2 ){break}
  p2=runif(6,-1,1)
}
H2=res2$hessian

K_grid=seq(K_low,K_up,length=10)
L2=DF2=NULL

for(k in 1:length(K_grid)){
  K2=K_grid[k]
  res2=nlm(l2.func,p=rep(0,6),hessian=TRUE)
  H2_PL=res2$hessian
  K2=K_low
  L2[k]=-l2.func(res2$estimate) 
  if( is.na(det(H2_PL))|det(H2_PL)==0 ){DF2[k]=10}else{
    if(det(H2_PL)<=0.000000000000001){H2_PL=H2_PL+diag(rep(0.0001,6))}
    DF2[k]=max( sum( diag(solve(H2_PL)%*%H2) ), 3)
  }
}

K2_est=K_grid[L2-DF2==max(L2-DF2)][1]

########## Plotting LCV ##########
if(LCV_plot==TRUE){
  par(mfrow=c(1,3))
  plot(K_grid,L1,xlab="K1",ylab="logL",type="b",lwd=3)
  plot(K_grid,DF1,xlab="K1",ylab="DF",,type="b",lwd=3)
  plot(K_grid,L1-DF1,xlab="K1",ylab="LCV=logL-DF",type="b",lwd=3)
  points(K1_est,max(L1-DF1),xlab="K1",col="red",pch=17,cex=2)

  plot(K_grid,L2,xlab="K2",ylab="logL",type="b",lwd=3)
  plot(K_grid,DF2,xlab="K2",ylab="DF",type="b",lwd=3)
  plot(K_grid,L2-DF2,xlab="K2",ylab="LCV=logL-DF",type="b",lwd=3)
  points(K2_est,max(L2-DF2),col="red",pch=17,cex=2)
}
############ Likelihood function ###############
l.func=function(phi){
  
  beta1=phi[1]
  beta2=phi[2]
  g1=exp(  pmin(phi[c(3,4,5,6,7)],500)  ) ## M-spline coefficients ##
  g2=exp(  pmin(phi[c(8,9,10,11,12)],500)  ) ## M-spline coefficients ##
  eta=exp(phi[13])
  theta=min( exp(phi[14]),exp(5) )
  
  l=-K1_est*t(g1)%*%Omega%*%g1-K2_est*t(g2)%*%Omega%*%g2
  
  for(i in 1:G){
    
    r1=as.vector( M_func(T1[[i]])%*%g1 )
    r2=as.vector( M_func(T2[[i]])%*%g2 )
    
    l=l+sum( d1[[i]]*(log(r1)+beta1*Z1[[i]]) )+sum( d2[[i]]*(log(r2)+beta2*Z2[[i]]) )
    
    m1=sum(d1[[i]])
    m2=sum(d2[[i]])
    m12=sum(d1[[i]]*d2[[i]])
    
    func1=function(u){
      R1=as.vector( I_func(T1[[i]])%*%g1 )
      R2=as.vector( I_func(T2[[i]])%*%g2 )
      S1=pmin( exp( theta*u%*%t( exp(beta1*Z1[[i]])*R1 ) ), exp(500) )
      S2=pmin( exp( theta*u^alpha%*%t( exp(beta2*Z2[[i]])*R2 ) ), exp(500) )
      A=(S1+S2-1)
      Eta1=apply((S1/A)[,as.logical(d1[[i]]),drop=FALSE],MARGIN=1,FUN=prod)
      Eta2=apply((S2/A)[,as.logical(d2[[i]]),drop=FALSE],MARGIN=1,FUN=prod)
      Psi=rowSums( (1/theta)*log(A) )
      D12=exp(-Psi+500)  ### Adjustment to avoid too small D12 ###
      u^(m1+alpha*m2)*Eta1*Eta2*D12*(1+theta)^m12*dgamma(u,shape=1/eta,scale=eta)
      
    }
    
    Int=try( integrate(func1,0.001,10,stop.on.error = FALSE) ) 
    if( class(Int)=="try-error" ){l=l-500000}else{
      if(Int$value==0){l=l-500000}else{
        l=l+log(Int$value)-500 ### Re-adjustment to avoid too small D12 ### 
        
      }
    }
    
  }
  
  -l  
}

repeat{
  res=nlm(l.func,p=p0,hessian=TRUE)
  H_PL=res$hessian
  if( min( eigen(H_PL)$values )>0 ){break}
  p0=runif(14,-1,1)
}
H_PL=res$hessian

temp=(det(H_PL)==0)|is.na(det(H_PL))
if(temp){V=solve( H_PL+diag(rep(0.0001,14)) )}else{V=solve(H_PL)}

conv=res$code
iteration=res$iterations
ML=-res$minimum

beta1_est=res$est[1]
beta2_est=res$est[2]
g_est=exp(res$est[c(3,4,5,6,7)])
h_est=exp(res$est[c(8,9,10,11,12)])
eta_est=exp(res$est[13])
theta_est=exp(res$est[14])
tau_est=theta_est/(theta_est+2)

beta1_se=sqrt(V[1,1])
beta2_se=sqrt(V[2,2])
eta_se=exp(res$est[13])*sqrt(V[13,13])
theta_se=exp(res$est[14])*sqrt(V[14,14])
tau_se=2/((theta_est+2)^2)*theta_se

g_var=diag(exp(res$est[c(3,4,5,6,7)]))%*%
  V[c(3,4,5,6,7),c(3,4,5,6,7)]%*%diag(exp(res$est[c(3,4,5,6,7)]))
h_var=diag(exp(res$est[c(8,9,10,11,12)]))%*%
  V[c(8,9,10,11,12),c(8,9,10,11,12)]%*%diag(exp(res$est[c(8,9,10,11,12)]))

beta1_res=c(estimate=beta1_est,SE=beta1_se,
  Low95=beta1_est-1.96*beta1_se,Up95=beta1_est+1.96*beta1_se)
beta2_res=c(estimate=beta2_est,SE=beta2_se,
            Low95=beta2_est-1.96*beta2_se,Up95=beta2_est+1.96*beta2_se)
eta_res=c(estimate=eta_est,SE=eta_se,
          Low95=eta_est*exp(-1.96*sqrt(V[13,13])),Up95=eta_est*exp(1.96*sqrt(V[13,13])))
theta_res=c(estimate=theta_est,SE=theta_se,
            Low95=theta_est*exp(-1.96*sqrt(V[14,14])),Up95=theta_est*exp(1.96*sqrt(V[14,14])))
tau_res=c(estimate=tau_est,tau_se=tau_se,
          Low95=tau_est-1.95*tau_se,Up95=tau_est+1.95*tau_se)

list(beta1=beta1_res,beta2=beta2_res,eta=eta_res,theta=theta_res,
     tau=tau_res,K1=K1_est,K2=K2_est,
     g=g_est,h=h_est,g_var=g_var,h_var=h_var,
     convergence_code=conv,iterations=iteration,gradient=res$gradient,ML=ML)

}
