splineCox.reg <-
function(t.event,event,Z,xi1=min( t.event ),xi3=max( t.event ),
         kappa=c(seq(10,1e+17,length=30)),LCV.plot=TRUE,p0=rep(0,5+p)){

d=event
Z=as.matrix(Z)
p=ncol(Z)

Omega=c(192,-132,24,12,0,
        -132,96,-24,-12,12,
        24,-24,24,-24,24,
        12,-12,-24,96,-132,
        0,12,24,-132,192)
Omega=matrix(Omega,5,5)/( (xi3-xi1)/2 )^5

## Log-likelihood function ##
l.func=function(phi){
  b=phi[(5+1):(5+p)]
  g=exp(  pmin(phi[1:5],500)  ) ## M-spline coefficients ##
  l=-K*t(g)%*%Omega%*%g
  r=as.vector( M.spline(t.event,xi1=xi1,xi3=xi3)%*%g )
  R=as.vector( I.spline(t.event,xi1=xi1,xi3=xi3)%*%g )
  bZ=as.numeric( Z%*%b )
  l=l+sum( d*(log(r)+bZ) )
  l=l-sum(  pmin( exp(bZ)*R, exp(500) )  )
  -l  
}

DF_upper=18+p

L=DF=NULL

for(k in 1:length(kappa)){
  K=kappa[k]
  res=nlm(l.func,p=p0,hessian=TRUE)
  D_PL=diag( c(1/exp(res$estimate[1:5]),rep(1,p)) )
  H_PL=-D_PL%*%res$hessian%*%D_PL
  H=H_PL
  H[1:5,1:5]=H[1:5,1:5]+2*K*Omega
  K=0
  L[k]=-l.func(res$estimate) 
  if( is.na(det(H_PL))|det(H_PL)==0 ){DF[k]=DF_upper}else{
    DF[k]=min( max( sum( diag(solve(H_PL,tol=10^(-50))%*%H) ), p+2) ,DF_upper)
  }
}

LCV=L-DF
K=K_est=kappa[which.max(LCV)]
DF_est=DF[which.max(LCV)]
LCV_est=LCV[which.max(LCV)]

########## Plotting LCV ##########
if(LCV.plot==TRUE){
  par(mfrow=c(1,3))
  plot(kappa,L,xlab="K",ylab="logL",type="b",lwd=3)
  plot(kappa,pmin(DF,10+p),xlab="K",ylab="DF",type="b",lwd=3)
  plot(kappa,LCV,xlab="K",ylab="LCV=logL-DF",type="b",lwd=3)
  points(K_est,LCV_est,col="red",pch=17,cex=2)
}

res=nlm(l.func,p=p0,hessian=TRUE)

beta_est=res$est[(5+1):(5+p)]
h_est=exp(res$est[1:5])

H_PL=-res$hessian
V=solve(-H_PL,tol=10^(-50))
beta_se=sqrt(diag(V)[(5+1):(5+p)])
h_var=diag(h_est)%*%V[1:5,1:5]%*%diag(h_est)

beta_res=c(estimate=beta_est,SE=beta_se,
        Lower=beta_est-1.96*beta_se,Upper=beta_est+1.96*beta_se)

list(beta=beta_res,h=h_est,h_var=h_var,
     kappa=K_est,DF=DF_est,LCV=LCV_est)

}
