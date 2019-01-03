Weibull.simu <-
function(G,N,scale1,scale2,shape1=1,shape2=1,
        beta1,beta2,eta=0.5,theta=2,alpha=1,C.max=5){
  
  X.vec=D.vec=C.vec=t.event=t.death=event=death=Z=group=NULL
  
  ij=0
  
  for(i in 1:G){
    
    u=rgamma(1,shape=1/eta,scale=eta)
    for(j in 1:N){
      
      ij=ij+1
      
      group[ij]=i
      Z[ij]=runif(1)
      r1=scale1*u*exp(beta1*Z[ij])
      r2=scale2*(u^alpha)*exp(beta2*Z[ij])
      V1=runif(1)
      V2=runif(1)
      X=( -1/r1*log(1-V1) )^(1/shape1)
      W=(1-V1)^(-theta)
      D=( 1/theta/r2*log(1-W+W*(1-V2)^(-theta/(theta+1))) )^(1/shape2)
      C=runif(1,min=0,max=C.max)
      X.vec[ij]=X
      D.vec[ij]=D
      C.vec[ij]=C
      t.event[ij]=min(X,D,C)
      t.death[ij]=min(D,C)
      event[ij]=as.numeric( t.event[ij]==X )
      death[ij]=as.numeric( t.death[ij]==D )
    }
    
  }
 
  data.frame(X=X.vec,D=D.vec,C=C.vec,
             t.event=t.event,event=event,
             t.death=t.death,death=death,
             group=group,Z=Z)
}

