Weibull.simu <-
function(G,N,scale1,scale2,shape1=1,shape2=1,
        beta1,beta2,eta=0.5,copula="Clayton",theta=2,
        alpha=1,beta12=0,C.max=5,
        Z.dist=runif,...){

  if((copula!="Clayton")&(copula!="Gumbel")){
    warning("The argument copula is wrong. It should either Clayton or Gumbel")
  }else{

  X.vec=D.vec=C.vec=t.event=t.death=event=death=Z=group=NULL

  ij=0

  for(i in 1:G){

    u=rgamma(1,shape=1/eta,scale=eta)
    for(j in 1:N){

      ij=ij+1

      group[ij]=i
      Z[ij]=Z.dist(1,...)
      r1=scale1*u*exp(beta1*Z[ij])
      r2=scale2*(u^alpha)*exp(beta2*Z[ij])
      theta12=theta*exp(beta12*Z[ij])

      if(copula=="Clayton"){
        V1=runif(1)
        V2=runif(1)
        X=( -1/r1*log(1-V1) )^(1/shape1)
        W=(1-V1)^(-theta12)
        D=( 1/theta12/r2*log(1-W+W*(1-V2)^(-theta12/(theta12+1))) )^(1/shape2)
      }

      if(copula=="Gumbel"){
        U=runif(1)
        W=runif(1)
        a=theta12
        func=function(v){
          A=(-log(U))^a/U*((-log(U))^(a+1)+(-log(v))^(a+1))^(-a/(a+1))
          A*exp(-((-log(U))^(a+1)+(-log(v))^(a+1))^(1/(a+1)))-W
        }
        V=uniroot(func,lower=0,upper=1)$root
        X=(-log(U)/r1)^(1/shape1)
        D=(-log(V)/r2)^(1/shape2)
      }

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
}

