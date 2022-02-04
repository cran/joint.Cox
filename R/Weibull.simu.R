Weibull.simu <-
function(G,N,scale1,scale2,shape1=1,shape2=1,
        beta1,beta2,eta=0.5,copula="Clayton",theta=2,d=0,
        alpha=1,beta12=0,C.max=5,cmprsk=FALSE,tau=FALSE,
        Z.dist=runif,...){

  if((copula!="Clayton")&(copula!="Gumbel")&(copula!="Frank")&(copula!="BB1")){
    warning("The input for copula is wrong. It should be Clayton, Frank, Gumbel, or BB1")
  }else{

  X.vec=D.vec=C.vec=t.event=t.death=event=death=Z=group=tau.vec=NULL

  ij=0

  for(i in 1:G){

    u=rgamma(1,shape=1/eta,scale=eta)
    for(j in 1:N){

      ij=ij+1

      group[ij]=i
      Z[ij]=Z.dist(1,...)
      r1=scale1*u*exp(beta1*Z[ij])
      r2=scale2*(u^alpha)*exp(beta2*Z[ij])

      if(copula=="Clayton"){
        a=theta*exp(beta12*Z[ij])
        V1=runif(1)
        V2=runif(1)
        X=( -1/r1*log(1-V1) )^(1/shape1)
        W=(1-V1)^(-a)
        D=( 1/a/r2*log(1-W+W*(1-V2)^(-a/(a+1))) )^(1/shape2)
        Tau=a/(a+2)
      }

      if(copula=="Frank"){
        U=runif(1)
        W=runif(1)
        a=theta+beta12*Z[ij]
        A=exp(-a*U)-W*exp(-a*U)+W*exp(-a)
        B=exp(-a*U)-W*exp(-a*U)+W
        V=-log(A/B)/a
        X=(-log(U)/r1)^(1/shape1)
        D=(-log(V)/r2)^(1/shape2)
        func1=function(x){x/(exp(x)-1)}
        Tau=1-4/a*(1-integrate(func1,0,a)$value/a)
      }

      if(copula=="Gumbel"){
        U=runif(1)
        W=runif(1)
        a=theta*exp(beta12*Z[ij])
        func2=function(v){
          A=(-log(U))^a/U*((-log(U))^(a+1)+(-log(v))^(a+1))^(-a/(a+1))
          A*exp(-((-log(U))^(a+1)+(-log(v))^(a+1))^(1/(a+1)))-W
        }
        V=uniroot(func2,lower=0.00000000001,upper=0.9999999999)$root
        X=(-log(U)/r1)^(1/shape1)
        D=(-log(V)/r2)^(1/shape2)
        Tau=a/(a+1)
      }

      if(copula=="BB1"){
        U=runif(1)
        W=runif(1)
        a=theta*exp(beta12*Z[ij])
        U1=U^(-a)-1
        func3=function(v){
          v1=v^(-a)-1
          A=(U1^(d+1)+v1^(d+1))^(1/(d+1))
          U^(-a-1)*U1^d*A^(-d)*(1+A)^(-(1+a)/a)-W
        }
        V=uniroot(func3,lower=0.00000000001,upper=0.9999999999)$root
        X=(-log(U)/r1)^(1/shape1)
        D=(-log(V)/r2)^(1/shape2)
        Tau=1-2/(d+1)/(a+2)
      }

      C=runif(1,min=0,max=C.max)
      X.vec[ij]=X
      D.vec[ij]=D
      C.vec[ij]=C
      t.event[ij]=min(X,D,C)
      t.death[ij]=min(D,C)
      event[ij]=as.numeric( t.event[ij]==X )
      death[ij]=as.numeric( t.death[ij]==D )
      tau.vec[ij]=Tau
    }

  }

  Dat=data.frame(X=X.vec,D=D.vec,C=C.vec,
             t.event=t.event,event=event,
             t.death=t.death,death=death,
             group=group,Z=Z)

  if(cmprsk==TRUE){
    event1=as.numeric(t.event==X.vec)
    event2=as.numeric(t.event==D.vec)
    Dat=data.frame(X=X.vec,D=D.vec,C=C.vec,
                   t.event=t.event,event1=event1,
                   event2=event2,group=group,Z=Z)
  }

  if(tau==TRUE){ Dat=cbind(Dat,tau=tau.vec) }
  Dat
  }
}

