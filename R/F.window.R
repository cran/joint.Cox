F.window <-
function(time,width,X,Z1,Z2,beta1,beta2,eta,theta,alpha,
                     g,h,xi1,xi3,Fplot=TRUE){

if(time<xi1){warning("time should be larger than xi1")}
if(time+width>=xi3){warning("out-of-prediction bound; time+width should be smaller than xi3")}

########### M-spline matrix #############
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

g1=g
g2=h
  
bZ1=t(Z1)%*%beta1
bZ2=t(Z2)%*%beta2

if(time<xi1){R1_t=0}else{R1_t=as.vector( I_func(time)%*%g1 )}
if(X<xi1){R1_X=0}else{R1_X=as.vector( I_func(X)%*%g1 )}
if(time<xi1){R2_t=0}else{R2_t=as.vector( I_func(time)%*%g2 )}
if(time+width<xi1){R2_tw=0}else{R2_tw=as.vector( I_func(time+width)%*%g2 )}

S_t_u=function(u){
  S1_t=exp( theta*u%*%t( exp(bZ1)*R1_t ) )
  S2_t=exp( theta*u^alpha%*%t( exp(bZ2)*R2_t ) )
  (S1_t+S2_t-1)^(-1/theta)
}

S_tw_u=function(u){
  S1_t=exp( theta*u%*%t( exp(bZ1)*R1_t ) )
  S2_tw=exp( theta*u^alpha%*%t( exp(bZ2)*R2_tw ) )
  (S1_t+S2_tw-1)^(-1/theta)
}

func1=function(u){  pmax(S_t_u(u),0,na.rm=TRUE)*dgamma(u,shape=1/eta,scale=eta) }
S_t=integrate(func1,0.001,10,stop.on.error = FALSE)$value

func2=function(u){  pmax(S_tw_u(u),0,na.rm=TRUE)*dgamma(u,shape=1/eta,scale=eta) }
S_tw=integrate(func2,0.001,10,stop.on.error = FALSE)$value

func3=function(u){
  f3=( exp(-u^alpha%*%t(exp(bZ2)*R2_t))-S_t_u(u) )*dgamma(u,shape=1/eta,scale=eta) 
  pmax(f3,0,na.rm=TRUE)
}
S_0t=integrate(func3,0.001,10,stop.on.error = FALSE)$value

func4=function(u){
  f4=( exp(-u^alpha%*%t(exp(bZ2)*R2_tw))-S_tw_u(u) )*dgamma(u,shape=1/eta,scale=eta) 
  pmax(f4,0,na.rm=TRUE)
}
S_0tw=integrate(func4,0.001,10,stop.on.error = FALSE)$value

func5=function(u){
  S1_X=exp( theta*u%*%t( exp(bZ1)*R1_X ) )
  S2_t=exp( theta*u^alpha%*%t( exp(bZ2)*R2_t ) )
  f5=u*S1_X*(S1_X+S2_t-1)^(-1/theta-1)*dgamma(u,shape=1/eta,scale=eta) 
  pmax(f5,0,na.rm=TRUE)
}
S_Xt=integrate(func5,0.001,10,stop.on.error = FALSE)$value

func6=function(u){
  S1_X=exp( theta*u%*%t( exp(bZ1)*R1_X ) )
  S2_tw=exp( theta*u^alpha%*%t( exp(bZ2)*R2_tw ) )
  f6=u*S1_X*(S1_X+S2_tw-1)^(-1/theta-1)*dgamma(u,shape=1/eta,scale=eta) 
  pmax(f6,0,na.rm=TRUE)
}
S_Xtw=integrate(func6,0.001,10,stop.on.error = FALSE)$value

if((S_tw<=0)|(S_t<=0)){F_noevent=1}else{F_noevent=1-S_tw/S_t}
if((S_0tw<=0)|(S_0t<=0)){F_event=F_noevent}else{F_event=1-(S_0tw)/(S_0t)}
if((S_Xtw<=0)|(S_Xt<=0)){F_event_at_X=F_noevent}else{F_event_at_X=1-(S_Xtw)/(S_Xt)}

if(Fplot==TRUE){
  num_grid=500
  x_grid=seq(X,time,length=num_grid)
  plot(x_grid,rep(0,num_grid),xlim=c(xi1,xi3),ylim=c(0,1),type="l",lwd=4,
       xlab="t",ylab="Probability of death between t and t+w",col="red")
  abline(h=0)
  abline(v=time,col="gray")
  abline(v=time+width,col="gray")
  points(X,0,lwd=4,col="red")
  p_grid=seq(time,time+width,length=num_grid)
  points(p_grid,rep(F_event_at_X,num_grid),col="red")
  points(p_grid,rep(F_noevent,num_grid),col="blue")
  text(X,0.05,"X",cex=1)
  text(time,0.05,"t",cex=1)
  text(time+width,0.05,"t+w",cex=1)
  text(time,F_event_at_X+0.05,"Event occurred at X (< t)",cex=1)
  text(time,F_noevent+0.05,"Event not occurred before t",cex=1)
}

c(t=time,w=width,X=X,F_event_at_X=F_event_at_X,F_event=F_event,F_noevent=F_noevent)
}
