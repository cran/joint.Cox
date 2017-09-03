F.prediction <-
function(time,widths,X,Z1,Z2,beta1,beta2,eta,theta,alpha,
                     g,h,xi1,xi3,Fplot=TRUE){
temp="F_noevent"
if(is.numeric(X)){ if(X<=time){temp="F_event_at_X"} }
  
F.pred=F.windows(time,widths,X,Z1,Z2,beta1,beta2,
                 eta,theta,alpha,g,h,xi1,xi3,Fplot=FALSE)[,temp]

if(Fplot==TRUE){   
  plot(time+widths,F.pred,xlim=c(0, max(time+widths)),ylim=c(-0.05,1.05),
       type="l",lwd=2,xlab="t+w",ylab="Probability of death in ( t, t+w )",col="blue")
  abline(h=0)
  abline(v=time,col="gray")
  if(X<=time){  points(X,0,lwd=3,col="red") }
  text(time,-0.05,"t",cex=1)
}

cbind(t=time,w=widths,X=X,F=F.pred)
}
