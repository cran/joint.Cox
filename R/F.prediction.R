F.prediction <-
function(time,widths,X,Z1,Z2,beta1,beta2,eta,theta,alpha,
                     g,h,xi1,xi3,Fplot=TRUE){
temp="F_noevent"
if(is.numeric(X)){ if(X<=time){temp="F_event_at_X"} }
  
F.pred=F.windows(time,widths,X,Z1,Z2,beta1,beta2,
                 eta,theta,alpha,g,h,xi1,xi3,Fplot=FALSE)[,temp]

if(Fplot==TRUE){   print("draw plot in future update")  }
cbind(t=time,w=widths,X=X,F=F.pred)
}
