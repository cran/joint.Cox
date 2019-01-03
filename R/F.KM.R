F.KM <-
function(time,widths,t.death,death){
m=length(widths)

res_KM=survfit(Surv(t.death, death)~1)
t.sort=summary(res_KM)$time
res_SKM=summary(res_KM)$surv
res_SKM[max(sum(t.sort<=time),1)]
F_widths=numeric(m)
S_LM=res_SKM[max(sum(t.sort<=time),1)]
for(i in 1:m){ 
  F_widths[i]=1-res_SKM[max(sum(t.sort<=time+widths[i]),1)]/S_LM
}

cbind(t=time,w=widths,F=F_widths)
}
  