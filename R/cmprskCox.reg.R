cmprskCox.reg = function (t.event,event1,event2,Z1,Z2,group,theta,alpha = 1,
                          kappa_grid = c(seq(10,1e+17,length = 30)),LCV_plot = TRUE,
                          Randomize_num = 10,Adj = 500,convergence.par = FALSE) {
  
  d1 = event1
  d2 = event2
  p1 = ncol(Z1)
  p2 = ncol(Z2)
  G_id = as.numeric((levels(factor(group))))
  G = length(G_id)
  
  n.event1 = tapply(d1,group,FUN = sum)
  n.event2 = tapply(d2,group,FUN = sum)
  n.censor = tapply(1-d1-d2,group,FUN = sum)
  count = cbind(table(group),n.event1,n.event2,n.censor)
  colnames(count) = c("No.of samples","No.of events","No.of deaths","No.of censors")
  
  xi1 = min(t.event)
  xi3 = max(t.event)
  res1 = splineCox.reg(t.event,event1,Z1,xi1 = xi1,xi3 = xi3, 
                       kappa_grid = kappa_grid,LCV_plot = LCV_plot)
  res2 = splineCox.reg(t.event,event2,Z2,xi1 = xi1,xi3 = xi3, 
                       kappa_grid = kappa_grid,LCV_plot = LCV_plot)
  
  K1_est = res1$kappa
  LCV1_res = c(K1 = K1_est,LCV1 = res1$LCV)
  K2_est = res2$kappa
  LCV2_res = c(K2 = K2_est,LCV2 = res2$LCV)
  
  Omega = c( 192,-132,  24,  12,   0,
            -132,  96, -24, -12,  12,
              24, -24,  24, -24,  24,
              12, -12, -24,  96,-132,
               0,  12,  24,-132, 192)
  Omega = matrix(Omega,5,5)/((xi3-xi1)/2)^5
  
  l.func = function(phi) {
    
    g1 = exp(pmax(pmin(phi[1:5],500),-500))
    g2 = exp(pmax(pmin(phi[6:10],500),-500))
    eta = exp(phi[11])
    beta1 = phi[(11+1):(11+p1)]
    beta2 = phi[(11+p1+1):(11+p1+p2)]
    
    l = -K1_est*t(g1)%*%Omega%*%g1-K2_est*t(g2)%*%Omega%*%g2
    
    bZ1 = as.vector(as.matrix(Z1)%*%beta1)
    bZ2 = as.vector(as.matrix(Z2)%*%beta2)
    r1 = as.vector(M.spline(t.event,xi1 = xi1,xi3 = xi3)%*%g1)
    r2 = as.vector(M.spline(t.event,xi1 = xi1,xi3 = xi3)%*%g2)
    R1 = as.vector(I.spline(t.event,xi1 = xi1,xi3 = xi3)%*%g1)
    R2 = as.vector(I.spline(t.event,xi1 = xi1,xi3 = xi3)%*%g2)
    
    l = l+sum(d1*(log(r1)+bZ1))+sum(d2*(log(r2)+bZ2))
    
    for (i in G_id) {
      
      Gi = c(group == i)
      m1 = sum(d1[Gi])
      m2 = sum(d2[Gi])
      EZ1 = exp(bZ1[Gi])*R1[Gi]
      EZ2 = exp(bZ2[Gi])*R2[Gi]
      D1 = as.logical(d1[Gi])
      D2 = as.logical(d2[Gi])
      
      func1 = function(u) {
        
        S1 = pmin(exp(theta*u%*%t(EZ1)),exp(500))
        S2 = pmin(exp(theta*u^alpha %*% t(EZ2)),exp(500))
        A = (S1+S2-1)
        E1 = apply((S1/A)[,D1,drop = FALSE],MARGIN = 1,FUN = prod)
        E2 = apply((S2/A)[,D2,drop = FALSE],MARGIN = 1,FUN = prod)
        Psi = rowSums((1/theta)*log(A))
        D12 = exp(-Psi+Adj)
        
        return(u^(m1+alpha*m2)*E1*E2*D12*dgamma(u,shape = 1/eta,scale = eta))
        
      }
      
      Int = try(integrate(func1,0.001,10,stop.on.error = FALSE))
      
      if (class(Int) == "try-error") {
        
        l = l-5e+05
        
      } else {
        
        if (Int$value == 0) {
          
          l = l-5e+05
          
        } else {
          
          l = l+log(Int$value)-Adj
          
        }
        
      }
      
    }
    
    return(-l)
    
  }
  
  p0 = rep(0,11+p1+p2)
  res = nlm(l.func,p = p0,hessian = TRUE)
  MPL = -res$minimum
  R_num = 0
  
  repeat {
    
    if (min(eigen(res$hessian)$values) > 0 & res$code == 1) {break}
    if (R_num >= Randomize_num) {break}
    
    R_num = R_num + 1
    p0_Rand = runif(11+p1+p2,-1,1)
    res_Rand = nlm(l.func,p = p0_Rand,hessian = TRUE)
    MPL_Rand = -res_Rand$minimum
    
    if (MPL_Rand > MPL) {
      
      res = res_Rand
      MPL = -res$minimum
      
    }
    
  }
  
  H_PL = -res$hessian
  DF_upper = 18+p1+p2
  temp = (det(H_PL) == 0) | is.na(det(H_PL))
  
  if (temp) {
    V = solve(-H_PL+diag(rep(1e-04,11+p1+p2)),tol = 1e-50)
  } else {V = solve(-H_PL,tol = 1e-50)}
  
  D_PL = diag(c(1/exp(res$estimate[1:11]),rep(1,p1 +p2)))
  H_PL = D_PL%*%H_PL%*%D_PL
  H = H_PL
  H[1:5,1:5] = H[1:5,1:5]+2*K1_est*Omega
  H[6:10,6:10] = H[6:10,6:10]+2*K2_est*Omega
  
  if (is.na(det(H_PL)) | det(H_PL) == 0) {  DF = DF_upper   } else {
    DF = min(max(sum(diag(solve(H_PL,tol = 1e-50)%*%H)),p1+p2+2),DF_upper)
  }
  
  K1_est = K2_est = 0
  LCV = -l.func(res$estimate)-DF
  convergence_res = c(MPL = MPL,DF = DF,LCV = LCV,code = res$code,
                      No.of.iterations = res$iterations,No.of.randomizations = R_num)
  beta1_est = res$est[(11+1):(11+p1)]
  beta2_est = res$est[(11+p1+1):(11+p1+p2)]
  g_est = exp(res$est[1:5])
  h_est = exp(res$est[6:10])
  eta_est = exp(res$est[11])
  beta1_se = sqrt(diag(V)[(11+1):(11+p1)])
  beta2_se = sqrt(diag(V)[(11+p1+1):(11+p1+p2)])
  eta_se = eta_est*sqrt(diag(V)[11])
  g_var = diag(g_est)%*%V[1:5,1:5]%*%diag(g_est)
  h_var = diag(h_est)%*%V[6:10,6:10]%*%diag(h_est)
  beta1_res = c(estimate = beta1_est,SE = beta1_se,
                Lower = beta1_est-1.96*beta1_se,
                Upper = beta1_est+1.96*beta1_se)
  beta2_res = c(estimate = beta2_est,SE = beta2_se,
                Lower = beta2_est-1.96*beta2_se,
                Upper = beta2_est+1.96*beta2_se)
  eta_res = c(estimate = eta_est,SE = eta_se,
              Lower = eta_est*exp(-1.96*sqrt(diag(V)[11])),
              Upper = eta_est*exp( 1.96*sqrt(diag(V)[11])))
  
  if (convergence.par == FALSE) {    convergence.parameters = NULL  } else {
    convergence.parameters = list(log_estimate = res$est,gradient = -res$gradient,log_var = V)
  }
  
  return(list(count = count,beta1 = beta1_res,beta2 = beta2_res, 
         eta = eta_res,theta=theta,tau=theta/(theta+2),LCV1 = LCV1_res,LCV2 = LCV2_res,
         g = g_est,h = h_est,g_var = g_var,h_var = h_var,
         convergence = convergence_res,convergence.parameters = convergence.parameters))
  
}