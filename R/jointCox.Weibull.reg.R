jointCox.Weibull.reg <- 
function (t.event,event,t.death,death,Z1,Z2,group,alpha=1,
          Randomize_num=10,Adj=500,convergence.par=FALSE){

  T1 = t.event
  T2 = t.death
  d1 = event
  d2 = death
  Z1 = as.matrix(Z1)
  Z2 = as.matrix(Z2)
  p1 = ncol(Z1)
  p2 = ncol(Z2)
  
  G_id = as.numeric((levels(factor(group))))
  G = length(G_id)
  
  N = length(t.event)
  ########### Summary ###########
  n.event=xtabs(d1~group)
  n.death=xtabs(d2~group)
  n.censor=xtabs(1-d2~group)
  count=cbind(table(group),n.event,n.death,n.censor)  
  colnames(count)=c("No.of samples","No.of events","No.of deaths","No.of censors")
 
  T1_mean = mean(T1)
  T2_mean = mean(T2)

  ## Likelihood function ##
  l.func = function(phi) {
    g1 = exp(pmax(pmin(phi[1:2], 500), -500))
    g2 = exp(pmax(pmin(phi[3:4], 500), -500))
    eta = exp(phi[5])
    theta = min(exp(phi[6]), exp(3))
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    l = 0
    bZ1 = as.vector(Z1 %*% beta1)
    bZ2 = as.vector(Z2 %*% beta2)

    ## Weibull baseline hazard ##
    weibull.h = function(x, scale0, shape0) {
      scale0 * shape0 * x^(shape0 - 1)
    }
    weibull.H = function(x, scale0, shape0) {
      scale0 * x^(shape0)
    }

    r1 = as.vector(weibull.h(T1, scale0 = g1[1], shape0 = g1[2]))
    R1 = as.vector(weibull.H(T1, scale0 = g1[1], shape0 = g1[2]))
    r2 = as.vector(weibull.h(T2, scale0 = g2[1], shape0 = g2[2]))
    R2 = as.vector(weibull.H(T2, scale0 = g2[1], shape0 = g2[2]))
    l = l + sum(d1 * (log(r1) + bZ1)) + sum(d2 * (log(r2) + bZ2))
    for (i in G_id) {
      Gi = c(group == i)
      m1 = sum(d1[Gi])
      m2 = sum(d2[Gi])
      m12 = sum(d1[Gi] * d2[Gi])
      EZ1 = exp(bZ1[Gi]) * R1[Gi]
      EZ2 = exp(bZ2[Gi]) * R2[Gi]
      D1 = as.logical(d1[Gi])
      D2 = as.logical(d2[Gi])

      ## integration , Adj to avoid too small value ##
      func1 = function(u) {
        S1 = pmin(exp(theta * u %*% t(EZ1)), exp(500))
        S2 = pmin(exp(theta * u^alpha %*% t(EZ2)), exp(500))
        A = (S1 + S2 - 1)
        E1 = apply(log(S1/A)[, D1, drop = FALSE], MARGIN = 1, FUN = sum)
        E2 = apply(log(S2/A)[, D2, drop = FALSE], MARGIN = 1, FUN = sum)
        Psi = rowSums((1/theta) * log(A))
        exp((m1 + alpha * m2) * log(u) + E1 + E2 - Psi + m12 * log(1 + theta) + log(dgamma(u, shape = 1/eta, scale = eta)) + Adj)
      }
      Int = try(integrate(func1, 0.001, 10, stop.on.error = FALSE))
      if (class(Int) == "try-error") {
        l = l - 5e+05
      }
      else {
        if (Int$value == 0) {
          l = l - 5e+05
        }
        else {
          l = l + log(Int$value) - Adj
        }
      }
    }
    -l
  }
  p0 = rep(0, 6 + p1 + p2)
  p0[c(1, 3)] = p0[c(1, 3)] - log(c(T1_mean, T2_mean))
  res = nlm(l.func, p = p0, hessian = TRUE)
  MPL = -res$minimum
  R_num = 0
  repeat {
    if ((min(eigen(res$hessian)$values) > 0) & (res$code == 1)) {
      break
    }
    if (R_num >= Randomize_num) {
      break
    }
    R_num = R_num + 1
    p0_Rand = runif(6 + p1 + p2, -1, 1)
    p0_Rand[c(1, 3)] = p0_Rand[c(1, 3)] - log(c(T1_mean, T2_mean))
    res_Rand = nlm(l.func, p = p0_Rand, hessian = TRUE)
    MPL_Rand = -res_Rand$minimum
    if (MPL_Rand > MPL) {
      res = res_Rand
      MPL = -res$minimum
    }
  }
  H_L = -res$hessian
  V = solve(-H_L, tol = 10^(-50))

  DF=6 + p1 + p2
  AIC = 2 * DF - 2 * (-l.func(res$estimate))
  BIC = log(N) * DF - 2 * (-l.func(res$estimate))
  convergence_res = c(MPL = MPL, DF = DF, AIC = AIC, BIC = BIC,
                      code = res$code, No.of.iterations = res$iterations, No.of.randomizations = R_num)
  est = c(exp(res$est[1:6]), res$est[(6 + 1):(6 + p1 + p2)])
  est_var = diag(c(est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(est[1:6], rep(1, p1 + p2)))

  beta1_est = res$est[(6 + 1):(6 + p1)]
  beta2_est = res$est[(6 + p1 + 1):(6 + p1 + p2)]
  g_est = exp(res$est[1:2])
  h_est = exp(res$est[3:4])
  eta_est = exp(res$est[5])
  theta_est = exp(res$est[6])
  tau_est = theta_est/(theta_est + 2)

  beta1_se = sqrt(diag(V)[(6 + 1):(6 + p1)])
  beta2_se = sqrt(diag(V)[(6 + p1 + 1):(6 + p1 + p2)])
  eta_se = eta_est * sqrt(diag(V)[5])
  theta_se = theta_est * sqrt(diag(V)[6])
  tau_se = 2/(theta_est + 2)^2 * theta_se
  g_var = diag(g_est) %*% V[1:2, 1:2] %*% diag(g_est)
  h_var = diag(h_est) %*% V[3:4, 3:4] %*% diag(h_est)
  g_se = sqrt(diag(g_var))
  h_se = sqrt(diag(h_var))

  
  beta1_res=c(estimate=beta1_est,SE=beta1_se,
              Lower=beta1_est-1.96*beta1_se,Upper=beta1_est+1.96*beta1_se)
  beta2_res=c(estimate=beta2_est,SE=beta2_se,
              Lower=beta2_est-1.96*beta2_se,Upper=beta2_est+1.96*beta2_se)
  eta_res = c(Estimate=eta_est,SE=eta_se,
              Lower=eta_est*exp(-1.96*sqrt(diag(V)[5])),
              Upper=eta_est*exp(1.96*sqrt(diag(V)[5])))
  theta_Lower=theta_est*exp(-1.96*sqrt(diag(V)[6]))
  theta_Upper=theta_est*exp(1.96*sqrt(diag(V)[6]))
  theta_res=c(Estimate=theta_est,SE=theta_se,Lower=theta_Lower,Upper=theta_Upper)
  tau_res=c(Estimate=tau_est,SE=tau_se,
            Lower=theta_Lower/(theta_Lower+2),Upper=theta_Upper/(theta_Upper+2))
  scale1=c(Estimate=g_est[1],SE=g_se[1],Lower=g_est[1]*exp(-1.96*sqrt(diag(V)[1:1])),
           Upper=g_est[1]*exp(+1.96*sqrt(diag(V)[1:1])))
  shape1=c(Estimate=g_est[2],SE=g_se[2],Lower=g_est[2]*exp(-1.96*sqrt(diag(V)[2:2])),
           Upper=g_est[2]*exp(+1.96*sqrt(diag(V)[2:2])))
  
  scale2=c(Estimate=h_est[1],SE=h_se[1],Lower=h_est[1]*exp(-1.96*sqrt(diag(V)[3:3])),
           Upper=h_est[1]*exp(+1.96*sqrt(diag(V)[3:3])))
  shape2=c(Estimate=h_est[2],SE=h_se[2],Lower=h_est[2]*exp(-1.96*sqrt(diag(V)[4:4])),
           Upper=h_est[2]*exp(+1.96*sqrt(diag(V)[4:4])))
 
  if (convergence.par==FALSE){convergence.parameters = NULL}else{
    convergence.parameters=list(log_estimate=res$est,gradient=-res$gradient,log_var=V)
  }
  list(count=count,
       alpha=alpha,beta1=beta1_res,beta2=beta2_res,
       eta=eta_res,theta=theta_res,tau=tau_res,
       scale1=scale1,shape1=shape1,scale2=scale2,shape2=shape2,
       convergence=convergence_res,convergence.parameters=convergence.parameters)
}
