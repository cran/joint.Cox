F.window.Weibull <-function (time, width, X, Z1, Z2, beta1, beta2, eta, theta, alpha, 
                              scale1, shape1,scale2,shape2, xi1, xi3, Fplot = TRUE){
    if (theta < 0) {
      warning("in general, theta should be positive or zero")
    }
    
    bZ1 = t(Z1) %*% beta1
    bZ2 = t(Z2) %*% beta2
    
    R1_t=scale1*time^shape1
    R1_X = scale1*X^shape1
    R2_t=scale2*time^shape2
    R2_tw= scale2*(time+width)^shape2
    
    S1_t_func = function(u) {
      exp(-u %*% t(exp(bZ1) * R1_t))
    }
    S1_X_func = function(u) {
      exp(-u %*% t(exp(bZ1) * R1_X))
    }
    S2_t_func = function(u) {
      exp(-u^alpha %*% t(exp(bZ2) * R2_t))
    }
    S2_tw_func = function(u) {
      exp(-u^alpha %*% t(exp(bZ2) * R2_tw))
    }
    if (theta > 0) {
      S_t_func = function(u) {
        A = pmax(S1_t_func(u)^(-theta) + S2_t_func(u)^(-theta) - 
                   1, 0)
        A^(-1/theta) * dgamma(u, shape = 1/eta, scale = eta)
      }
      S_t = integrate(S_t_func, 0.001, 10, stop.on.error = FALSE)$value
      S_tw_func = function(u) {
        A = pmax(S1_t_func(u)^(-theta) + S2_tw_func(u)^(-theta) - 
                   1, 0)
        A^(-1/theta) * dgamma(u, shape = 1/eta, scale = eta)
      }
      S_tw = integrate(S_tw_func, 0.001, 10, stop.on.error = FALSE)$value
      S_Xt_func = function(u) {
        A = pmax(S1_X_func(u)^(-theta) + S2_t_func(u)^(-theta) - 
                   1, 0)
        B = u * S1_X_func(u)^(-theta) * A^(-1/theta - 1) * 
          dgamma(u, shape = 1/eta, scale = eta)
        pmax(B, 0, na.rm = TRUE)
      }
      S_Xt = integrate(S_Xt_func, 0.001, 10, stop.on.error = FALSE)$value
      S_Xtw_func = function(u) {
        A = pmax(S1_X_func(u)^(-theta) + S2_tw_func(u)^(-theta) - 
                   1, 0)
        B = u * S1_X_func(u)^(-theta) * A^(-1/theta - 1) * 
          dgamma(u, shape = 1/eta, scale = eta)
        pmax(B, 0, na.rm = TRUE)
      }
      S_Xtw = integrate(S_Xtw_func, 0.001, 10, stop.on.error = FALSE)$value
    }
    else {
      S_t_func = function(u) {
        S1_t_func(u) * S2_t_func(u) * dgamma(u, shape = 1/eta, 
                                             scale = eta)
      }
      S_t = integrate(S_t_func, 0.001, 10, stop.on.error = FALSE)$value
      S_tw_func = function(u) {
        S1_t_func(u) * S2_tw_func(u) * dgamma(u, shape = 1/eta, 
                                              scale = eta)
      }
      S_tw = integrate(S_tw_func, 0.001, 10, stop.on.error = FALSE)$value
      S_Xt_func = function(u) {
        u * S1_X_func(u) * S2_t_func(u) * dgamma(u, shape = 1/eta, 
                                                 scale = eta)
      }
      S_Xt = integrate(S_Xt_func, 0.001, 10, stop.on.error = FALSE)$value
      S_Xtw_func = function(u) {
        u * S1_X_func(u) * S2_tw_func(u) * dgamma(u, shape = 1/eta, 
                                                  scale = eta)
      }
      S_Xtw = integrate(S_Xtw_func, 0.001, 10, stop.on.error = FALSE)$value
    }
    if ((S_tw <= 0) | (S_t <= 0)) {
      F_noevent = 1
    }
    else {
      F_noevent = 1 - S_tw/S_t
    }
    if ((S_Xtw <= 0) | (S_Xt <= 0)) {
      F_event_at_X = F_noevent
    }
    else {
      F_event_at_X = 1 - (S_Xtw)/(S_Xt)
    }
    if (Fplot == TRUE) {
      num_grid = 500
      x_grid = seq(X, time, length = num_grid)
      plot(x_grid, rep(0, num_grid), xlim = c(xi1, xi3), ylim = c(-0.05, 
                                                                  1.05), type = "l", lwd = 3, xlab = "time", ylab = "Probability of death in ( t, t+w )", 
           col = "red")
      polygon(c(time, time + width, time + width, time), c(0, 
                                                           0, 1.092, 1.092), col = gray(0.95))
      abline(h = 0)
      abline(v = time, col = "gray")
      abline(v = time + width, col = "gray")
      points(X, 0, lwd = 5, col = "red")
      p_grid = seq(time, time + width, length = num_grid)
      points(p_grid, rep(F_event_at_X, num_grid), col = "red")
      points(p_grid, rep(F_noevent, num_grid), col = "blue")
      text(X, 0.05, "X", cex = 1)
      text(time, -0.05, "t", cex = 1)
      text(time + width, -0.05, "t+w", cex = 1)
      text(time, F_event_at_X + 0.03, "Event at X (< t)", cex = 1)
      text(time, F_noevent + 0.03, "No event (X>t)", cex = 1)
    }
    c(t = time, w = width, X = X, F_event_at_X = F_event_at_X, 
      F_noevent = F_noevent)
}
