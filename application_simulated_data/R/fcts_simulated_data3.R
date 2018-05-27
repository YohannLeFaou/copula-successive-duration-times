
noyau_quadratique = function(x){return(15/16 * (1-x^2)^2 * ( abs(x) <= 1) )}

sigmoid = function(x){exp(x)/(1 + exp(x))}

calc_rate_cens = function(lambda, x){
  sum( exp(rexp(length(x), rate = lambda)) < x ) / length(x)
}

obj_fct_exp = function(lambda, x, target){
  return(abs( calc_rate_cens(lambda = lambda, x = x) - target) )
}


loglike_copula = function(u1, u2, family, par, poids){
  - sum( poids * log( BiCopPDF(u1 = u1, u2 = u2, family = family, par = par)))
}

make_KM_weights = function(vect_y, vect_delta, vect_c = NULL){
  # function which computes the weights of the observations in the Kaplan Meier case
  ## vect_c : values taken by C in the case C is never censored
  
  is_C_censored = is.null(vect_c)
  # is_C_censored = FALSE if C is not censored (in some situations, we allways observe C)
  # is_C_censored = TRUE if C is censored (by the time of interest)
  
  if(is_C_censored){
    # S_C is estimated by Kaplan Meier estimator
    delta_c = 1 - vect_delta
    data_c = data.frame(y = vect_y, delta_c = delta_c)
    km_c = survival::survfit(survival::Surv(time = y, event = delta_c) ~ 1,
                             data = data_c,
                             type = "kaplan-meier")
  } else{
    # S_C is estimated by the empirical estimator of CDF
    delta_c = rep(1, length(vect_y))
    data_c = data.frame(y = vect_c,
                        delta_c = delta_c)
    km_c = survival::survfit(survival::Surv(time = y , event = delta_c) ~ 1,
                             data = data_c,
                             type = "kaplan-meier")
  }
  c_time = c(0,km_c$time)
  c_surv = c(1,km_c$surv)
  N = length(vect_y)
  w = stats::approx(x = c_time, # weights are estimated by interpolation of S_C estimator
                    y = c_surv,
                    xout = vect_y,
                    method = "linear",
                    rule = 2)
  weights = ifelse(vect_delta == 0, 0, 1/N * 1/w$y)
  return(list(weights = weights,
              censoring_survfit = km_c))
}


test_conditional_copula2 = function(n, fct_tau, family, rate_cens,
                                    noyau, h, grid, marginal){
  
  ##
  # n : nombres d'oservations dans l'?chantillon simul?
  # dim : dimension du vecteur X des covariables (on peut faire plus que 1 dimension)
  # fonction_theta : la fonction qui donne le param?tre theta de la copule en fonction de X
  # family : famille de copules utilis?e
  # 1 : Gaussian copula
  # 3 : Clayton copula
  # 4 : Gumbel copula
  # 5 : Frank copula
  # rate_cens : proportion d'observation censur?es dans l'?chantillon simul?
  # noyau : noyau K utilis?
  # h : parametre h qui intervient dans la fonction noyau K
  # grid_pas : pas de la grille (de valeurs de X) pour laquelle on effectue les calculs
  # en chacun des points
  # borne_inf : Borne inf?rieure des valeurs de x de la grille
  # borne_sup : Borne sup?rieure des valeurs de x de la grille
  
  marginal = match.arg(as.character(marginal), c("unif", "lnorm"))
  
  # simulate T and U
  X = runif(n = n, min = 0, max = 1)
  tau = fct_tau(X)
  par = BiCopTau2Par(family = family, tau = tau)
  data = t(sapply(X = par, FUN =  BiCopSim, N = 1, family = family))
  if (marginal == "lnorm") data = apply(X = data, MARGIN = 2, FUN = qlnorm, meanlog = 0, sdlog = 1)
  data = as.data.frame(data)
  colnames(data) = c("T", "U")
  
  # simulate C
  if(rate_cens > 0){
    rate_cens_simul = -1
    while (abs(rate_cens_simul - rate_cens) > 0.05){
      res_optimize = optimize(f =  obj_fct_exp, lower = 0, upper = 10, 
                              x = data$T + data$U , target = rate_cens, tol = 0.02)
      data$C = exp(rexp(n = n, rate = res_optimize$minimum))
      
      # simulate observable var.
      data$eta = 1.0 * (data$T <= data$C)
      data$gamma = 1.0 * ((data$U + data$T) <= data$C)
      data$Y = pmin(data$T, data$C)
      data$Z = pmin(data$U, data$C - data$Y)
      # ifelse((data$U + data$T) <= data$C, pmin(data$U, data$C - data$T), 0)
      
      # verif. of censoring rates :
      rate_cens_simul = mean(data$U + data$T > data$C)
      print(c(paste0("rate_cens : ", rate_cens) ,
              paste0("rate_cens_simul : ", rate_cens_simul)))
    }
  } else {
    data$eta = 1.0
    data$gamma = 1.0
    data$Y = data$T
    data$Z = data$U
  }
  
  # pseudo-observations
  KM_T = survfit(formula = Surv(time = Y, event = eta) ~ 1,
                 data = data)
  data$A = 1 - approx(x = c(0, KM_T$time),
                      y = c(1, KM_T$surv),
                      xout = data$Y,
                      method = "constant",
                      rule = 2)$y
  
  KM_U = survfit(formula = Surv(time = Z, event = gamma) ~ 1,
                 data = data)
  data$B = 1 - approx(x = c(0, KM_U$time),
                      y = c(1, KM_U$surv),
                      xout = data$Z,
                      method = "constant",
                      rule = 2)$y
  
  
  # IPCW weights
  if(rate_cens > 0){
    poids = make_KM_weights(vect_y = data$Y + data$Z,
                            vect_delta = data$gamma)$weights
  } else {
    poids = rep(1, n)
  }
  
  # estimation
  
  est_par = c()
  est_tau = c()
  
  cat("Nb. of steps : ", length(grid), "\n")
  for (i in 1:length(grid)){
    cat("  ", i)
    #poids_nw = sapply(X = X, FUN = function(x){return(1 / h * noyau((x - grid[i]) / h))}) # pas besoin de sapply
    poids_nw = noyau((X - grid[i]) / h) #1 / h * 
    
    # if (family == 1 | family == 3){
    #   res = try(expr = BiCopEst(u1 = data$A,
    #                             u2 = data$B,
    #                             family = family,
    #                             method = "mle",
    #                             weights = poids * poids_nw / sum(poids * poids_nw)),
    #             silent = T)
    #   
    #   if (length(res) == 1){
    #     est_par[i] = NA
    #     est_tau[i] = NA
    #   } else {
    #     est_par[i] = res$par
    #     est_tau[i] = BiCopPar2Tau(family = family, par = res$par)
    #   }
    # }
    
    #if (family == 1) interval = c(0, 1)
    #if (family == 3) interval = c(0, 8)
    #if (family == 4) interval = c(1, 5)
    #if (family == 5) interval = c(0, 19)
    
    if (family == 1) {
      res = try(expr = optimize(f = loglike_copula,
                                interval = c(0, 1),
                                u1 = data$A,
                                u2 = data$B,
                                family = family,
                                poids = poids * poids_nw / sum(poids * poids_nw)),
                silent = T)
      
      if (length(res) == 1){
        est_par[i] = NA
        est_tau[i] = NA
      } else {
        est_par[i] = res$minimum
        est_tau[i] = BiCopPar2Tau(family = family, par = res$minimum)
      }
    }
    
    
    if (family == 3) {
      res = try(expr = optimize(f = loglike_copula,
                                interval = c(0, 8),
                                u1 = data$A,
                                u2 = data$B,
                                family = family,
                                poids = poids * poids_nw / sum(poids * poids_nw)),
                silent = T)
      
      if (length(res) == 1){
        est_par[i] = NA
        est_tau[i] = NA
      } else {
        est_par[i] = res$minimum
        est_tau[i] = BiCopPar2Tau(family = family, par = res$minimum)
      }
    }
    

    if (family == 4){
      res = try(expr = optimize(f = loglike_copula,
                                interval = c(1,5),
                                u1 = data$A,
                                u2 = data$B,
                                family = family,
                                poids = poids * poids_nw / sum(poids * poids_nw)),
                silent = T)
      
      if (length(res) == 1){
        est_par[i] = NA
        est_tau[i] = NA
      } else {
        est_par[i] = res$minimum
        est_tau[i] = BiCopPar2Tau(family = family, par = res$minimum)
      }
    }
    
    if (family == 5){
      res = try(expr = optimize(f = loglike_copula,
                                interval = c(0,19),
                                u1 = data$A,
                                u2 = data$B,
                                family = family,
                                poids = poids * poids_nw / sum(poids * poids_nw)),
                silent = T)
      
      if (length(res) == 1){
        est_par[i] = NA
        est_tau[i] = NA
      } else {
        est_par[i] = res$minimum
        est_tau[i] = BiCopPar2Tau(family = family, par = res$minimum)
      }
    }
  }
  
  exact_tau = fct_tau(grid)
  exact_par = BiCopTau2Par(family = family, tau = exact_tau)
  
  # output
  result = data.frame(
    x = grid,
    est_par = est_par,
    exact_par = exact_par,
    error_par = abs(est_par - exact_par),
    est_tau = est_tau,
    exact_tau = exact_tau,
    error_tau = abs(est_tau - exact_tau)
  )
  return(result)
}


make_result = function(v_family,
                        v_rate_cens,
                        v_n,
                        v_h,
                        n_repet,
                        grid,
                        fct_tau,
                        noyau,
                        marginal){
  
  param = expand.grid(iter = 1:n_repet,
                      h = v_h,
                      n = v_n,
                      rate_cens = v_rate_cens,
                      family = v_family)
  
  n_cores = detectCores() - 1
  cl = makeCluster(n_cores)
  registerDoParallel(cl)
  
  result = foreach(i = 1:dim(param)[1],
                   .combine = rbind,
                   .packages = c("copula",
                                 "VineCopula",
                                 "survival"),
                   .export= c("test_conditional_copula2",
                              "sigmoid",
                              "calc_rate_cens",
                              "obj_fct_exp",
                              "make_KM_weights",
                              "loglike_copula"),
                   .verbose = F) %dopar% {
                     res = test_conditional_copula2(n = param$n[i],
                                                    fct_tau = fct_tau,
                                                    family = param$family[i],
                                                    rate_cens = param$rate_cens[i],
                                                    noyau = noyau,
                                                    h = param$h[i],
                                                    grid = grid,
                                                    marginal)
                     return(cbind(family = param$family[i],
                                  rate_cens = param$rate_cens[i],
                                  n = param$n[i],
                                  h = param$h[i],
                                  iter = param$iter[i],
                                  res))
                   }
  stopCluster(cl)
  return(result)
}


