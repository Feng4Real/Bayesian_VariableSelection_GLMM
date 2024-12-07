library(MASS)
library(gridExtra)
library(runjags)
library(ggmcmc)
library(coda)
library(knitr)
library(ggplot2)



for (sim_ind in 1:100) {
  n = 120 ## Number of groups
  n_i = 15 ## Number of samples per group
  
  
  # Sim fixed effects -------------------------------------------------------
  
  l = 10 ## Number of fixed effect parameters
  beta = numeric(l) ## Betas
  true_beta_length = 5 ## Number of non-zero fixed effects
  non_zero_beta = c(2, 5, 7, 8, 10)
  beta[non_zero_beta] = runif(true_beta_length, -0.4, 0.4)
  beta[1] = 2 ## Fix intercept to be non-zero
  
  
  
  # Sim random effects ------------------------------------------------------
  
  q = l ## Number of random effect parameters
  num_non_zero_Q = 3 ## Number of non-zero random effects
  non_zero_Q = c(1, 5, 10)

  Z_cov_mat = matrix(0, nrow = q, ncol = q)
  diag(Z_cov_mat)[non_zero_Q] = c(0.08, 0.15, 0.06)
  Z_cov_mat[1,5] = Z_cov_mat[5,1] = 0.04
  Z_cov_mat[1,10] = Z_cov_mat[10,1] = 0.02
  Z_cov_mat[5,10] = Z_cov_mat[10,5] = 0.09
  
  Q_mu = rep(0, q)
  Q = matrix(NA, nrow = n, ncol = q)
  ind = 1
  
  index = rep(NA, n * n_i)
  X = matrix(NA, nrow = 0, ncol = l)
  Y = matrix(NA, nrow = 0, ncol = 1)
  Z = matrix(NA, nrow = 0, ncol = q)
  Z.rho = X.beta = matrix(NA, nrow = 0, ncol = 1)
  
  ## i is group id
  for (i in 1:n) {
    ## Qi are random effects for group i
    Q[i, ] = mvrnorm(1, Q_mu, Z_cov_mat)
    
    for (j in 1:n_i) {
      Xij = c(1, rnorm(l - 1, 0, 1)) ## Simulate covs
      Zij = Xij[1:q] ## Grab covs for random effects too
      yij = rpois(1, exp(Xij %*% beta + Zij %*% Q[i, ])) ## Simulate response
      index[ind] = i
      X = rbind(X, Xij)
      Z = rbind(Z, Zij)
      Y = rbind(Y, yij)
      X.beta = rbind(X.beta, Xij %*% beta)
      Z.rho = rbind(Z.rho, Zij %*% Q[i, ])
      ind = ind + 1
    }
  }
  

  Y = as.numeric(Y)

  r_dim = q * (q - 1) / 2 ## Dimension of off-diagonal
  cov_mat_R = diag(1, nrow = r_dim)
  
  dat.list <- list(
    q = q,
    l = l,
    r_dim = r_dim,
    nobs = length(Y),
    n = n,
    cov_mat = cov_mat_R,
    X = X,
    Z = Z,
    Y = Y,
    h = 0.1,
    v = 0.01,
    nu = 0.01,
    index = index,
    beta = beta,
    Z_diag = Z_diag,
    Q = Q
  )
  
  
  saveRDS(dat.list, file = paste0("../../Data/Sims/Poisson_case_1_sim_", sim_ind, ".rds"))
  
  print(sim_ind)
}


