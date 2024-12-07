library(MASS)
library(gridExtra)
library(runjags)
library(ggmcmc)
library(coda)
library(knitr)
library(ggplot2)


args = commandArgs(trailingOnly = TRUE)
ind = as.numeric(args[1])

list2env(readRDS(paste0("./Data/Sims/Poisson_case_1_sim_", ind, ".rds")), globalenv())


double_pareto = "model {


# Random Effects ----------------------------------------------------------

  for (j in 1:q){
    for(i in 1:n){ # n is number of groups
        m_k[i, j] ~ dgamma(1,1)
        m_k2[i, j] = pow(m_k[i, j], 2) / 2
        tau_k[i, j] ~ dexp(m_k2[i, j])
        Xi[i,j] ~ dnorm(0, 1/ tau_k[i, j]) # need to use precision here
    }
  }
  # priors for cov matrices
  # truncated normal
  h_sq = pow(h,2)
  for (k in 1:q){
    gamma_phik2[k] ~ dgamma(0.5, 0.5)
    phik2[k] = pow(gamma_phik2[k], 2) # no need to inverse here b/c dnorm uses precision
    side1_unscale[k] ~  dnorm(0, phik2[k]*h_sq) T(0,)
    side1[k] = side1_unscale[k]*0.5
    side2[k] = 0
    pk0[k] ~ dbern(0.5)
    indicator1[k] <- ifelse(pk0[k] > 0, 1, 0)
    indicator2[k] <- ifelse(pk0[k] <= 0, 1, 0)
    lambda[k] <- (side1[k] * indicator1[k]) + (side2[k] * indicator2[k])
    indicator3[k] = ifelse(lambda[k] == 0, 0, 1)
  }

# Fixed Effects -----------------------------------------------------------

  g ~ dgamma(0.5, nobs/2)
  non_scale_sigma2 ~ dgamma(0.5, 0.5)
  sigma2 = pow(non_scale_sigma2, -1)

  for (i in 1:l){
    beta_phi_k[i] ~ dgamma(v, nu)
    beta_phi_k_2[i] = pow(beta_phi_k[i], 2)
    beta_K_k[i] ~ dexp( beta_phi_k_2[i] )

    beta_side1[i] ~ dnorm(0, sigma2/(g*beta_K_k[i]) )

    beta_side2[i] = 0
    latent_J[i] ~ dbern(0.5)
    latent_indicator1[i] <- ifelse(latent_J[i] > 0, 1, 0)
    latent_indicator2[i] <- ifelse(latent_J[i] <= 0, 1, 0)
    beta[i] <- (beta_side1[i] * latent_indicator1[i]) + (beta_side2[i] * latent_indicator2[i])
  }


# Likelihood --------------------------------------------------------------

  for (i in 1:nobs) {
    for(j in 1:q){
      rho[i,j] = lambda[j] * Xi[index[i], j]
    }
  
    # Y[i] ~ dpois( exp(X[i, ] %*% beta[]) )
    # Y[i] ~ dpois( exp(X[i, ] %*% beta[] + Z[i,] %*% Lambda_mat %*% Xi[index[i], ]) )
    # Y[i] ~ dpois( exp(X[i, ] %*% beta[] + Z[i, ] %*% Xi[index[i], ] ))
    Y[i] ~ dpois( exp(X[i, ] %*% beta[] + Z[i, ] %*% rho[i,] ))
  }
}
"

cat(double_pareto, file = "double_pareto_ssvs.txt")




hyper.param = expand.grid(c(0.01, 1, 5), c(0.1, 1, 10))

for(hyper.ind in 1:nrow(hyper.param)){
  double_pareto_dat <- list(
    q = q,
    l = l,
    nobs = length(index),
    n = n,
    X = X,
    Z = Z,
    Y = Y,
    h = hyper.param[hyper.ind, 2],
    v = hyper.param[hyper.ind, 1],
    nu = hyper.param[hyper.ind, 1],
    index = index
  )
  vars <- c("beta", "lambda", "Xi", "g", "sigma2", "rho", "latent_indicator1")
  double_pareto_out = run.jags(
    "double_pareto_ssvs.txt",
    vars,
    data = double_pareto_dat,
    n.chains = 3,
    adapt = 1000,
    burnin = 3000,
    sample = 1000,
    method = "parallel"
  )

  saveRDS(double_pareto_out, paste0("./Results/Sims/Diag/Case1/Poisson_diag_hyper_",hyper.ind,"_ind_",ind,".rds"))
}


