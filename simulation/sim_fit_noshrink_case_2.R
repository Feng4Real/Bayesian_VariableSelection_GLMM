library(MASS)
library(gridExtra)
library(runjags)
library(ggmcmc)
library(coda)
library(knitr)
library(ggplot2)


args = commandArgs(trailingOnly = TRUE)
ind = as.numeric(args[1])


list2env(readRDS(paste0("./Data/Sims/Poisson_case_2_sim_", ind, ".rds")), globalenv())


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
    lambda[k] ~  dnorm(0, phik2[k]*h_sq) T(0,)
  }

# Fixed Effects -----------------------------------------------------------

  non_scale_sigma2 ~ dgamma(0.5, 0.5)
  sigma2 = pow(non_scale_sigma2, -1)

  for (i in 1:l){
    beta[i] ~ dnorm(0, sigma2)
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

hyper.noshrink.ind = c(1, 4, 7)

for(hyper.ind in hyper.noshrink.ind){
  double_pareto_dat <- list(
    q = q,
    l = l,
    nobs = length(index),
    n = n,
    X = X,
    Z = Z,
    Y = Y,
    h = hyper.param[hyper.ind, 2],
    index = index
  )
  vars <- c("beta", "lambda", "Xi", "sigma2", "rho")
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
  
  saveRDS(double_pareto_out, paste0("./Results/Sims/Noshrink/Case2/Poisson_noshrink_hyper_",hyper.ind,"_ind_",ind,".rds"))
}

