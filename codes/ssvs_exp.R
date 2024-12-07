library(MASS)
library(gridExtra)
library(runjags)
library(ggmcmc)
library(coda)
library(knitr)
library(ggplot2)

# change the beta values because poission might have extreme values
# change true beta values from 2 to 1 significantly changed the behavior of the results
# when beta=2, the three chains tends to disagree when burn-in is low, didn't test high burn in cases
# when beta=1, the three chains agrees more, the samples are more stable, the estimates of the fixed effects are more concentrated.

# negative binomial, mean, over-dispersion
# gamma regression, mean, 
# Beta Regression 
# baye tuning parameter selection hyperparameter
# 

sample_generate = function()
{
  # 100 subjects
  n = 100
  n_i = 10
  sigma = 0.5
  beta = c(1,1,1,0,0,0,0,0)
  cov_mat = matrix(c(0,0,0,0,0,0.9,0.48,0.06,0,0.48,0.4,0.1,0,0.06,0.1,0.1), ncol=4)
  #cov_mat = matrix(c(1,0,0,0,0,0.9,0.48,0.06,0,0.48,0.4,0.1,0,0.06,0.1,0.1 ),ncol=4)
  Q_mu = c(0,0,0,0)
  samples = matrix(nrow = 1000, ncol = 10)
  Qi = matrix(NA, nrow = n, ncol = 4)
  index = 1
  
  ## i is group id
  for (i in 1:n) {
    ## Qi are random effects for group i
    Qi[i,] = mvrnorm(1, Q_mu, cov_mat) 
    
    for (j in 1:n_i) {
      Xij = c(1,runif(7,-10,10))
      Zij = Xij[1:4]
      ratee = exp(Xij%*%beta+Zij %*% Qi[i,])
      yij = rexp(1, ratee)
      samples[index,1] = i
      samples[index,2:9] = Xij
      samples[index, 10] = yij
      index = index + 1
    } 
  }
  colnames(samples) = c("id","X1","X2","X3","X4","X5","X6","X7","X8","Y")
  return(list(samples = samples, re = Qi))
}




samples = sample_generate()
Y = samples$samples[,10]
X= samples$samples[,2:9]
## Confirm different relationships
samples.df = data.frame(samples$samples)


# double pareto model for the fixed effects

double_pareto = "model {
  # priors for r
  # dim = q(q-1)/2
  mu = rep(0, 28)
  rprior ~ dmnorm(mu, cov_mat)

  # priors for random effects
  for (j in 1:p){
    for(i in 1:n){ # n is number of groups
        Xi[i, j] ~ dnorm(0, 1)
    }
  }
  # priors for cov matrices
  # truncated normal
  h_sq = pow(h,2)
  for (k in 1:p){
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
  
  # constrains the R's to be 0 if lambda is 0.

  for (i in 1:p) {
     Lambda_mat[i,i] <- lambda[i]
     for (j in (i+1):p) {
        Lambda_mat[i,j] <- 0
        Lambda_mat[j,i] <- Lambda_mat[i,j]
     }
  }

  for (i in 1:p) {
     offset[i] <- (i - 1) * (2*p - i)/2
     Gamma_mat[i,i] <- 1
     for (j in (i+1):p) {
        Gamma_mat[j,i] <- ifelse(indicator3[i] == indicator3[j], rprior[offset[i] + j - i], 0) 
        Gamma_mat[i,j] <- 0
     }
  }
  
  # priors for fixed effects
  
  g ~ dgamma(0.5, nobs/2)
  non_scale_sigma2 ~ dgamma(0.5, 0.5)
  sigma2 = pow(non_scale_sigma2, -1)
  
  for (i in 1:l){
    #beta_phi_k[i] ~ dgamma(v, nu)
    #beta_phi_k_2[i] = pow(beta_phi_k[i], 2)
    beta_K_k[i] ~ dexp( 1 )
    
    beta_side1[i] ~ dnorm(0, sigma2/(g*beta_K_k[i]) ) 
    
    beta_side2[i] = 0
    latent_J[i] ~ dbern(0.5)
    latent_indicator1[i] <- ifelse(latent_J[i] > 0, 1, 0)
    latent_indicator2[i] <- ifelse(latent_J[i] <= 0, 1, 0)
    beta[i] <- (beta_side1[i] * latent_indicator1[i]) + (beta_side2[i] * latent_indicator2[i])
  }

  # likelihood
  for (i in 1:nobs) {
     #ratee[i] = exp(X[i, ] %*% beta[]+ X[i,] %*% Lambda_mat %*% Gamma_mat %*%  Xi[index[i], ]  )
     #ratee[i] = exp( X[i,] %*% beta[])
     ratee[i] = exp( X[i,] %*% Lambda_mat %*% Gamma_mat %*%  Xi[index[i], ]  )
     #ratee[i] = ifelse( X[i,] %*% Lambda_mat %*% Gamma_mat %*%  Xi[index[i], ] >= 10, 10, X[i,] %*% Lambda_mat %*% Gamma_mat %*%  Xi[index[i], ] )
     true_rate[i] = ifelse( ratee[i] > 10, 10, ratee[i] )
     Y[i] ~ dexp( true_rate[i] )
    #Y[i] ~ dpois( exp(X[i, ] %*% beta[]+ X[i,] %*% Lambda_mat %*% Gamma_mat %*%  Xi[index[i], ]) )
  }
}
"

# prior for lower-triangular random effect matrix
cov_mat_R = diag(1,nrow = 28)

# adapt is chooing the variance of the proposal normal

cat(double_pareto, file="double_pareto_ssvs.txt")

double_pareto_dat <- list(p=8, l=8, nobs=1000, n = 100, cov_mat= cov_mat_R, X=X, Y=Y, h= 0.1, v=0.01, nu=0.01,
                          index = samples.df$id)

vars <- c("beta", "Gamma_mat", "Lambda_mat", "Xi", "g", "sigma2", "beta_K_k")


double_pareto_out <- run.jags("double_pareto_ssvs.txt", vars, data=double_pareto_dat, n.chains=3,
                              adapt=100, burnin=100, sample=100)


summary(double_pareto_out)


