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


sample_generate = function()
{
  # 100 subjects
  n = 100
  n_i = 10
  sigma = 0.5
  beta = c(1,1,1,0,0,0,0,0)
  cov_mat = matrix(c(0,0,0,0,0,0.9,0.48,0.06,0,0.48,0.4,0.1,0,0.06,0.1,0.1 ),ncol=4)
  #cov_mat = matrix(c(1,0,0,0,0,0.9,0.48,0.06,0,0.48,0.4,0.1,0,0.06,0.1,0.1 ),ncol=4)
  Q_mu = c(0,0,0,0)
  samples = matrix(nrow = 1000, ncol = 10)
  # yij = Xij^T beta + Zij^T Q + epslion
  
  Qi = matrix(NA, nrow = n, ncol = 4)
  index = 1
  for (i in 1:n) { ## i is group id
    Qi[i,] = mvrnorm(1, Q_mu, cov_mat) ## Qi are random effects for group i
    print(Qi[i,])
    for (j in 1:n_i) {
      Xij = c(1,runif(7,-1,1))
      Zij = Xij[1:4]
      #epsij = rnorm(n=1,mean=0,sd=sigma)
      yij = rpois(1,exp(Xij%*%beta+Zij %*% Qi[i,]))
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

double_pareto = "model {
  # priors for r
  # dim = q(q-1)/2
  mu = rep(0, 6)
  rprior ~ dmnorm(mu, cov_mat)

  # priors for random effects
  for (j in 1:p){
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

  # likelihood
  for (i in 1:nobs) {
    #Y[i] ~ dnorm(X[i, ] %*% beta[]+ X[i,1:4] %*% Lambda_mat %*% Gamma_mat %*%  Xi[index[i], ], 2) # We can estimate the variance
    Y[i] ~ dpois( exp(X[i, ] %*% beta[]+ X[i,1:4] %*% Lambda_mat %*% Gamma_mat %*%  Xi[index[i], ]) )
  }
}
"

horse_shoe = "model {
  # priors for r
  # dim = q(q-1)/2
  mu = rep(0, 6)
  rprior ~ dmnorm(mu, cov_mat)

  # priors for random effects
  for (j in 1:p){
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
    nu2[i] ~ dt(0,1,1) T(0,) # cauchy is special T
    beta_K_k[i] ~ dt(0,nu2[i],1) T(0,) # cauchy is special T
    

    beta_side1[i] ~ dnorm(0, sigma2/(g*beta_K_k[i]) ) 
    beta_side2[i] = 0
    latent_J[i] ~ dbern(0.5)
    latent_indicator1[i] <- ifelse(latent_J[i] > 0, 1, 0)
    latent_indicator2[i] <- ifelse(latent_J[i] <= 0, 1, 0)
    beta[i] <- (beta_side1[i] * latent_indicator1[i]) + (beta_side2[i] * latent_indicator2[i])
  }

  # likelihood
  for (i in 1:nobs) {
    #Y[i] ~ dnorm(X[i, ] %*% beta[]+ X[i,1:4] %*% Lambda_mat %*% Gamma_mat %*%  Xi[index[i], ], 2) # We can estimate the variance
    Y[i] ~ dpois( exp(X[i, ] %*% beta[]+ X[i,1:4] %*% Lambda_mat %*% Gamma_mat %*%  Xi[index[i], ]) )
  }
}
"

cat(horse_shoe, file="ssvs.txt")
cov_mat_R = diag(1,nrow = 6)

#dat <- list(p=4, l=8, nobs=1000, n = 100, cov_mat= cov_mat_R, X=X, Y=Y, h= 0.1, v=0.01, nu=0.01,
#            index = samples.df$id)
dat <- list(p=4, l=8, nobs=1000, n = 100, cov_mat= cov_mat_R, X=X, Y=Y, h= 0.1, index=samples.df$id)
vars <- c("beta", "Gamma_mat", "Lambda_mat", "Xi", "g", "sigma2", "beta_K_k")
out <- run.jags("ssvs.txt", vars, data=dat, n.chains=3,
                adapt=100, burnin=300, sample=100)
summary(out) 

par.est = data.frame(as.mcmc(out))
plot(par.est$beta.1.)
plot(par.est$beta.2.)
plot(par.est$beta.3.)
plot(par.est$beta.4.)
plot(par.est$beta.5.)





names(par.est)[9]
names(par.est)[24]
names(par.est)[25]
names(par.est)[40]

## see if we can reconstruct the true covariance matrix
Xi = matrix(colMeans( par.est[,41:440] ), ncol = 4, nrow =100 )
Gamma_mat = matrix(colMeans( par.est[,9:24] ), ncol = 4, nrow =4 )
Lambda_mat = matrix(colMeans( par.est[,25:40] ), ncol = 4, nrow =4 )
par.est
Gamma_mat
Lambda_mat

round(Lambda_mat %*% Gamma_mat %*% t(Gamma_mat) %*% (Lambda_mat),4)

true_cov = matrix(c(0,0,0,0,0,0.9,0.48,0.06,0,0.48,0.4,0.1,0,0.06,0.1,0.1 ),ncol=4)
chol(true_cov, pivot = T) # # chol decomp
decomp_mat = matrix(chol(true_cov,pivot = T)[1:16],nrow = 4)
t(decomp_mat) %*% decomp_mat
t(decomp_mat)
Lambda_mat %*% Gamma_mat




plot(par.est$g)
plot(par.est$sigma2)




plot(par.est$Lambda_mat.1.1.)
plot(par.est$Lambda_mat.2.2.)
plot(par.est$Lambda_mat.3.3.)
plot(par.est$Lambda_mat.4.4.)


for(i in 1:nrow(par.est)){
  lambda.mat.tmp = diag(c(par.est$Lambda_mat.1.1.[i],
                          par.est$Lambda_mat.2.2.[i],
                          par.est$Lambda_mat.3.3.[i],
                          par.est$Lambda_mat.4.4.[i]))
  gamma.mat.tmp = matrix(par.est[i,9:24], ncol = 4)
}

hist(rgamma(1000,0.1,0.1))





