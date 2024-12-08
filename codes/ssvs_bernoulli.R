library(MASS)
library(gridExtra)
library(runjags)
library(ggmcmc)
library(coda)
library(knitr)
library(ggplot2)
library(reshape2)


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
      Xij = c(1,runif(7,-1,1))
      Zij = Xij[1:4]
      mean_link = exp(Xij%*%beta+Zij %*% Qi[i,]) / (1+exp(Xij%*%beta+Zij %*% Qi[i,]))
      yij = rbinom(1,1, mean_link)
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
  mu = rep(0, 28)
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
    link_mean[i] = exp(X[i, ] %*% beta[]+ X[i,] %*% Lambda_mat %*% Gamma_mat %*%  Xi[index[i], ])/(1+exp(X[i, ] %*% beta[]+ X[i,] %*% Lambda_mat %*% Gamma_mat %*%  Xi[index[i], ]))
    Y[i] ~ dbern( link_mean[i] )
  }
}
"

horse_shoe = "model {
  # priors for r
  # dim = q(q-1)/2
  mu = rep(0, 28)
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
    link_mean[i] = exp(X[i, ] %*% beta[]+ X[i,] %*% Lambda_mat %*% Gamma_mat %*%  Xi[index[i], ])/(1+exp(X[i, ] %*% beta[]+ X[i,] %*% Lambda_mat %*% Gamma_mat %*%  Xi[index[i], ]))
    Y[i] ~ dbern( link_mean[i] )
  }
}
"
cov_mat_R = diag(1,nrow = 28)

cat(horse_shoe, file="horse_shoe_ssvs.txt")
cat(double_pareto, file="double_pareto_ssvs.txt")

double_pareto_dat <- list(p=8, l=8, nobs=1000, n = 100, cov_mat= cov_mat_R, X=X, Y=Y, h= 0.1, v=0.01, nu=0.01,
                          index = samples.df$id)
horse_shoe_dat <- list(p=8, l=8, nobs=1000, n = 100, cov_mat= cov_mat_R, X=X, Y=Y, h= 0.1, index=samples.df$id)
vars <- c("beta", "Gamma_mat", "Lambda_mat", "Xi", "g", "sigma2", "beta_K_k",
          "latent_J", "indicator1")
horse_shoe_out <- run.jags("horse_shoe_ssvs.txt", vars, data=horse_shoe_dat, n.chains=3,
                           adapt=100, burnin=300, sample=100)
double_pareto_out <- run.jags("double_pareto_ssvs.txt", vars, data=double_pareto_dat, n.chains=3,
                              adapt=100, burnin=300, sample=100)


summary(horse_shoe_out) 
summary(double_pareto_out)

horse_shoe_par.est = data.frame(as.mcmc(horse_shoe_out))
double_pareto_par.est = data.frame(as.mcmc(double_pareto_out))







plot(horse_shoe_par.est$beta.1.)
plot(horse_shoe_par.est$beta.2.)
plot(horse_shoe_par.est$beta.3.)
plot(horse_shoe_par.est$beta.4.)
plot(horse_shoe_par.est$beta.5.)
plot(horse_shoe_par.est$beta.6.)
plot(horse_shoe_par.est$beta.7.)
plot(horse_shoe_par.est$beta.8.)

# the random effects are not being shrinking to 0,
# this means we are over-selecting?

horse_shoe_Lambda_mat = matrix(colMeans( horse_shoe_par.est[,73:136] ), ncol = 8, nrow =8 )
horse_shoe_Lambda_mat

double_pareto_Lambda_mat = matrix(colMeans( double_pareto_par.est[,73:136] ), ncol = 8, nrow =8 )
double_pareto_Lambda_mat









# Fixed Effects -----------------------------------------------------------

hs.beta.est = data.frame(as.mcmc(horse_shoe_out, vars = "beta[1:8]"))
pareto.beta.est = data.frame(as.mcmc(double_pareto_out, vars = "beta[1:8]"))
truth.beta = data.frame(matrix(c(1,1,1,0,0,0,0,0), nrow = 1))
names(truth.beta) = names(hs.beta.est)

hs.beta.melt = cbind(melt(hs.beta.est), model = "Horsehoe")
pareto.beta.melt = cbind(melt(pareto.beta.est), model = "Pareto")
truth.beta.melt = cbind(melt(truth.beta), model = "Truth")

beta.df = rbind(hs.beta.melt, pareto.beta.melt)


ggplot(beta.df, aes(x = value, col = model)) + geom_density(linewidth = 1.3) +
  facet_wrap(~variable, scales = "free") +
  geom_vline(data = truth.beta.melt, aes(xintercept = value))




# Random Effects ----------------------------------------------------------







# Model Selection ---------------------------------------------------------


hs.var.est = data.frame(as.mcmc(horse_shoe_out, vars = c("latent_J", "indicator1")))
hs.var.est = cbind(hs.var.est[,c(1:8)], "__", hs.var.est[,c(9:16)])
hs.var.cvec = apply(hs.var.est, 1, paste0, collapse = "")

pareto.var.est = data.frame(as.mcmc(double_pareto_out, vars = c("latent_J", "indicator1")))
pareto.var.est = cbind(pareto.var.est[,c(1:8)], "__", pareto.var.est[,c(9:16)])
pareto.var.cvec = apply(pareto.var.est, 1, paste0, collapse = "")




sort(table(hs.var.cvec), decreasing = T)[1:5]
sort(table(pareto.var.cvec), decreasing = T)[1:5]
## Compare to truth
paste0(c(as.numeric(truth.beta), "__", c(0, 1, 1, 1, rep(0, 4))), collapse = "") ## Ugly code, I know
       