library(MASS)
library(gridExtra)
library(runjags)
library(ggmcmc)
library(coda)
library(knitr)
library(ggplot2)
library(readxl)



mcmc.iter = 1000


# this is the final model we are shooting for
# glmmTMB(freq~ MAT+ TSD+ PCV+X.Sand+X.N+BA+(1+VH+LS |site),family =nbinom1(link="log"), offset = log(quads), data = data)

# brm.full.fit = brm(
#   freq ~ MAT + MAP + TSD + PCV + X.Sand + X.N + BA +
#     I(MAP^2) + I(MAT^2) + I(TSD^2) + I(X.Sand^2) + I(X.N^2) +
#     I(BA^2) * (VH + LS + LMA + LCC) +
#     I(VH^2) + I(LS^2) + I(LMA^2) + I(LCC^2) +
#     (1 + VH + LS + LMA + LCC | site) +
#     (1 + MAP + MAT + TSD + X.Sand + X.N + BA | sp) +
#     offset(log(quads)),
#   family = "negbinomial",
#   cores = 3,
#   chains = 3,
#   control = list(adapt_delta = 0.95, max_treedepth = 12),
#   data = df
# )

data = data.frame(read_xlsx("./Data/forest_data.xlsx"))
data$LMA = 1 / data$SLA ## leaf mass area is inverse of specific leaf area
df = data

Y = df$freq

## Extra model, added PCV^2
X = model.matrix(
  ~ -1 + MAT + MAP + TSD + PCV + X.Sand + X.N + BA +
    I(MAT ^ 2) + I(MAP ^ 2) + I(TSD ^ 2) + I(PCV ^ 2) + I(X.Sand ^ 2) + I(X.N ^ 2) +
    I(BA ^ 2) * (VH + LS + LMA + LCC) +
    I(VH ^ 2) + I(LS ^ 2) + I(LMA ^ 2) + I(LCC ^ 2),
  data = df
)

saved_names = c(colnames(X))
X = as.matrix(scale(X))
colnames(X) = saved_names
head(X)

## For site
Z1 = cbind(1, X[, c("VH", "LS", "LMA", "LCC")])
colnames(Z1) = c("intercept", colnames(X[, c("VH", "LS", "LMA", "LCC")]))

## For sp
Z2 = cbind(1, X[, c("MAT", "MAP", "TSD", "PCV", "X.Sand", "X.N", "BA")]) ## Added PCV
colnames(Z2) = c("intercept", colnames(X[, c("MAT", "MAP", "TSD", "PCV", "X.Sand", "X.N", "BA")]))

offset = log(df$quads)

q1 = ncol(Z1)
q2 = ncol(Z2)
l = ncol(X)
nobs = nrow(df)
index1 = as.numeric(as.factor(data$site))
index2 = as.numeric(as.factor(data$sp))
n1 = length(unique(index1))
n2 = length(unique(index2))



double_pareto = "model {


  h_sq = pow(h,2)

# Random Effects ----------------------------------------------------------
## Z1

  for (j in 1:q1){
    for(i in 1:n1){ # n is number of groups
        m_k_Z1[i, j] ~ dgamma(1,1)
        m_k2_Z1[i, j] = pow(m_k_Z1[i, j], 2) / 2
        tau_k_Z1[i, j] ~ dexp(m_k2_Z1[i, j])
        Xi_Z1[i,j] ~ dnorm(0, 1/ tau_k_Z1[i, j]) # need to use precision here
    }
  }
  # priors for cov matrices
  # truncated normal
  for (k in 1:q1){
    gamma_phik2_Z1[k] ~ dgamma(0.5, 0.5)
    phik2_Z1[k] = pow(gamma_phik2_Z1[k], 2) # no need to inverse here b/c dnorm uses precision
    side1_unscale_Z1[k] ~  dnorm(0, phik2_Z1[k]*h_sq) T(0,)
    side1_Z1[k] = side1_unscale_Z1[k]*0.5
    side2_Z1[k] = 0
    pk0_Z1[k] ~ dbern(0.5)
    indicator1_Z1[k] <- ifelse(pk0_Z1[k] > 0, 1, 0)
    indicator2_Z1[k] <- ifelse(pk0_Z1[k] <= 0, 1, 0)
    lambda_Z1[k] <- (side1_Z1[k] * indicator1_Z1[k]) + (side2_Z1[k] * indicator2_Z1[k])
    indicator3_Z1[k] = ifelse(lambda_Z1[k] == 0, 0, 1)
  }


# Random Effects ----------------------------------------------------------
## Z2

  for (j in 1:q2){
    for(i in 1:n2){ # n is number of groups
        m_k_Z2[i, j] ~ dgamma(1,1)
        m_k2_Z2[i, j] = pow(m_k_Z2[i, j], 2) / 2
        tau_k_Z2[i, j] ~ dexp(m_k2_Z2[i, j])
        Xi_Z2[i,j] ~ dnorm(0, 1/ tau_k_Z2[i, j]) # need to use precision here
    }
  }
  # priors for cov matrices
  # truncated normal
  for (k in 1:q2){
    gamma_phik2_Z2[k] ~ dgamma(0.5, 0.5)
    phik2_Z2[k] = pow(gamma_phik2_Z2[k], 2) # no need to inverse here b/c dnorm uses precision
    side1_unscale_Z2[k] ~  dnorm(0, phik2_Z2[k]*h_sq) T(0,)
    side1_Z2[k] = side1_unscale_Z2[k]*0.5
    side2_Z2[k] = 0
    pk0_Z2[k] ~ dbern(0.5)
    indicator1_Z2[k] <- ifelse(pk0_Z2[k] > 0, 1, 0)
    indicator2_Z2[k] <- ifelse(pk0_Z2[k] <= 0, 1, 0)
    lambda_Z2[k] <- (side1_Z2[k] * indicator1_Z2[k]) + (side2_Z2[k] * indicator2_Z2[k])
    indicator3_Z2[k] = ifelse(lambda_Z2[k] == 0, 0, 1)
  }


# Fixed Effects -----------------------------------------------------------

  g ~ dgamma(0.5, nobs/2)
  non_scale_sigma2 ~ dgamma(0.5, 0.5)
  sigma2 = pow(non_scale_sigma2, -1)
  beta_int ~ dnorm(0, 0.01)

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



# Other param -------------------------------------------------------------

  r ~ dunif(0, 100)


# Likelihood --------------------------------------------------------------

  for (i in 1:nobs) {
    for(j in 1:q1){
      rho_Z1[i,j] = lambda_Z1[j] * Xi_Z1[index1[i], j]
    }
    for(j in 1:q2){
      rho_Z2[i,j] = lambda_Z2[j] * Xi_Z2[index2[i], j]
    }

    log(mu[i]) = exp(beta_int + X[i, ] %*% beta[] + Z1[i, ] %*% rho_Z1[i,] + Z2[i, ] %*% rho_Z2[i,] + offset[i])
    p[i] = r / (r + mu[i])

    # Y[i] ~ dpois( exp(X[i, ] %*% beta[]) )
    # Y[i] ~ dpois( exp(X[i, ] %*% beta[] + Z[i,] %*% Lambda_mat %*% Xi[index[i], ]) )
    # Y[i] ~ dpois( exp(X[i, ] %*% beta[] + Z[i, ] %*% Xi[index[i], ] ))
    # Y[i] ~ dpois(exp(beta_int + X[i, ] %*% beta[] + Z1[i, ] %*% rho_Z1[i,] + Z2[i, ] %*% rho_Z2[i,] ))
    Y[i] ~ dnegbin(p[i], r)
  }
}
"


cat(double_pareto, file = "double_pareto_ssvs.txt")

double_pareto_dat <- list(
  q1 = q1,
  q2 = q2,
  l = l,
  nobs = nobs,
  n1 = n1,
  n2 = n2,
  X = X,
  Z1 = Z1,
  Z2 = Z2,
  Y = Y,
  h = 1,
  v = 0.01,
  nu = 0.01,
  index1 = index1,
  index2 = index2,
  offset = offset
)
vars <- c(
  "beta_int",
  "beta",
  "g",
  "sigma2",
  "latent_indicator1",
  "r",
  "indicator3_Z1",
  "lambda_Z1",
  "indicator3_Z2",
  "lambda_Z2"
)


double_pareto_out = run.jags(
  "double_pareto_ssvs.txt",
  vars,
  data = double_pareto_dat,
  n.chains = 3,
  adapt = mcmc.iter,
  burnin = mcmc.iter,
  sample = mcmc.iter * 3,
  method = "parallel"
)


save.image("./Results/real_fit_diag.RData")
saveRDS(double_pareto_out, "./Results/real_fit_diag.rds")
