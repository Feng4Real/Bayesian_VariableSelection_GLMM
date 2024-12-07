library(MASS)
library(gridExtra)
library(runjags)
library(ggmcmc)
library(coda)
library(knitr)
library(ggplot2)
library(readxl)
library(dplyr)
library(purrr)
library(tidyr)

load("./Results/real_fit_diag.RData")


# Beta --------------------------------------------------------------------

beta.est = as.mcmc(double_pareto_out, vars = "beta")
beta.nz = data.frame(beta.est) %>% mutate(across(everything(), ~ if_else(. != 0, 1, 0)))
colnames(beta.est) = c("Intercept", colnames(X))
plot(beta.est[, 1])
plot(beta.est[, 2])
plot(beta.est[, 3])


## Find most frequent betas
beta.nz.mean = apply(beta.est, 2, function(x) {
  sum(x != 0)
})
names(beta.nz.mean) = c("Intercept", colnames(X))
sort(beta.nz.mean, decreasing = T)

beta.nz.collapse = apply(beta.nz, 1, paste0, collapse = "")
beta.nz.collapse.tab = table(beta.nz.collapse)
sort(beta.nz.collapse.tab, decreasing = T)[1:5]




# Lambda site ------------------------------------------------------------------

lambda.Z1.est = as.mcmc(double_pareto_out, vars = "lambda_Z1")
lambda.Z1.nz = data.frame(lambda.Z1.est) %>% mutate(across(everything(), ~ if_else(. != 0, 1, 0)))
colnames(lambda.Z1.est) = colnames(Z1)

plot(lambda.Z1.est[, 1])
plot(lambda.Z1.est[, 2])
plot(lambda.Z1.est[, 3])

colMeans(lambda.Z1.est)


## Find most frequent lambdas
lambda.Z1.nz.mean = apply(lambda.Z1.est, 2, function(x) {
  sum(x != 0)
})
names(lambda.Z1.nz.mean) = colnames(Z1)
sort(lambda.Z1.nz.mean, decreasing = T)

lambda.Z1.nz.collapse = apply(lambda.Z1.nz, 1, paste0, collapse = "")
lambda.Z1.nz.collapse.tab = table(lambda.Z1.nz.collapse)
sort(lambda.Z1.nz.collapse.tab, decreasing = T)[1:5]





# Lambda sp ------------------------------------------------------------------

lambda.Z2.est = as.mcmc(double_pareto_out, vars = "lambda_Z2")
lambda.Z2.nz = data.frame(lambda.Z2.est) %>% mutate(across(everything(), ~ if_else(. != 0, 1, 0)))
colnames(lambda.Z2.est) = colnames(Z2)

plot(lambda.Z2.est[, 1])
plot(lambda.Z2.est[, 2])
plot(lambda.Z2.est[, 3])

colMeans(lambda.Z2.est)


## Find most frequent lambdas
lambda.Z2.nz.mean = apply(lambda.Z2.est, 2, function(x) {
  sum(x != 0)
})
names(lambda.Z2.nz.mean) = colnames(Z2)
sort(lambda.Z2.nz.mean, decreasing = T)

lambda.Z2.nz.collapse = apply(lambda.Z2.nz, 1, paste0, collapse = "")
lambda.Z2.nz.collapse.tab = table(lambda.Z2.nz.collapse)
sort(lambda.Z2.nz.collapse.tab, decreasing = T)[1:5]


# Lambdas ------------------------------------------------

betalambda.nz = cbind(lambda.Z1.nz, lambda.Z2.nz)
colnames(betalambda.nz) = c(paste0("RE.Site.", colnames(Z1)),
                            paste0("RE.sp.", colnames(Z2)))


## Find most frequent lambda and beta
betalambda.nz.collapse = apply(betalambda.nz, 1, paste0, collapse = "")
betalambda.nz.collapse.tab = table(betalambda.nz.collapse)
sort(betalambda.nz.collapse.tab, decreasing = T)[1:5]


# Lambda and Beta together ------------------------------------------------

betalambda.nz = cbind(beta.nz, lambda.Z1.nz, lambda.Z2.nz)
colnames(betalambda.nz) = c(paste0("FE.", colnames(X)),
                            paste0("RE.Site.", colnames(Z1)),
                            paste0("RE.sp.", colnames(Z2)))


## Find most frequent lambda and beta
betalambda.nz.collapse = apply(betalambda.nz, 1, paste0, collapse = "")
betalambda.nz.collapse.tab = table(betalambda.nz.collapse)
sort(betalambda.nz.collapse.tab, decreasing = T)[1:5]



