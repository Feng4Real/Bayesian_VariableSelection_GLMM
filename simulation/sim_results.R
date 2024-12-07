
library(runjags)
library(ggplot2)
library(reshape2)
library(dplyr)
library(coda)


runjags.options(silent.runjags = TRUE, silent.jags = TRUE)




# Noshrink ----------------------------------------------------------------


hyper.param = expand.grid(c(0.01, 1, 5), c(0.1, 1, 10))
hyper.noshrink.ind = c(1, 4, 7)

beta.est.noshrink = beta.inprob.est.noshrink =
  beta.best.noshrink = beta.est.resid.noshrink =
  array(NA, dim = c(100, 10, 2, 3))

for (case in 1:2) {
  for (hyper.ind in 1:3) {
    bad.ind.noshrink = c()
    
    
    # Full Model --------------------------------------------------------------
    
    beta.est = beta.inprob.est = beta.best = matrix(NA, nrow = 100, ncol = 10)
    for (ind in 1:100) {
      tmp.noshrink = try(readRDS(
        file = paste0(
          "../Results/Sims/Noshrink/Case",case,"/Poisson_noshrink_hyper_",
          hyper.noshrink.ind[hyper.ind],
          "_ind_",
          ind,
          ".rds"
        )
      ))
      
      list2env(readRDS(paste0("../Data/Sims/Poisson_case_",case,"_sim_", ind, ".rds")), globalenv())
      if (!("try-error" %in% class(tmp.noshrink))) {
        ## Look at beta estimates
        beta.est.noshrink[ind, , case, hyper.ind] = summary(tmp.noshrink, vars = c("beta"))[, 'Mean']
        
        beta.est.resid.noshrink[ind, , case, hyper.ind] = (beta.est.noshrink[ind, , case, hyper.ind] - beta)^2
        
      } else{
        bad.ind.noshrink = c(bad.ind.noshrink, ind)
      }
      print(ind)
    }
    
    
    
    
    cat(paste0("\n\n\nDone with Hyper ind ", hyper.ind))
  }
  
}


mean(colMeans(beta.est.resid.noshrink[, , 1, 1], na.rm = T))
mean(colMeans(beta.est.resid.noshrink[, , 2, 1], na.rm = T))

mean(colMeans(beta.est.resid.noshrink[, , 1, 2], na.rm = T))
mean(colMeans(beta.est.resid.noshrink[, , 2, 2], na.rm = T))

mean(colMeans(beta.est.resid.noshrink[, , 1, 3], na.rm = T))
mean(colMeans(beta.est.resid.noshrink[, , 2, 3], na.rm = T))













# Case 1, basic hyperparam ------------------------------------------------



hyper.param = expand.grid(c(0.01, 1, 5), c(0.1, 1, 10))

beta.est.full = beta.inprob.est.full = beta.best.full = beta.est.resid.full =
  alpha.est.full = alpha.inprob.est.full = alpha.best.full = alpha.est.resid.full =
  array(NA, dim = c(100, 10, 2))

beta.est.diag = beta.inprob.est.diag = beta.best.diag = alpha.est.diag = beta.est.resid.diag =
  alpha.inprob.est.diag = alpha.best.diag = alpha.est.resid = alpha.est.resid.diag = 
  array(NA, dim = c(100, 10, 2))
full.best = diag.best = array(NA, dim = c(100, 20, 2))


for (case in 1) {
  for (hyper.ind in 2) {
    bad.ind.full = bad.ind.diag = c()
    
    
    # Full Model --------------------------------------------------------------
    
    beta.est = beta.inprob.est = beta.best = matrix(NA, nrow = 100, ncol = 10)
    alpha.est = alpha.inprob.est = alpha.best = matrix(NA, nrow = 100, ncol = 10)
    # full.best = matrix(NA, nrow = 100, ncol = 20)
    for (ind in 1:100) {
      tmp.full = try(readRDS(
        file = paste0(
          "../Results/Sims/Full/Case",case,"/Poisson_full_hyper_",
          hyper.ind,
          "_ind_",
          ind,
          ".rds"
        )
      ))
      
      list2env(readRDS(paste0("../Data/Sims/Poisson_case_1_sim_", ind, ".rds")), globalenv())
      if (!("try-error" %in% class(tmp.full))) {
        ## Look at beta estimates
        beta.est.full[ind, , case] = summary(tmp.full, vars = c("beta"))[, 'Mean']
        
        beta.est.resid.full[ind, , case] = (beta.est.full[ind, , case] - beta)^2
        
        ## Look at prob of inclusion for beta
        beta.inprob.est.full[ind, , case] = summary(tmp.full, vars = c("latent_indicator1"))[, 'Mean']
        
        ## Find optimal model
        beta.in.full = as.mcmc(tmp.full, vars = "beta")
        beta.nz.full = data.frame(beta.in.full) %>% mutate(across(everything(), ~ if_else(. != 0, 1, 0)))

        ## Look at alpha estimates
        alpha.est.full[ind, , case] = summary(tmp.full, vars = c("lambda"))[, 'Mean']
        
        ## Find optimal model
        alpha.in.full = as.mcmc(tmp.full, vars = "lambda")
        alpha.nz.full = data.frame(alpha.in.full) %>% mutate(across(everything(), ~ if_else(. != 0, 1, 0)))

        ## Find full optimal model
        both.nz.collapse = apply(cbind(beta.nz.full, alpha.nz.full), 1, paste0, collapse = "")
        both.nz.collapse.tab = table(both.nz.collapse)
        digits = as.numeric(strsplit(as.character(names(sort(both.nz.collapse.tab, decreasing = T)[1])), "")[[1]])
        full.best[ind, , case] = digits

        
      } else{
        bad.ind.full = c(bad.ind.full, ind)
      }
      print(ind)
    }
    
    
    
    
    # Diag Model --------------------------------------------------------------
    
    
    for (ind in 1:100) {
      tmp.diag = try(readRDS(
        file = paste0(
          "../Results/Sims/Diag/Case",case,"/Poisson_diag_hyper_",
          hyper.ind,
          "_ind_",
          ind,
          ".rds"
        )
      ))
      
      list2env(readRDS(paste0("../Data/Sims/Poisson_case_1_sim_", ind, ".rds")), globalenv())
      if (!("try-error" %in% class(tmp.diag))) {
        ## Look at beta estimates
        beta.est.diag[ind, , case] = summary(tmp.diag, vars = c("beta"))[, 'Mean']
        beta.est.resid.diag[ind, , case] = (beta.est.diag[ind, , case] - beta)^2
        
        ## Look at prob of inclusion for beta
        beta.inprob.est.diag[ind, , case] = summary(tmp.diag, vars = c("latent_indicator1"))[, 'Mean']

        
        ## Find optimal model
        beta.in.diag = as.mcmc(tmp.diag, vars = "beta")
        beta.nz.diag = data.frame(beta.in.diag) %>% mutate(across(everything(), ~ if_else(. != 0, 1, 0)))
        
        beta.nz.collapse = apply(beta.nz.diag, 1, paste0, collapse = "")
        beta.nz.collapse.tab = table(beta.nz.collapse)
        digits = as.numeric(strsplit(as.character(names(sort(beta.nz.collapse.tab, decreasing = T)[1])), "")[[1]])
        
        beta.best.diag[ind, , case] = digits
        
        ## Look at alpha estimates
        alpha.est.diag[ind, , case] = summary(tmp.diag, vars = c("lambda"))[, 'Mean']
        
        ## Find optimal model
        alpha.in.diag = as.mcmc(tmp.diag, vars = "lambda")
        alpha.nz.diag = data.frame(alpha.in.diag) %>% mutate(across(everything(), ~ if_else(. != 0, 1, 0)))
        
        alpha.nz.collapse = apply(alpha.nz.diag, 1, paste0, collapse = "")
        alpha.nz.collapse.tab = table(alpha.nz.collapse)
        digits = as.numeric(strsplit(as.character(names(sort(alpha.nz.collapse.tab, decreasing = T)[1])), "")[[1]])
        alpha.best.diag[ind, , case] = digits
       
        ## Find diag optimal model
        both.nz.collapse = apply(cbind(beta.nz.diag, alpha.nz.diag), 1, paste0, collapse = "")
        both.nz.collapse.tab = table(both.nz.collapse)
        digits = as.numeric(strsplit(as.character(names(sort(both.nz.collapse.tab, decreasing = T)[1])), "")[[1]])
        diag.best[ind, , case] = digits
      
      } else{
        bad.ind.diag = c(bad.ind.diag, ind)
      }
      print(ind)
    }
    
    
    cat(paste0("\n\n\nDone with Hyper ind ", hyper.ind))
  }
  
}

full.best.collapse = apply(full.best[,,case], 1, paste0, collapse = "")
sort(table(full.best.collapse), decreasing = T)[1:5]

diag.best.collapse = apply(diag.best[,,case], 1, paste0, collapse = "")
sort(table(diag.best.collapse), decreasing = T)[1:5]


mean(colMeans(beta.est.resid.full[, , case], na.rm = T))
mean(colMeans(beta.est.resid.diag[, , case], na.rm = T))






# Case 2, basic hyperparam ----------------------------------------------------


hyper.param = expand.grid(c(0.01, 1, 5), c(0.1, 1, 10))

beta.est.full = beta.inprob.est.full = beta.best.full = beta.est.resid.full =
  alpha.est.full = alpha.inprob.est.full = alpha.best.full = alpha.est.resid.full =
  array(NA, dim = c(100, 10, 2))

beta.est.diag = beta.inprob.est.diag = beta.best.diag = alpha.est.diag = beta.est.resid.diag =
  alpha.inprob.est.diag = alpha.best.diag = alpha.est.resid = alpha.est.resid.diag = 
  array(NA, dim = c(100, 10, 2))
full.best = diag.best = array(NA, dim = c(100, 20, 2))


for (case in 2) {
  for (hyper.ind in 2) {
    bad.ind.full = bad.ind.diag = c()
    
    
    # Full Model --------------------------------------------------------------
    
    beta.est = beta.inprob.est = beta.best = matrix(NA, nrow = 100, ncol = 10)
    alpha.est = alpha.inprob.est = alpha.best = matrix(NA, nrow = 100, ncol = 10)
    # full.best = matrix(NA, nrow = 100, ncol = 20)
    for (ind in 1:100) {
      tmp.full = try(readRDS(
        file = paste0(
          "../Results/Sims/Full/Case",case,"/Poisson_full_hyper_",
          hyper.ind,
          "_ind_",
          ind,
          ".rds"
        )
      ))
      
      list2env(readRDS(paste0("../Data/Sims/Poisson_case_2_sim_", ind, ".rds")), globalenv())
      if (!("try-error" %in% class(tmp.full))) {
        ## Look at beta estimates
        beta.est.full[ind, , case] = summary(tmp.full, vars = c("beta"))[, 'Mean']
        
        beta.est.resid.full[ind, , case] = (beta.est.full[ind, , case] - beta)^2
        
        ## Look at prob of inclusion for beta
        beta.inprob.est.full[ind, , case] = summary(tmp.full, vars = c("latent_indicator1"))[, 'Mean']
        
        ## Find optimal model
        beta.in.full = as.mcmc(tmp.full, vars = "beta")
        beta.nz.full = data.frame(beta.in.full) %>% mutate(across(everything(), ~ if_else(. != 0, 1, 0)))

        ## Look at alpha estimates
        alpha.est.full[ind, , case] = summary(tmp.full, vars = c("lambda"))[, 'Mean']
        
        ## Find optimal model
        alpha.in.full = as.mcmc(tmp.full, vars = "lambda")
        alpha.nz.full = data.frame(alpha.in.full) %>% mutate(across(everything(), ~ if_else(. != 0, 1, 0)))

        ## Find full optimal model
        both.nz.collapse = apply(cbind(beta.nz.full, alpha.nz.full), 1, paste0, collapse = "")
        both.nz.collapse.tab = table(both.nz.collapse)
        digits = as.numeric(strsplit(as.character(names(sort(both.nz.collapse.tab, decreasing = T)[1])), "")[[1]])
        full.best[ind, , case] = digits

      } else{
        bad.ind.full = c(bad.ind.full, ind)
      }
      print(ind)
    }
    
    
    
    
    # Diag Model --------------------------------------------------------------
    
    
    for (ind in 1:100) {
      tmp.diag = try(readRDS(
        file = paste0(
          "../Results/Sims/Diag/Case",case,"/Poisson_diag_hyper_",
          hyper.ind,
          "_ind_",
          ind,
          ".rds"
        )
      ))
      list2env(readRDS(paste0("../Data/Sims/Poisson_case_2_sim_", ind, ".rds")), globalenv())
      if (!("try-error" %in% class(tmp.diag))) {
        ## Look at beta estimates
        beta.est.diag[ind, , case] = summary(tmp.diag, vars = c("beta"))[, 'Mean']
        beta.est.resid.diag[ind, , case] = (beta.est.diag[ind, , case] - beta)^2
        
        ## Look at prob of inclusion for beta
        beta.inprob.est.diag[ind, , case] = summary(tmp.diag, vars = c("latent_indicator1"))[, 'Mean']
        
        
        ## Find optimal model
        beta.in.diag = as.mcmc(tmp.diag, vars = "beta")
        beta.nz.diag = data.frame(beta.in.diag) %>% mutate(across(everything(), ~ if_else(. != 0, 1, 0)))
        
        beta.nz.collapse = apply(beta.nz.diag, 1, paste0, collapse = "")
        beta.nz.collapse.tab = table(beta.nz.collapse)
        digits = as.numeric(strsplit(as.character(names(sort(beta.nz.collapse.tab, decreasing = T)[1])), "")[[1]])
        
        beta.best.diag[ind, , case] = digits
        
        ## Look at alpha estimates
        alpha.est.diag[ind, , case] = summary(tmp.diag, vars = c("lambda"))[, 'Mean']
        
        ## Find optimal model
        alpha.in.diag = as.mcmc(tmp.diag, vars = "lambda")
        alpha.nz.diag = data.frame(alpha.in.diag) %>% mutate(across(everything(), ~ if_else(. != 0, 1, 0)))
        
        alpha.nz.collapse = apply(alpha.nz.diag, 1, paste0, collapse = "")
        alpha.nz.collapse.tab = table(alpha.nz.collapse)
        digits = as.numeric(strsplit(as.character(names(sort(alpha.nz.collapse.tab, decreasing = T)[1])), "")[[1]])
        alpha.best.diag[ind, , case] = digits

        ## Find diag optimal model
        both.nz.collapse = apply(cbind(beta.nz.diag, alpha.nz.diag), 1, paste0, collapse = "")
        both.nz.collapse.tab = table(both.nz.collapse)
        digits = as.numeric(strsplit(as.character(names(sort(both.nz.collapse.tab, decreasing = T)[1])), "")[[1]])
        diag.best[ind, , case] = digits

        
      } else{
        bad.ind.diag = c(bad.ind.diag, ind)
      }
      print(ind)
    }
    
    
    cat(paste0("\n\n\nDone with Hyper ind ", hyper.ind))
  }
  
}



full.best.collapse = apply(full.best[,,2], 1, paste0, collapse = "")
sort(table(full.best.collapse), decreasing = T)[1:5]

diag.best.collapse = apply(diag.best[,,2], 1, paste0, collapse = "")
sort(table(diag.best.collapse), decreasing = T)[1:5]



mean(colMeans(beta.est.resid.full[, , case], na.rm = T))
mean(colMeans(beta.est.resid.diag[, , case], na.rm = T))


























# Hyperparameter check ----------------------------------------------------




hyper.param = expand.grid(c(0.01, 1, 5), c(0.1, 1, 10))

beta.est.full = beta.inprob.est.full = beta.best.full = beta.est.resid.full =
  alpha.est.full = alpha.inprob.est.full = alpha.best.full = alpha.est.resid.full =
  array(NA, dim = c(100, 10, 2, nrow(hyper.param)))

full.best = array(NA, dim = c(100, 20, 2, nrow(hyper.param)))


for (case in 1:2) {
  for (hyper.ind in 1:nrow(hyper.param)) {
    bad.ind.full = bad.ind.diag = c()
    
    
    # Full Model --------------------------------------------------------------
    
    beta.est = beta.inprob.est = beta.best = matrix(NA, nrow = 100, ncol = 10)
    alpha.est = alpha.inprob.est = alpha.best = matrix(NA, nrow = 100, ncol = 10)
    for (ind in 1:100) {
      tmp.full = try(readRDS(
        file = paste0(
          "../Results/Sims/Full/Case",case,"/Poisson_full_hyper_",
          hyper.ind,
          "_ind_",
          ind,
          ".rds"
        )
      ))
      
      list2env(readRDS(paste0("../Data/Sims/Poisson_case_1_sim_", ind, ".rds")), globalenv())
      if (!("try-error" %in% class(tmp.full))) {
        ## Look at beta estimates
        beta.est.full[ind, , case, hyper.ind] = summary(tmp.full, vars = c("beta"))[, 'Mean']
        
        beta.est.resid.full[ind, , case, hyper.ind] = (beta.est.full[ind, , case, hyper.ind] - beta)^2
        
        ## Look at prob of inclusion for beta
        beta.inprob.est.full[ind, , case, hyper.ind] = summary(tmp.full, vars = c("latent_indicator1"))[, 'Mean']
        
        ## Find optimal model
        beta.in.full = as.mcmc(tmp.full, vars = "beta")
        beta.nz.full = data.frame(beta.in.full) %>% mutate(across(everything(), ~ if_else(. != 0, 1, 0)))
        
        ## Look at alpha estimates
        alpha.est.full[ind, , case, hyper.ind] = summary(tmp.full, vars = c("lambda"))[, 'Mean']
        
        ## Find optimal model
        alpha.in.full = as.mcmc(tmp.full, vars = "lambda")
        alpha.nz.full = data.frame(alpha.in.full) %>% mutate(across(everything(), ~ if_else(. != 0, 1, 0)))
        
        ## Find full optimal model
        both.nz.collapse = apply(cbind(beta.nz.full, alpha.nz.full), 1, paste0, collapse = "")
        both.nz.collapse.tab = table(both.nz.collapse)
        digits = as.numeric(strsplit(as.character(names(sort(both.nz.collapse.tab, decreasing = T)[1])), "")[[1]])
        full.best[ind, , case, hyper.ind] = digits

        
      } else{
        bad.ind.full = c(bad.ind.full, ind)
      }
      print(ind)
    }
    
    
    cat(paste0("\n\n\nDone with Hyper ind ", hyper.ind))
  }
  
}

for(hyper.ind in 1:nrow(hyper.param)){
  full.best.collapse = apply(full.best[,,1, hyper.ind], 1, paste0, collapse = "")
  print(cbind(hyper.param[hyper.ind,], sort(table(full.best.collapse), decreasing = T)[1]))
  # print(sort(table(full.best.collapse), decreasing = T)[1])
}

for(hyper.ind in 1:nrow(hyper.param)){
  print(cbind(hyper.param[hyper.ind,], round(mean(colMeans(beta.est.resid.full[, , 1, hyper.ind], na.rm = T)), digits = 6)))
  # print(sort(table(full.best.collapse), decreasing = T)[1])
}


full.best.collapse = apply(full.best[,,case], 1, paste0, collapse = "")
sort(table(full.best.collapse), decreasing = T)[1:5]

diag.best.collapse = apply(diag.best[,,case], 1, paste0, collapse = "")
sort(table(diag.best.collapse), decreasing = T)[1:5]


mean(colMeans(beta.est.resid.full[, , case], na.rm = T))
mean(colMeans(beta.est.resid.diag[, , case], na.rm = T))




















