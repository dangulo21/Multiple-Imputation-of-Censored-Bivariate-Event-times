######################################################################################
######################################################################################
##################### Imputation of Bivariate Censored Outcomes ######################
######################################################################################
######################################################################################

library(dplyr)
library(survival)
library(mhazard)
library(MASS)
library(Rfast)
library(mgcv)


  ######################################################################################
  ################################## READ DATA #########################################
  ######################################################################################
  
  method = "exponential"
  
  rho = 0.5
  censoring = 40
  
  path_data = paste0("data/",method,"/real_data_corr_",rho,"_cens_",censoring,".RData")
  load(path_data)
  
  load(paste0("data/",method,"/true_probs_corr_",rho,".RData"))

  #######################################################################################
  ############################### SIMULATION SET UP #####################################
  #######################################################################################
  
  load(paste0("data/",method,"/parameters.RData"))
  
  n = 250

  sim = 121
  set.seed(sim+19700620)
  
  data_sim <- sample(seq(1,N), n, replace = FALSE)
  data_sim <- real_data[data_sim,] |> dplyr::select(x1,x2,delta1,delta2)
  remove(real_data)
  
  #########################################################################################
  ############################## Dabrowska and Variance ###################################
  #########################################################################################
  
  Dabrowska <- 
    npSurv2(
      Y1 = data_sim$x1,Y2 = data_sim$x2, 
      Delta1 = data_sim$delta1, Delta2 = data_sim$delta2,
      estimator = "dabrowska",
      conf.int = TRUE,
      R = 1000)
  
  
  ft.1 <- sort(unique(data_sim[data_sim$delta1 == TRUE,]$x1))
  ft.2 <- sort(unique(data_sim[data_sim$delta2 == TRUE,]$x2))
  
  
  boot.S <- function(data_sim,index){
    
    Y1 <- data_sim[index,1]
    Y2 <- data_sim[index,2]
    Delta1 <- data_sim[index,3]
    Delta2 <- data_sim[index,4]
    ft1 <- unique(sort(Y1[Delta1==T]))
    ft2 <- unique(sort(Y2[Delta2==T]))
    
    dabrowska <- 
      npSurv2(
        Y1,Y2,Delta1, Delta2,
        estimator = "dabrowska"
      )
    
    index.i <- sapply(ft.1, function(x1, x2) {sum(x2<=x1)}, x2 = ft1)
    index.j <- sapply(ft.2, function(x1, x2) {sum(x2<=x1)}, x2 = ft2)
    
    return(dabrowska$Fhat[c(1,index.i+1),c(1,index.j+1)])
  }
  
  shat.boot <- boot::boot(data_sim, statistic = boot.S, R = 1000)
  
  #Getting the bootstrap reults 
  quants = c(0.025,0.975)
  CI.boot <- apply(shat.boot$t, 2, function(x) quantile(x,quants)) 
  sd.boot <- apply(shat.boot$t, 2, sd)
  
  
  df <- 
    data.frame(
      expand.grid(t1 = c(0,ft.1), t2 = c(0,ft.2)),
      s.hat = c(Dabrowska$Fhat),
      lb = CI.boot[1,],
      ub = CI.boot[2,],
      sd.boot = sd.boot)
  
  
  
  for(i in 1:nrow(df_comp)){
    
    val.t1 = max(c(0,ft.1)[c(0,ft.1) <= df_comp$t1[i]])
    val.t2 = max(c(0,ft.2)[c(0,ft.2) <= df_comp$t2[i]])
    
    df_comp$prob.T12[i] = df[df$t1 == val.t1 & df$t2 == val.t2,]$s.hat
    df_comp$se.T12[i] = df[df$t1 == val.t1 & df$t2 == val.t2,]$sd.boot
    df_comp$LB.T12[i] = df[df$t1 == val.t1 & df$t2 == val.t2,]$lb
    df_comp$UB.T12[i] = df[df$t1 == val.t1 & df$t2 == val.t2,]$ub

    
    df_comp$prob.T1[i] = df[df$t1 == val.t1 & df$t2 == 0,]$s.hat
    df_comp$se.T1[i] = df[df$t1 == val.t1 & df$t2 == 0,]$sd.boot
    df_comp$LB.T1[i] = df[df$t1 == val.t1 & df$t2 == 0,]$lb
    df_comp$UB.T1[i] = df[df$t1 == val.t1 & df$t2 == 0,]$ub
    
    
    df_comp$prob_T2[i] = df[df$t1 == 0 & df$t2 == val.t2,]$s.hat
    df_comp$se.T2[i] = df[df$t1 == 0 & df$t2 == val.t2,]$sd.boot
    df_comp$LB.T2[i] = df[df$t1 == 0 & df$t2 == val.t2,]$lb
    df_comp$UB.T2[i] = df[df$t1 == 0 & df$t2 == val.t2,]$ub
  }
  
  
  df_res <-
    df_comp |>
    mutate(n = n, sim = sim)
  
  ####################################################################################
  #################################### SAVE ##########################################
  ####################################################################################
  
  name_all = paste0("rho_",rho,"_cens_",censoring,"_n_",n,"_sim_",sim,".RData")
  name = paste0("results-all/",name_all)
  save(df_res, file = name)
  
  












