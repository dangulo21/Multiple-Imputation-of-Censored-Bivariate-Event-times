#####################################################################################################
#####################################################################################################
##################### Imputation of Bivariate Censored Outcomes #####################################
#####################################################################################################
#####################################################################################################
library(dplyr)
library(survival)
library(mhazard)
library(MASS)
library(Rfast)
library(mgcv)
  
  ##################################################################################################
  ################################## READ DATA #####################################################
  ##################################################################################################
  
  method = "exponential-uniform"
  
  rho = 0.5
  censoring = 50
  
  # Load the real bivariate time-to-event and censored data.
  # Change the directory to where the data is located. 
  load(paste0("../data/",method,"/real_data_corr_",rho,"_cens_",censoring,".RData"))
  
  # Load the true probabilities of bivariate time-to-event, across a grid of (t1,t2) values.
  # Change the directory to where the data is located. 
  load(paste0("../data/",method,"/true_probs_corr_",rho,".RData"))
  
  ##################################################################################################
  ############################### SIMULATION SET UP ################################################
  ##################################################################################################
  
  # Load the simulation parameters, such as N and tau (for tau restricted-event-times).
  # Change the directory to where the data is located.
  load(paste0("../data/",method,"/parameters.RData"))
  
  #Set the sample size
  n = 100
  
  # Checking for simulation seeds
  sim = 487
  set.seed(sim + 19700620)
  
  data_sim <- sample(seq(1,N), n, replace = FALSE)
  data_sim <- real_data[data_sim,] |> dplyr::select(x1,x2,delta1,delta2)

  ####################################################################################
  #################################### Estimation ####################################
  ####################################################################################
  
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
        estimator = "dabrowska")
    
    index.i <- sapply(ft.1, function(x1, x2) {sum(x2<=x1)}, x2 = ft1)
    index.j <- sapply(ft.2, function(x1, x2) {sum(x2<=x1)}, x2 = ft2)
    
    #Note: c(Matrix) returns a vector reading the matrix each column at a time. 
    return(c(dabrowska$Fhat[c(1,index.i+1),c(1,index.j+1)]))
  }
  
  shat.boot <- boot::boot(data_sim, statistic = boot.S, R = 1000)
  
  #Getting the bootstrap results 
  CI.boot <- apply(shat.boot$t, 2, function(x) quantile(x, c(0.025,0.975))) 
  sd.boot <- apply(shat.boot$t, 2, sd)
  
  
  # Dabrowska's mhzard functions returns a matrix with nrow = t1 and ncol = t2
  # With an added zero at the beginning. Therefore c(Fhat) reads the survival probabilities
  # at: all t1, with t2 = 0, followed by all t1, with t2 = 1, and son on. 
  # shat.boot$t is a matrix of [R = 1000 x (t1xt2)]. Where the order of the columns
  # is all t1, with t2 = 0, followed by all r1, with t2 = 1 and so on. 
  
  
  # df_dab has the order: all t1, with t2= 0, followed by all t1, with t2 = 1, and so on
  df_dab <- 
    data.frame(
      expand.grid(t1 = c(0,ft.1), t2 = c(0,ft.2)),
      s.hat = c(Dabrowska$Fhat),
      lb = CI.boot[1,],
      ub = CI.boot[2,],
      sd.boot = sd.boot)
  
  # Mapping survival probabilities to the t1xt2 grid of values in df_comp
  for(i in 1:nrow(df_comp)){
    
    val.t1 = max(c(0,ft.1)[c(0,ft.1) <= df_comp$t1[i]])
    val.t2 = max(c(0,ft.2)[c(0,ft.2) <= df_comp$t2[i]])
    
    #Bivariate
    df_comp$shat_T12[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == val.t2,]$s.hat
    df_comp$se_T12[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == val.t2,]$sd.boot
    df_comp$LB_T12[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == val.t2,]$lb
    df_comp$UB_T12[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == val.t2,]$ub
    
    
    #Marginal T1
    df_comp$shat_T1[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == 0,]$s.hat
    df_comp$se_T1[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == 0,]$sd.boot
    df_comp$LB_T1[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == 0,]$lb
    df_comp$UB_T1[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == 0,]$ub
    
    
    #Marginal T2
    df_comp$shat_T2[i] = df_dab[df_dab$t1 == 0 & df_dab$t2 == val.t2,]$s.hat
    df_comp$se_T2[i] = df_dab[df_dab$t1 == 0 & df_dab$t2 == val.t2,]$sd.boot
    df_comp$LB_T2[i] = df_dab[df_dab$t1 == 0 & df_dab$t2 == val.t2,]$lb
    df_comp$UB_T2[i] = df_dab[df_dab$t1 == 0 & df_dab$t2 == val.t2,]$ub}
  

  # Estimating c-loglog CP and length of 95% CI.                 
  df_dab <-
    df_comp |>
    mutate(
      # Bias Confidence interval estimation 
      
      # Bivariate
      shat_T12_log = ifelse(shat_T12 == 0, 0.001, shat_T12),
      bias_T12 = true_prob - shat_T12,
      coverage_T12 = (true_prob >= LB_T12)*(true_prob <= UB_T12),
      length_T12 = UB_T12 - LB_T12,
      log_LB_T12 = shat_T12_log^exp(-qnorm(0.975)*(se_T12/(log(shat_T12_log)*shat_T12_log))),
      log_UB_T12 = shat_T12_log^exp(+qnorm(0.975)*(se_T12/(log(shat_T12_log)*shat_T12_log))),
      coverage_log_T12 = (true_prob >= log_LB_T12)*(true_prob <= log_UB_T12),
      length_log_T12 = log_UB_T12-log_LB_T12,
      
      # Marginal T1
      shat_T1_log = ifelse(shat_T1 == 0, 0.001, shat_T1),
      bias_T1 = true_prob_T1 - shat_T1,
      coverage_T1 = (true_prob_T1 >= LB_T1)*(true_prob_T1 <= UB_T1),
      length_T1 = UB_T1 - LB_T1,
      log_LB_T1 = shat_T1_log^exp(-qnorm(0.975)*(se_T1/(log(shat_T1_log)*shat_T1_log))),
      log_UB_T1 = shat_T1_log^exp(+qnorm(0.975)*(se_T1/(log(shat_T1_log)*shat_T1_log))),
      coverage_log_T1 = (true_prob_T1 >= log_LB_T1)*(true_prob_T1 <= log_UB_T1),
      length_log_T1 = log_UB_T1-log_LB_T1,
      
      
      # Marginal T2
      shat_T2_log = ifelse(shat_T2 == 0, 0.001, shat_T2),
      bias_T2 = true_prob_T2 - shat_T2,
      coverage_T2 = (true_prob_T2 >= LB_T2)*(true_prob_T2 <= UB_T2),
      length_T2 = UB_T2 - LB_T2,
      log_LB_T2 = shat_T2_log^exp(-qnorm(0.975)*(se_T2/(log(shat_T2_log)*shat_T2_log))),
      log_UB_T2 = shat_T2_log^exp(+qnorm(0.975)*(se_T2/(log(shat_T2_log)*shat_T2_log))),
      coverage_log_T2 = (true_prob_T2 >= log_LB_T2)*(true_prob_T2 <= log_UB_T2),
      length_log_T2 = log_UB_T2-log_LB_T2) |>
    mutate(
      # Simulation configuration
      correlation = rho,
      censoring = censoring,
      n = n, 
      sim = sim,
      method = "Dabrowska")
  
  df_dab <-
    df_dab |>
    dplyr::select(correlation:method,
                  t1:shat_T12, bias_T12, se_T12, coverage_T12, length_T12, coverage_log_T12, length_log_T12,
                  shat_T1, bias_T1, se_T1, coverage_T1, length_T1, coverage_log_T1, length_log_T1,
                  shat_T2, bias_T2, se_T2, coverage_T2, length_T2, coverage_log_T2, length_log_T2)
  
  
  ####################################################################################
  #################################### SAVE ##########################################
  ####################################################################################

  name_all = paste0("rho_",rho,"_cens_",censoring,"_n_",n,"_sim_",sim,".RData")
  name = paste0("MI/results-all/",name_all)
  save(df_dab, file = name)
  


  





