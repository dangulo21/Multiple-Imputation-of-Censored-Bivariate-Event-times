######################################################################################
######################################################################################
############################ APPLICATION - HPV/AIDS ##################################
######################################################################################
######################################################################################
library(dplyr)
library(tidyr)
library(survival)
library(mhazard)
library(MASS)
library(Rfast)
library(mgcv)
library(survSpearman)


  source("../scripts/00_helpers.R")
  
  ######################################################################################
  ################################## READ DATA #########################################
  ######################################################################################
  
  df <- survSpearman::CCASAnetData
  
  # Define tau
  tau = 15
 
  # Create ta-restricted times-to-event
  df <- 
    df |>
    mutate(
      timeToRegimenChange = pmin(timeToRegimenChange, tau),
      timeToViralFailure = pmin(timeToViralFailure, tau),
      regimenChange = ifelse(timeToRegimenChange >= tau, 1, regimenChange),
      viralFailure = ifelse(timeToViralFailure >= tau, 1, viralFailure))
  
  data_sim <- df |> dplyr::select(-site)
  n = nrow(data_sim)
  
  # Organize the data
  colnames(data_sim) <- c("x1","delta1","x2","delta2")
  
  #Define the burn-in parameters, L, and number completed data sets M, for IT-MI
  L = 5
  M = 5
  
  # Generate a common grid of values to evaluate estimates probabilites of IT-MI and Dabrowska's estimator
  
  df_comp <- 
    data.frame(
      expand_grid(t1 = seq(0.0,tau-0.5,0.5), t2 = seq(0.0,tau-0.5,0.5)))

  
  ######################################################################################
  ############################## MULTIPLE IMPUTATION ###################################
  ######################################################################################
  
  df_res <- data.frame()
  
  #create imputation vectors
  df_imp <-
    data_sim |>
    mutate(
      impute.T1 = case_when(
        delta1 == T ~ x1,
        delta1 == F ~ NA),
      impute.T2 = case_when(
        delta2 == T ~ x2,
        delta2 == F ~ NA))
  
  ############################################################################################
  ############################ Initialization Stage ##########################################
  ############################################################################################
  
  # Dimension to impute first (the one with more observed events)
  time.first <- ifelse(sum(df_imp$delta1) > sum(df_imp$delta2), "x1", "x2")
  delta.first <- ifelse(sum(df_imp$delta1) > sum(df_imp$delta2), "delta1", "delta2")
  impute.first <- ifelse(sum(df_imp$delta1) > sum(df_imp$delta2), "impute.T1", "impute.T2") 
  
  
  for(i in which(!df_imp[,delta.first])){
 
    # Risk set individuals with restriction of being at risk
    risk.set <- which(df_imp[,time.first] > df_imp[i,time.first])
    
    # Variables related to the first dimension to impute for those in the risk set
    df_risk_set <- df_imp[risk.set, c(time.first, delta.first)]
    
    # Estimate the KM among those in the risk set.
    # Note that the KM, by construction of the sample, will go to zero
    # The last observed time (which will be in the risk set always), is an observed time
    km <- survfit(Surv(eval(parse(text = time.first)), eval(parse(text = delta.first))) ~ 1, data = df_risk_set)
    km <- summary(km)
    
    # Generate an uniform distribution
    u <- runif(1)
    
    # Inverse transform theorem
    impute.index = which.max(km$surv <= u)
    
    # Assign the imputation
    df_imp[i,impute.first] = km$time[impute.index]
    
  }
  
  
  ############################################################################################
  ############################ NON PARAMETRIC GIBBS SAMPLER ##################################
  ############################################################################################
  
  # (L+M)*2 is because in each iteration we are imputing one event at a time
  
  iter = 0
  df_mi <- data.frame()
  
  while(iter < (L+M)*2){
    iter = iter + 1
    
    # Select dimension to impute based one dimension imputed previously
    time.second <- ifelse(time.first == "x1","x2","x1")
    delta.second <- ifelse(delta.first == "delta1","delta2","delta1")
    impute.second <- ifelse(impute.first == "impute.T1","impute.T2","impute.T1")
    
    # Clean dimension to impute
    # For those who need an imputation, their impute.second variable should be NA. 
    df_imp[which(!df_imp[,delta.second]),impute.second] <- NA
    
    # Imputation based on restricted marginal (using two restrictions)
    for(i in which(!df_imp[,delta.second])){
      res_marg = restricted_marginal(i, time.second, delta.second, impute.first, delta.first, df_imp)
      df_imp[i,impute.second] = res_marg$impute
    }
    
    # Once the iteration is complete, we switch the labels such that the next iteration imputes 
    # the remaining time-to-event
    
    time.first = time.second
    delta.first = delta.second
    impute.first = impute.second
    
    ####################################################################################
    ################### Analysis of Augmented Dataset Stage ############################
    ####################################################################################
    
    # Save the results after the burning period is concluded: L*2
    # Since each iteration in the outer loop is for one time-to-event at a time
    # We only save the iteration after we have updated both time-to-events
    
    if(iter > L*2 & iter %% 2 != 0){
      
      df_iter <- df_comp
      
      df_iter$imp = iter
      
      # For that completed data set (df_imp), get:
      # i. the estimated probabilities across the grid of values (df_iter).
      
        for(i in 1:nrow(df_iter)){
          t1 = df_iter[i,]$t1
          t2 = df_iter[i,]$t2
          
          # estimated bivariate survival probabilities
          df_iter$shat_T12[i] = nrow(df_imp[df_imp$impute.T1 >= t1 & df_imp$impute.T2 >= t2,])/n
          df_iter$shat_T1[i] = nrow(df_imp[df_imp$impute.T1 >= t1,])/n
          df_iter$shat_T2[i] = nrow(df_imp[df_imp$impute.T2 >= t2,])/n}
        
        # With the estimated probabilities, get variability estimates
        df_iter$var_T12 = df_iter$shat_T12*(1-df_iter$shat_T12)/n
        df_iter$var_T1 = df_iter$shat_T1*(1-df_iter$shat_T1)/n
        df_iter$var_T2 = df_iter$shat_T2*(1-df_iter$shat_T2)/n
        
        # Complementary log-log transformation
        # Check: Prostate cancer: net survival and cause-specific survival rates after multiple imputation
        df_iter$log_T12 = log(-log(1-df_iter$shat_T12)) 
        df_iter$var_log_T12 = df_iter$var_T12/((log(1-df_iter$shat_T12)*(1-df_iter$shat_T12))^2)
        #
        df_iter$log_T1 = log(-log(1-df_iter$shat_T1))
        df_iter$var_log_T1 = df_iter$var_T1/((log(1-df_iter$shat_T1)*(1-df_iter$shat_T1))^2)
        
        df_iter$log_T2 = log(-log(1-df_iter$shat_T2))
        df_iter$var_log_T2 = df_iter$var_T2/((log(1-df_iter$shat_T2)*(1-df_iter$shat_T2))^2)
        
        # save i. 
        df_mi <- rbind(df_mi,df_iter)
    }
  }
  
  ####################################################################################
  ####################### Multiple Imputation Calculations ###########################
  ####################################################################################
  
  df_mi <- 
    df_mi |>
    group_by(t1, t2) |>
    summarise(
      # MI across the M imputed datasets
      # Bivariate
      shat_T12_bar = mean(shat_T12),
      u_T12 = mean(var_T12),
      b_T12 = var(shat_T12),
      
      # log(log(1-S(t_1,t_2)))
      log_T12_bar = mean(log_T12),
      u_log_T12 = mean(var_log_T12),
      b_log_T12 = var(log_T12),
      
      # First Marginal
      shat_T1_bar = mean(shat_T1),
      u_T1 = mean(var_T1),
      b_T1 = var(shat_T1),
      # log(log(1-S(t_1,0)))
      log_T1_bar = mean(log_T1),
      u_log_T1 = mean(var_log_T1),
      b_log_T1 = var(log_T1),
      
      # Second Marginal
      shat_T2_bar = mean(shat_T2),
      u_T2 = mean(var_T2),
      b_T2 = var(shat_T2),
      # log(log(1-S(0,t_2)))
      log_T2_bar = mean(log_T2),
      u_log_T2 = mean(var_log_T2),
      b_log_T2 = var(log_T2)) |>
    rename(
      # Rename so it is easier
      shat_T12 = shat_T12_bar, log_T12 = log_T12_bar,
      shat_T1 = shat_T1_bar, log_T1 = log_T1_bar,
      shat_T2 = shat_T2_bar, log_T2 = log_T2_bar) |>
    mutate(
      # Auxiliary variables for CI estimation
      # Bivariate
      se_T12 = sqrt(u_T12 + b_T12*(1 + 1/M)),
      v_T12 = (M-1)*((1+(M/(M+1))*u_T12/b_T12)^2),
      se_log_T12 = sqrt(u_log_T12 + b_log_T12*(1 + 1/M)),
      v_log_T12 = (M-1)*((1+(M/(M+1))*u_log_T12/b_log_T12)^2),
      
      # First Marginal
      se_T1 = sqrt(u_T1+b_T1*(1+1/M)),
      v_T1 = (M-1)*((1+(M/(M+1))*u_T1/b_T1)^2),
      se_log_T1 = sqrt(u_log_T1 + b_log_T1*(1+1/M)),
      v_log_T1 = (M-1)*((1+(M/(M+1))*u_log_T1/b_log_T1)^2),
      
      # Second Marginal
      se_T2 = sqrt(u_T2+b_T2*(1+1/M)),
      v_T2 = (M-1)*((1+(M/(M+1))*u_T2/b_T2)^2),
      se_log_T2 = sqrt(u_log_T2+b_log_T2*(1+1/M)),
      v_log_T2 = (M-1)*((1+(M/(M+1))*u_log_T2/b_log_T2)^2)) |>
    mutate(
      # Bias Confidence interval estimation 
      
      # Bivariate
      LB_T12 = shat_T12 - qt(0.975,v_T12)*se_T12,
      UB_T12 = shat_T12 + qt(0.975,v_T12)*se_T12,
      length_T12 = UB_T12 - LB_T12,
      log_LB_T12 = 1 - exp(-exp(log_T12 - qt(0.975,v_log_T12)*se_log_T12)),
      log_UB_T12 = 1 - exp(-exp(log_T12 + qt(0.975,v_log_T12)*se_log_T12)),
      length_log_T12 = log_UB_T12-log_LB_T12,
      
      
      # Marginal T1
      LB_T1 = shat_T1 - qt(0.975,v_T1)*se_T1,
      UB_T1 = shat_T1 + qt(0.975,v_T1)*se_T1,
      length_T1 = UB_T1 - LB_T1,
      log_LB_T1 = 1 - exp(-exp(log_T1 - qt(0.975,v_log_T1)*se_log_T1)),
      log_UB_T1 = 1 - exp(-exp(log_T1 + qt(0.975,v_log_T1)*se_log_T1)),
      length_log_T1 = log_UB_T1-log_LB_T1,
      
      
      # Marginal T2
      LB_T2 = shat_T2 - qt(0.975,v_T2)*se_T2,
      UB_T2 = shat_T2 + qt(0.975,v_T2)*se_T2,
      length_T2 = UB_T2 - LB_T2,
      log_LB_T2 = 1 - exp(-exp(log_T2 - qt(0.975,v_log_T2)*se_log_T2)),
      log_UB_T2 = 1 - exp(-exp(log_T2 + qt(0.975,v_log_T2)*se_log_T2)),
      length_log_T2 = log_UB_T2 - log_LB_T2) |>
    mutate(
      # Simulation configuration
      n = n,
      method = "IT-MI")
  
  df_mi <-
    df_mi |>
    dplyr::select(n:method,
                  t1:shat_T12, se_T12, LB_T12:length_log_T12,
                  shat_T1, se_T1, length_T1, length_log_T1,
                  shat_T2, se_T2, length_T2, length_log_T2)
  
  ####################################################################################################
  ####################################################################################################
  ########################################## DABROWSKA ###############################################
  ####################################################################################################
  ####################################################################################################
  
  
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
    Y2 <- data_sim[index,3]
    Delta1 <- data_sim[index,2]
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
  
  
  # df_dab has the order: all t1, with t2= 0, followed by all t1, with t2 = 1, and so on
  # Dabrowska's mhzard functions returns a matrix with nrow = t1 and ncol = t2
  # With an added zero at the beginning Therefore c(Fhat) reads the survival probabilities
  # at: all t1, with t2 = 0, followed by all t1, with t2 = 1, and son on. 
  # shat.boot$t is a matrix of 1000 x (t1xt2). Where the order of the columns
  # is all t1, with t2 = 0, followed by all r1, with t2 = 1 and so on (why? check return in function). 
  
  
  df_dab <- 
    data.frame(
      expand.grid(t1 = c(0,ft.1), t2 = c(0,ft.2)),
      s.hat = c(Dabrowska$Fhat),
      lb = CI.boot[1,],
      ub = CI.boot[2,],
      sd.boot = sd.boot)
  
  
  for(i in 1:nrow(df_comp)){
    
    val.t1 = max(c(0,ft.1)[c(0,ft.1) <= df_comp$t1[i]])
    val.t2 = max(c(0,ft.2)[c(0,ft.2) <= df_comp$t2[i]])
    
    # Bivariate
    df_comp$shat_T12[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == val.t2,]$s.hat
    df_comp$se_T12[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == val.t2,]$sd.boot
    df_comp$LB_T12[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == val.t2,]$lb
    df_comp$UB_T12[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == val.t2,]$ub
    
    # Marginal T1
    df_comp$shat_T1[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == 0,]$s.hat
    df_comp$se_T1[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == 0,]$sd.boot
    df_comp$LB_T1[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == 0,]$lb
    df_comp$UB_T1[i] = df_dab[df_dab$t1 == val.t1 & df_dab$t2 == 0,]$ub
    
    # Marginal T2
    df_comp$shat_T2[i] = df_dab[df_dab$t1 == 0 & df_dab$t2 == val.t2,]$s.hat
    df_comp$se_T2[i] = df_dab[df_dab$t1 == 0 & df_dab$t2 == val.t2,]$sd.boot
    df_comp$LB_T2[i] = df_dab[df_dab$t1 == 0 & df_dab$t2 == val.t2,]$lb
    df_comp$UB_T2[i] = df_dab[df_dab$t1 == 0 & df_dab$t2 == val.t2,]$ub
  }
  
  
  df_dab <-
    df_comp |>
    mutate(
      # Bias Confidence interval estimation 
      
      # Bivariate
      shat_T12_log = ifelse(shat_T12 == 0, 0.001, shat_T12),
      length_T12 = UB_T12 - LB_T12,
      log_LB_T12 = shat_T12_log^exp(-qnorm(0.975)*(se_T12/(log(shat_T12_log)*shat_T12_log))),
      log_UB_T12 = shat_T12_log^exp(+qnorm(0.975)*(se_T12/(log(shat_T12_log)*shat_T12_log))),
      length_log_T12 = log_UB_T12-log_LB_T12,
      
      # Marginal T1
      shat_T1_log = ifelse(shat_T1 == 0, 0.001, shat_T1),
      length_T1 = UB_T1 - LB_T1,
      log_LB_T1 = shat_T1_log^exp(-qnorm(0.975)*(se_T1/(log(shat_T1_log)*shat_T1_log))),
      log_UB_T1 = shat_T1_log^exp(+qnorm(0.975)*(se_T1/(log(shat_T1_log)*shat_T1_log))),
      length_log_T1 = log_UB_T1-log_LB_T1,
      
      
      # Marginal T2
      shat_T2_log = ifelse(shat_T2 == 0, 0.001, shat_T2),
      length_T2 = UB_T2 - LB_T2,
      log_LB_T2 = shat_T2_log^exp(-qnorm(0.975)*(se_T2/(log(shat_T2_log)*shat_T2_log))),
      log_UB_T2 = shat_T2_log^exp(+qnorm(0.975)*(se_T2/(log(shat_T2_log)*shat_T2_log))),
      length_log_T2 = log_UB_T2-log_LB_T2) |>
    mutate(
      # Simulation configuration
      n = n, 
      method = "Dabrowska")
  
  
  
  df_dab <-
    df_dab |>
    dplyr::select(n:method,
                  t1:shat_T12, se_T12, LB_T12:length_log_T12,
                  shat_T1, se_T1, length_T1, length_log_T1,
                  shat_T2, se_T2, length_T2, length_log_T2) |>
    dplyr::select(-LB_T1,-UB_T1,-LB_T2,-UB_T2)
  
  
  ####################################################################################
  #################################### SAVE ##########################################
  ####################################################################################
  
  res <- rbind(df_mi,df_dab)
  name = paste0("../results/application/application.RData")
  save(res, file = name)
   






