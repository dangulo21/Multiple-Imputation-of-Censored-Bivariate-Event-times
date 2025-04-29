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
library(MIWilson)


  source("../scripts/00_helpers.R")

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
  
  #Define the burn-in parameters, L, and number completed data sets M
  L = 2
  M = 5

  #Set the sample size
  n = 100
  
  # Checking for simulation seeds
  sim = 487
  set.seed(sim + 19700620)
  
  data_sim <- sample(seq(1,N), n, replace = FALSE)
  data_sim <- real_data[data_sim,] |> dplyr::select(x1,x2,delta1,delta2)
  
  # The maximum observed time should be an event (to make sure we have correct imputations, if not, change seeds)
  max.t1.delta = ifelse(TRUE %in% data_sim[which(data_sim$x1 == max(data_sim$x1)),]$delta1, TRUE, FALSE)
  max.t2.delta = ifelse(TRUE %in% data_sim[which(data_sim$x2 == max(data_sim$x2)),]$delta2, TRUE, FALSE)
  
  ####################################################################################################
  ####################################################################################################
  ########################################## IT - MI #################################################
  ####################################################################################################
  ####################################################################################################
  
  # Allocating the results of the MI procedure
  df_mi <- data.frame()
  
  # Create imputation vectors
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
  
  
  # Marginal Imputation
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
  avg.epsilon = NA
  
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
      avg.epsilon = c(avg.epsilon, res_marg$epsilon)
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
    
    k = qnorm(0.975)
    
    if(iter > L*2 & iter %% 2 != 0){
       
      df_iter <- df_comp
      df_iter$imp = iter
      
      # For that completed data set (df_imp), get the estimated probabilities across the grid of values. 
      for(i in 1:nrow(df_iter)){
        t1 = df_iter[i,]$t1
        t2 = df_iter[i,]$t2
        
        df_iter$shat_T12[i] = nrow(df_imp[df_imp$impute.T1 >= t1 & df_imp$impute.T2 >= t2,])/n
        df_iter$shat_T1[i] = nrow(df_imp[df_imp$impute.T1 >= t1,])/n
        df_iter$shat_T2[i] = nrow(df_imp[df_imp$impute.T2 >= t2,])/n
        
        df_iter$shat_T12_ac[i] = (nrow(df_imp[df_imp$impute.T1 >= t1 & df_imp$impute.T2 >= t2,]) + k^2/2)/(n+k^2)}
      
      # With the estimated probabilities, get variability estimates
      df_iter$var_T12 = df_iter$shat_T12*(1-df_iter$shat_T12)/n
      df_iter$var_T1 = df_iter$shat_T1*(1-df_iter$shat_T1)/n
      df_iter$var_T2 = df_iter$shat_T2*(1-df_iter$shat_T2)/n
      
      # Complementary log-log transformation
      # Check: Prostate cancer: net survival and cause-specific survival rates after multiple imputation
      df_iter$shat_T12_log = ifelse(df_iter$shat_T12 <= 0.0001, 0.0001, df_iter$shat_T12)
      df_iter$shat_T12_log = ifelse(df_iter$shat_T12_log >= 0.9999, 0.9999, df_iter$shat_T12_log)
      df_iter$log_T12 = log(-log(1-df_iter$shat_T12_log))
      df_iter$var_log_T12 = df_iter$var_T12/((log(1-df_iter$shat_T12_log)*(1-df_iter$shat_T12_log))^2)
      
      df_iter$shat_T1_log = ifelse(df_iter$shat_T1 <= 0.0001, 0.0001, df_iter$shat_T1)
      df_iter$shat_T1_log = ifelse(df_iter$shat_T1_log >= 0.9999, 0.9999, df_iter$shat_T1_log)
      df_iter$log_T1 = log(-log(1-df_iter$shat_T1_log))
      df_iter$var_log_T1 = df_iter$var_T1/((log(1-df_iter$shat_T1_log)*(1-df_iter$shat_T1_log))^2)

      df_iter$shat_T2_log = ifelse(df_iter$shat_T2 <= 0.0001, 0.0001, df_iter$shat_T2)
      df_iter$shat_T2_log = ifelse(df_iter$shat_T2_log >= 0.9999, 0.9999, df_iter$shat_T2_log)
      df_iter$log_T2 = log(-log(1-df_iter$shat_T2_log))
      df_iter$var_log_T2 = df_iter$var_T2/((log(1-df_iter$shat_T2_log)*(1-df_iter$shat_T2_log))^2)
      
      # Agrestti coulli CI
      df_iter$var_T12_ac = df_iter$shat_T12_ac*(1-df_iter$shat_T12_ac)/(n+k^2)
      
      df_mi <- rbind(df_mi,df_iter)
    }
  }
  
  ep = quantile(avg.epsilon, na.rm = T)
  
  
  ####################################################################################
  ####################### Multiple Imputation Calculations ###########################
  ####################################################################################
  
  df_mi <- 
    df_mi |>
    group_by(t1, t2, true_prob, true_prob_T1, true_prob_T2) |>
     summarise(
      # MI across the M imputed datasets #
      #### Bivariate ####
        shat_T12_bar = mean(shat_T12),
        u_T12 = mean(var_T12),
        b_T12 = var(shat_T12),
  
        # log(log(1-S(t_1,t_2)))
        log_T12_bar = mean(log_T12),
        u_log_T12 = mean(var_log_T12),
        b_log_T12 = var(log_T12),
        
        #Wilson CI
        wilson_LB_T12 = mi_wilson_phat(shat_T12, n, summaries = F)[1],
        wilson_UB_T12 = mi_wilson_phat(shat_T12, n, summaries = F)[2],
      
        #Wald CI
        wald_LB_T12 = mi_wald_phat(shat_T12, n, summaries = F)[1],
        wald_UB_T12 = mi_wald_phat(shat_T12, n, summaries = F)[2],
      
        #AC CI
        shat_T12_bar_ac = mean(shat_T12_ac),
        u_T12_ac = mean(var_T12_ac),
        b_T12_ac = var(shat_T12_ac),
      
      #### First Marginal ####
        shat_T1_bar = mean(shat_T1),
        u_T1 = mean(var_T1),
        b_T1 = var(shat_T1),
        # log(log(1-S(t_1,0)))
        log_T1_bar = mean(log_T1),
        u_log_T1 = mean(var_log_T1),
        b_log_T1 = var(log_T1),

      #### Second Marginal ####
        shat_T2_bar = mean(shat_T2),
        u_T2 = mean(var_T2),
        b_T2 = var(shat_T2),
        # log(log(1-S(0,t_2)))
        log_T2_bar = mean(log_T2),
        u_log_T2 = mean(var_log_T2),
        b_log_T2 = var(log_T2)) |>
    mutate(
      #Make sure there is no problems when dividing by zero
      #If b_T12 (or any other) is zero, i.e., always same phat for all M datsets
      #change it to 0.00001
      b_T12 = ifelse(b_T12 <= 0, 0.00001, b_T12),
      b_log_T12 = ifelse(b_log_T12 <= 0, 0.00001, b_log_T12),
      
      b_T12_ac = ifelse(b_T12_ac <= 0, 0.00001, b_T12_ac),
      
      b_T1 = ifelse(b_T1 <= 0, 0.00001, b_T1),
      b_log_T1 = ifelse(b_log_T1 <= 0, 0.00001, b_log_T1),
      
      b_T2 = ifelse(b_T2 <= 0, 0.00001, b_T2),
      b_log_T2 = ifelse(b_log_T2 <= 0, 0.00001, b_log_T2)) |>
    rename(
      # Rename so it is easier
      shat_T12 = shat_T12_bar, log_T12 = log_T12_bar, shat_T12_ac = shat_T12_bar_ac,
      shat_T1 = shat_T1_bar, log_T1 = log_T1_bar,
      shat_T2 = shat_T2_bar, log_T2 = log_T2_bar) |>
    mutate(
      # Auxiliary variables for CI estimation
      # Bivariate
      se_T12 = sqrt(u_T12 + b_T12*(1 + 1/M)),
      v_T12 = (M-1)*((1+(M/(M+1))*u_T12/b_T12)^2),
      se_log_T12 = sqrt(u_log_T12 + b_log_T12*(1 + 1/M)),
      v_log_T12 = (M-1)*((1+(M/(M+1))*u_log_T12/b_log_T12)^2),
      se_T12_ac = sqrt(u_T12_ac + b_T12_ac*(1 + 1/M)),
      v_T12_ac = (M-1)*((1+(M/(M+1))*u_T12_ac/b_T12_ac)^2),

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
      
      #### Bivariate ####
      bias_T12 = true_prob - shat_T12,
      LB_T12 = shat_T12 - qt(0.975,v_T12)*se_T12,
      UB_T12 = shat_T12 + qt(0.975,v_T12)*se_T12,
      coverage_T12 = (true_prob >= LB_T12)*(true_prob <= UB_T12),
      length_T12 = UB_T12 - LB_T12,
      
      log_LB_T12 = 1 - exp(-exp(log_T12 - qt(0.975,v_log_T12)*se_log_T12)),
      log_UB_T12 = 1 - exp(-exp(log_T12 + qt(0.975,v_log_T12)*se_log_T12)),
      coverage_log_T12 = (true_prob >= log_LB_T12)*(true_prob <= log_UB_T12),
      length_log_T12 = log_UB_T12 - log_LB_T12,
      
      coverage_wilson_T12 = (true_prob >= wilson_LB_T12)*(true_prob <= wilson_UB_T12),
      length_wilson_T12 = wilson_UB_T12 - wilson_LB_T12,
      
      coverage_wald_T12 = (true_prob >= wald_LB_T12)*(true_prob <= wald_UB_T12),
      length_wald_T12 = wald_UB_T12 - wald_LB_T12,
      
      LB_T12_ac = shat_T12_ac - qt(0.975,v_T12_ac)*se_T12_ac,
      UB_T12_ac = shat_T12_ac + qt(0.975,v_T12_ac)*se_T12_ac,
      coverage_T12_ac = (true_prob >= LB_T12_ac)*(true_prob <= UB_T12_ac),
      length_T12_ac = UB_T12_ac - LB_T12_ac,
      
      #### Marginal T1 ####
      bias_T1 = true_prob_T1 - shat_T1,
      LB_T1 = shat_T1 - qt(0.975,v_T1)*se_T1,
      UB_T1 = shat_T1 + qt(0.975,v_T1)*se_T1,
      coverage_T1 = (true_prob_T1 >= LB_T1)*(true_prob_T1 <= UB_T1),
      length_T1 = UB_T1 - LB_T1,
      log_LB_T1 = 1 - exp(-exp(log_T1 - qt(0.975,v_log_T1)*se_log_T1)),
      log_UB_T1 = 1 - exp(-exp(log_T1 + qt(0.975,v_log_T1)*se_log_T1)),
      coverage_log_T1 = (true_prob_T1 >= log_LB_T1)*(true_prob_T1 <= log_UB_T1),
      length_log_T1 = log_UB_T1-log_LB_T1,
      
      
      #### Marginal T2 ####
      bias_T2 = true_prob_T2 - shat_T2,
      LB_T2 = shat_T2 - qt(0.975,v_T2)*se_T2,
      UB_T2 = shat_T2 + qt(0.975,v_T2)*se_T2,
      coverage_T2 = (true_prob_T2 >= LB_T2)*(true_prob_T2 <= UB_T2),
      length_T2 = UB_T2 - LB_T2,
      log_LB_T2 = 1 - exp(-exp(log_T2 - qt(0.975,v_log_T2)*se_log_T2)),
      log_UB_T2 = 1 - exp(-exp(log_T2 + qt(0.975,v_log_T2)*se_log_T2)),
      coverage_log_T2 = (true_prob_T2 >= log_LB_T2)*(true_prob_T2 <= log_UB_T2),
      length_log_T2 = log_UB_T2 - log_LB_T2,
      
      #Epsilon
      epsilon.mean = mean(avg.epsilon, na.rm = T),
      epsilon.min = ep[1],
      epsilon.q1 = ep[2],
      epsilon.q2 = ep[3],
      epsilon.q3 = ep[4],
      epsilon.max = ep[5]) |>
    mutate(
      # Simulation configuration
      correlation = rho,
      censoring = censoring,
      n = n, 
      sim = sim,
      method = "IT-MI")
  

  df_mi <-
    df_mi |>
    dplyr::select(correlation:method,
      t1:shat_T12, bias_T12, se_T12, coverage_T12, length_T12, coverage_log_T12, length_log_T12, 
      coverage_wilson_T12, length_wilson_T12, coverage_wald_T12, length_wald_T12, coverage_T12_ac, length_T12_ac, 
      shat_T1, bias_T1, se_T1, coverage_T1, length_T1, coverage_log_T1, length_log_T1,
      shat_T2, bias_T2, se_T2, coverage_T2, length_T2, coverage_log_T2, length_log_T2, 
      epsilon.mean:epsilon.max)

  
  ####################################################################################
  #################################### SAVE ##########################################
  ####################################################################################
  
  name_all = paste0("rho_",rho,"_cens_",censoring,"_n_",n,"_sim_",sim,".RData")
  name = paste0("MI/results-all/",name_all)
  save(df_mi, file = name)
  
  


  





