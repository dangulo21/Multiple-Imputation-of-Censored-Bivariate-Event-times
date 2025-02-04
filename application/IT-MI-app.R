######################################################################################
######################################################################################
############################ APPLICATION - HPV/AIDS ##################################
######################################################################################
######################################################################################

library(dplyr)
library(tidyr)
library(survival)
library(MASS)
library(Rfast)
library(mgcv)
library(survSpearman)

  source("code/00_helpers.R")
  
  ######################################################################################
  ################################## READ DATA #########################################
  ######################################################################################
  
  df <- survSpearman::CCASAnetData
  
  tau = 15
 
  df <- 
    df |>
    mutate(
      timeToRegimenChange = pmin(timeToRegimenChange, tau),
      timeToViralFailure = pmin(timeToViralFailure, tau),
      regimenChange = ifelse(timeToRegimenChange >= tau, 1, regimenChange),
      viralFailure = ifelse(timeToViralFailure >= tau, 1, viralFailure))
 
  data_sim <- 
    df |>
    dplyr::select(-site)
  
  
  colnames(data_sim) <- c("x1","delta1","x2","delta2")
  
  L = 10
  M = 10
  n = nrow(data_sim)
  
  ## df_comp
  df_comp <- 
    data.frame(
      expand_grid(t1 = seq(0.0,tau-0.5,0.5), t2 = seq(0.0,tau-0.5,0.5)))
  
  set.seed(19980305)
  epsilon = 0.75
  
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
        delta1 == F ~ NA,
      ),
      impute.T2 = case_when(
        delta2 == T ~ x2,
        delta2 == F ~ NA,
      )
    )
  
  ######################################################################################
  ############################ Initialization Stage ####################################
  ######################################################################################
  
  #Dimension to impute first
  time.first <- ifelse(sum(df_imp$delta1) > sum(df_imp$delta2), "x1", "x2")
  delta.first <- ifelse(sum(df_imp$delta1) > sum(df_imp$delta2), "delta1", "delta2")
  impute.first <- ifelse(sum(df_imp$delta1) > sum(df_imp$delta2), "impute.T1", "impute.T2") 
  
  
  #Marginal imputation
  for(i in which(!df_imp[,delta.first])){
    
    risk.set <- which(df_imp[,time.first] > df_imp[i,time.first])
    
    df_risk_set <- df_imp[risk.set,c(time.first,delta.first)]
    
    if(length(risk.set) == 0){
      next
    }
    
    km <- survfit(Surv(eval(parse(text = time.first)), eval(parse(text = delta.first))) ~ 1, data = df_risk_set)
    km <- summary(km)
    
    u <- runif(1)
    if(u <= min(km$surv)){
      df_imp[i,impute.first] = tau
    }else{
      impute.index = which.max(km$surv <= u)
      df_imp[i,impute.first] = km$time[impute.index]
    }
    
  }
  
  df_imp[which(is.na(get.var(impute.first,df_imp))),impute.first] <- df_imp[which(is.na(get.var(impute.first,df_imp))),time.first] 
  
  
  ##########################################################################################
  ############################ Iterative Imputation Stage ##################################
  ##########################################################################################
  
  iter = 0
  while(iter < (L+M)*2){
    iter = iter + 1
    
    #Select dimension to impute
    time.second <- ifelse(time.first == "x1","x2","x1")
    delta.second <- ifelse(delta.first == "delta1","delta2","delta1")
    impute.second <- ifelse(impute.first == "impute.T1","impute.T2","impute.T1")
    
    #clean dimension to impute
    df_imp[which(!df_imp[,delta.second]),impute.second] <- NA
    
    #Restricted marginal
    for(i in which(!df_imp[,delta.second])){
      df_imp[i,impute.second] = restricted_marginal(i,time.second,delta.second,impute.first,delta.first,epsilon)
    }
    
    df_imp[which(is.na(get.var(impute.second,df_imp))),impute.second] <- df_imp[which(is.na(get.var(impute.second,df_imp))),time.second] 
    
    #Change labels
    time.first = time.second
    delta.first = delta.second
    impute.first = impute.second
    
    ####################################################################################
    ################### Analysis of Augmented Dataset Stage ############################
    ####################################################################################
    
    if(iter > L*2 & iter %% 2 != 0){
       
      df_res_imp <- df_comp
      df_res_imp$imp = iter
      
      for(i in 1:nrow(df_res_imp)){
        t1 = df_comp[i,]$t1
        t2 = df_comp[i,]$t2
        
        df_res_imp$est_prob[i] = nrow(df_imp[df_imp$impute.T1 > t1 & df_imp$impute.T2 > t2,])/n
        df_res_imp$est_prob_T1[i] = nrow(df_imp[df_imp$impute.T1 > t1,])/n
        df_res_imp$est_prob_T2[i] = nrow(df_imp[df_imp$impute.T2 > t2,])/n
      }
      
      df_res_imp$var_est_prob = df_res_imp$est_prob*(1-df_res_imp$est_prob)/n
      df_res_imp$var_est_prob_T1 = df_res_imp$est_prob_T1*(1-df_res_imp$est_prob_T1)/n
      df_res_imp$var_est_prob_T2 = df_res_imp$est_prob_T2*(1-df_res_imp$est_prob_T2)/n
      
      #Complementary log-log transform
      df_res_imp$est_log_log = log(-log(df_res_imp$est_prob)) 
      df_res_imp$var_log_log = df_res_imp$var_est_prob/((log(df_res_imp$est_prob)*df_res_imp$est_prob)^2)
      
      df_res_imp$est_log_log_T1 = log(-log(df_res_imp$est_prob_T1)) 
      df_res_imp$var_log_log_T1 = df_res_imp$var_est_prob_T1/((log(df_res_imp$est_prob_T1)*df_res_imp$est_prob_T1)^2)
      
      df_res_imp$est_log_log_T2 = log(-log(df_res_imp$est_prob_T2)) 
      df_res_imp$var_log_log_T2 = df_res_imp$var_est_prob_T2/((log(df_res_imp$est_prob_T2)*df_res_imp$est_prob_T2)^2)
      
      #Correlation
      df_res_imp$est_corr = cor(df_imp$impute.T1,df_imp$impute.T2)
      
      df_res <- rbind(df_res,df_res_imp)
      remove(df_res_imp)
      
    }
    
  }
  
  
  
  ####################################################################################
  ####################### Multiple Imputation Calculations ###########################
  ####################################################################################
  
  df_res <- 
    df_res |>
    group_by(t1,t2) |>
     summarise(
      #Bivariate
      prob_T12 = mean(est_prob),
      u_T12 = mean(var_est_prob),
      b_T12 = var(est_prob),
      loglog_prob_T12 = mean(est_log_log),
      u_loglog_T12 = mean(var_log_log),
      b_loglog_T12 = var(est_log_log),
      
      #First Marginal
      prob_T1 = mean(est_prob_T1),
      u_T1 = mean(var_est_prob_T1),
      b_T1 = var(est_prob_T1),
      loglog_prob_T1 = mean(est_log_log_T1),
      u_loglog_T1 = mean(var_log_log_T1),
      b_loglog_T1 = var(est_log_log_T1),
      
      #Second Marginal
      prob_T2 = mean(est_prob_T2),
      u_T2 = mean(var_est_prob_T2),
      b_T2 = var(est_prob_T2),
      loglog_prob_T2 = mean(est_log_log_T2),
      u_loglog_T2 = mean(var_log_log_T2),
      b_loglog_T2 = var(est_log_log_T2),
      
      #correlation
      corr = mean(est_corr),
      
      M = n()) |>
    mutate(
      #Bivariate
      se_prob_T12 = sqrt(u_T12+b_T12*(1+1/M)),
      v_T12 = (M-1)*((1+(M/(M+1))*u_T12/b_T12)^2),
      se_loglog_prob_T12 = sqrt(u_loglog_T12+b_loglog_T12*(1+1/M)),
      v_loglog_T12 = (M-1)*((1+(M/(M+1))*u_loglog_T12/b_loglog_T12)^2),
      
      #First Marginal
      se_prob_T1 = sqrt(u_T1+b_T1*(1+1/M)),
      v_T1 = (M-1)*((1+(M/(M+1))*u_T1/b_T1)^2),
      se_loglog_prob_T1 = sqrt(u_loglog_T1+b_loglog_T1*(1+1/M)),
      v_loglog_T1 = (M-1)*((1+(M/(M+1))*u_loglog_T1/b_loglog_T1)^2),
      
      #Second Marginal
      se_prob_T2 = sqrt(u_T2+b_T2*(1+1/M)),
      v_T2 = (M-1)*((1+(M/(M+1))*u_T2/b_T2)^2),
      se_loglog_prob_T2 = sqrt(u_loglog_T2+b_loglog_T2*(1+1/M)),
      v_loglog_T2 = (M-1)*((1+(M/(M+1))*u_loglog_T2/b_loglog_T2)^2),
      
      n = n)
  
  ####################################################################################
  #################################### SAVE ##########################################
  ####################################################################################
   
  name = paste0("results-final/application/IT-MI-app.RData")
  save(df_res, file = name)
  












