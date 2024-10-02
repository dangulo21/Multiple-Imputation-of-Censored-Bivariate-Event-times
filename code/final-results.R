library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(latex2exp)


#########################################################################
########################### Scenario ####################################
#########################################################################

cens_type = "exponential"
rho = 0.5
censoring = 50
n = 250

#########################################################################
############################## IT-MI ####################################
#########################################################################

#Change the path to where the 1000 replication results using the IT-MI method
#for the above scenario are:
path = "results-all/..."
files = list.files(path)
  
for(file in files){
  load(paste0(path,"/",file))
  
  #Calculate some useful variables
  results_combine <-
    results_combine |>
    mutate(
      #Bias
      bias.T12 = true_prob - prob_T12,
      bias.T1 = true_prob_T1 - prob_T1,
      bias.T2 = true_prob_T2 - prob_T2,
      
      #BIVARIATE
      LB.T12 = prob_T12 - qt(0.975,v_T12)*se_prob_T12,
      UB.T12 = prob_T12 + qt(0.975,v_T12)*se_prob_T12,
      coverage.T12 = (true_prob >= LB.T12)*(true_prob <= UB.T12),
      length.T12 = UB.T12-LB.T12,
      log.log.LB = prob_T12^(exp(+qt(0.975,v_loglog_T12)*se_loglog_prob_T12)),
      log.log.UB = prob_T12^(exp(-qt(0.975,v_loglog_T12)*se_loglog_prob_T12)),
      coverage.loglog.T12 = (true_prob >= log.log.LB)*(true_prob <= log.log.UB),
      length.loglog.T12 = log.log.UB-log.log.LB,
      
      #MARGINAL T1
      LB.T1 = prob_T1 - qt(0.975,v_T1)*se_prob_T1,
      UB.T1 = prob_T1 + qt(0.975,v_T1)*se_prob_T1,
      coverage.T1 = (true_prob_T1 >= LB.T1)*(true_prob_T1 <= UB.T1),
      length.T1 = UB.T1-LB.T1,
      log.log.LB.T1 = prob_T1^(exp(+qt(0.975,v_loglog_T1)*se_loglog_prob_T1)),
      log.log.UB.T1 = prob_T1^(exp(-qt(0.975,v_loglog_T1)*se_loglog_prob_T1)),
      coverage.loglog.T1 = (true_prob_T1 >= log.log.LB.T1)*(true_prob_T1 <= log.log.UB.T1),
      length.loglog.T1 = log.log.UB.T1-log.log.LB.T1,
      
      
      #MARGINAL T2
      LB.T2 = prob_T2 - qt(0.975,v_T2)*se_prob_T2,
      UB.T2 = prob_T2 + qt(0.975,v_T2)*se_prob_T2,
      coverage.T2 = (true_prob_T2 >= LB.T2)*(true_prob_T2 <= UB.T2),
      length.T2 = UB.T2-LB.T2,
      log.log.LB.T2 = prob_T2^(exp(+qt(0.975,v_loglog_T2)*se_loglog_prob_T2)),
      log.log.UB.T2 = prob_T2^(exp(-qt(0.975,v_loglog_T2)*se_loglog_prob_T2)),
      coverage.loglog.T2 = (true_prob_T2 >= log.log.LB.T2)*(true_prob_T2 <= log.log.UB.T2),
      length.loglog.T2 = log.log.UB.T2-log.log.LB.T2,
    )
  
  
  #Summarize results across simulations
  results_combine <-
    results_combine |>
    group_by(
      correlation, censoring, n, t1, t2,
      true_prob,true_prob_T1,true_prob_T2) |>
    summarize(
      
      #Both Endpoints
      avg.prob.T12 = mean(prob_T12, na.rm = T),
      avg.bias.T12 = mean(bias.T12, na.rm = T),
      avg.se.T12 = mean(se_prob_T12, na.rm = T),
      sd.T12 = sd(prob_T12, na.rm = T),
      cov.prob.T12 = mean(coverage.T12, na.rm = T),
      avg.length.T12 = mean(length.T12, na.rm = T),
      cov.prob.loglog.T12 = mean(coverage.loglog.T12, na.rm = T),
      avg.length.loglog.T12 = mean(length.loglog.T12, na.rm = T),
      
      #Endpoint T1
      avg.prob.T1 = mean(prob_T1, na.rm = T),
      avg.bias.T1 = mean(bias.T1, na.rm = T),
      avg.se.T1 = mean(se_prob_T1, na.rm = T),
      sd.T1 = sd(prob_T1, na.rm = T),
      cov.prob.T1 = mean(coverage.T1, na.rm = T),
      avg.length.T1 = mean(length.T1, na.rm = T),
      cov.prob.loglog.T1 = mean(coverage.loglog.T1, na.rm = T),
      avg.length.loglog.T1 = mean(length.loglog.T1, na.rm = T),
      
      #Endpoint T2
      avg.prob.T2 = mean(prob_T2, na.rm = T),
      avg.bias.T2 = mean(bias.T2, na.rm = T),
      avg.se.T2 = mean(se_prob_T2, na.rm = T),
      sd.T2 = sd(prob_T2, na.rm = T),
      cov.prob.T2 = mean(coverage.T2, na.rm = T),
      avg.length.T2 = mean(length.T2, na.rm = T),
      cov.prob.loglog.T2 = mean(coverage.loglog.T2, na.rm = T),
      avg.length.loglog.T2 = mean(length.loglog.T2, na.rm = T))
  
  results = rbind(results,results_combine)
  remove(results_combine)
}
  
results$method = "IT-MI"


#########################################################################
########################### Benchmarking ################################
#########################################################################

#Change the path to where the 1000 replication results using the Dabrowska
#bivariate survival estimator for the above scenario are:
path = "results-all/..."
files = list.files(path)
  
for(file in files){
  load(paste0(path,"/",file))
  
  #Calculate some useful variables
  results_combine <-
    results_combine |>
    mutate(
      #Bias
      bias.T12 = true_prob - prob.T12,
      bias.T1 = true_prob_T1 - prob.T1,
      bias.T2 = true_prob_T2 - prob_T2,
      
      #BIVARIATE
      coverage.T12 = (true_prob >= LB.T12)*(true_prob <= UB.T12),
      length.T12 = UB.T12-LB.T12,
      log.log.LB = prob.T12^exp(-qnorm(0.975)*(se.T12/(log(prob.T12)*prob.T12))),
      log.log.UB = prob.T12^exp(+qnorm(0.975)*(se.T12/(log(prob.T12)*prob.T12))),
      coverage.loglog.T12 = (true_prob >= log.log.LB)*(true_prob <= log.log.UB),
      length.loglog.T12 = log.log.UB-log.log.LB,
      
      #MARGINAL T1
      coverage.T1 = (true_prob_T1 >= LB.T1)*(true_prob_T1 <= UB.T1),
      length.T1 = UB.T1-LB.T1,
      log.log.LB.T1 = prob.T1^exp(-qnorm(0.975)*(se.T1/(log(prob.T1)*prob.T1))),
      log.log.UB.T1 = prob.T1^exp(+qnorm(0.975)*(se.T1/(log(prob.T1)*prob.T1))),
      coverage.loglog.T1 = (true_prob_T1 >= log.log.LB.T1)*(true_prob_T1 <= log.log.UB.T1),
      length.loglog.T1 = log.log.UB.T1-log.log.LB.T1,

      #MARGINAL T2
      coverage.T2 = (true_prob_T2 >= LB.T2)*(true_prob_T2 <= UB.T2),
      length.T2 = UB.T2-LB.T2,
      log.log.LB.T2 = prob_T2^exp(-qnorm(0.975)*(se.T2/(log(prob_T2)*prob_T2))),
      log.log.UB.T2 = prob_T2^exp(+qnorm(0.975)*(se.T2/(log(prob_T2)*prob_T2))),
      coverage.loglog.T2 = (true_prob_T2 >= log.log.LB.T2)*(true_prob_T2 <= log.log.UB.T2),
      length.loglog.T2 = log.log.UB.T2 - log.log.LB.T2)
  
  
  #Summarize results across simulations
  results_combine <-
    results_combine |>
    group_by(
      correlation, censoring, n, t1, t2,
      true_prob,true_prob_T1,true_prob_T2) |>
    summarize(
      #Both Endpoints
      avg.prob.T12 = mean(prob.T12, na.rm = T),
      avg.bias.T12 = mean(bias.T12, na.rm = T),
      avg.se.T12 = mean(se.T12, na.rm = T),
      sd.T12 = sd(prob.T12, na.rm = T),
      cov.prob.T12 = mean(coverage.T12, na.rm = T),
      avg.length.T12 = mean(length.T12, na.rm = T),
      cov.prob.loglog.T12 = mean(coverage.loglog.T12, na.rm = T),
      avg.length.loglog.T12 = mean(length.loglog.T12, na.rm = T),

      #Endpoint T1
      avg.prob.T1 = mean(prob.T1, na.rm = T),
      avg.bias.T1 = mean(bias.T1, na.rm = T),
      sd.T1 = sd(prob.T1, na.rm = T),
      cov.prob.T1 = mean(coverage.T1, na.rm = T),
      avg.length.T1 = mean(length.T1, na.rm = T),
      avg.se.T1 = mean(se.T1, na.rm = T),
      cov.prob.loglog.T1 = mean(coverage.loglog.T1, na.rm = T),
      avg.length.loglog.T1 = mean(length.loglog.T1, na.rm = T),

      #Endpoint T2
      avg.prob.T2 = mean(prob_T2, na.rm = T),
      avg.bias.T2 = mean(bias.T2, na.rm = T),
      sd.T2 = sd(prob_T2, na.rm = T),
      cov.prob.T2 = mean(coverage.T2, na.rm = T),
      avg.length.T2 = mean(length.T2, na.rm = T),
      avg.se.T2 = mean(se.T2, na.rm = T),
      cov.prob.loglog.T2 = mean(coverage.loglog.T2, na.rm = T),
      avg.length.loglog.T2 = mean(length.loglog.T2, na.rm = T))
  
  results2 = rbind(results2,results_combine)
  remove(results_combine)
} 

results2$method = "Dabrowska"


