restricted_marginal <- function(id,
                                var.time, var.delta,
                                covar.time, covar.delta){
    
  #Set starting epsilon
  epsilon = 1.0
  max.epsilon = 5
  cut.off = 5
      
  ### Set the risk set ###
  risk.set <- which(df_imp[,var.time] > df_imp[id,var.time] & 
                abs(df_imp[,covar.time] - df_imp[id,covar.time]) <= epsilon)
  df_risk_set <- df_imp[risk.set,]
      
  while(epsilon < max.epsilon){
      if(length(risk.set) >= cut.off){break}
      epsilon = epsilon + 0.001
      risk.set <- which(df_imp[,var.time] > df_imp[id,var.time] & 
                       abs(df_imp[,covar.time] - df_imp[id,covar.time]) <= epsilon)
            
      df_risk_set <- df_imp[risk.set,]
  }
    
  if(length(risk.set) < cut.off){return(NA)}
      
  km <- survfit(Surv(eval(parse(text = var.time)), eval(parse(text = var.delta))) ~ 1, data = df_risk_set)
  km <- summary(km)
  
  u <- runif(1)
  if(u <= min(km$surv)){
    return(tau)
  }else{
    impute.index = which.max(km$surv <= u)
    return(km$time[impute.index]) 
  }
  

}

















