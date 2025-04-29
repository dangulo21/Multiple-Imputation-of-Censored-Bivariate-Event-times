restricted_marginal <- function(id,
                                var.time,var.delta,
                                covar.time,covar.delta,
                                df_imp){
   
  # Euclidean parameter distance, epsilon
  epsilon = 0.2
  max.epsilon = tau  #the maximum difference is tau (someone with X=tau and someone with X=0, theoretically)
  
  ### Set the risk set ###
  
  # To ensure we have a proper imputation, in the risk set the largest observed time has to be and event
  # If that does not happen, increase delta until it does. 
  # By definition, since our RVs are min(T,tau) - we checked all simulation seeds - in every simulated dataset the last observed time for both endpoints is an event.
  # How much can epsilon increase until the last observed time is an event? tau 
  
  # Initial risk set: i) individuals at risk for the event being imputed, ii) close euclidean distance for the other endpoint

  # The minimal size of the risk set cn be set by: i. minimmum number of events, and ii. minimmum number of individuals
  # The first one i. minimmum number of events, is set up with:
  cut.off = min(sum(df_imp[which(df_imp[,var.time] > df_imp[id,var.time]),var.delta]) , 10)

  # The second one ii. minimmum number of individuals, is set up with:
  #cut.off = min(length(which(df_imp[,var.time] > df_imp[id,var.time])), 10)
  
  risk.set <- which(df_imp[,var.time] > df_imp[id,var.time] & 
                      abs(df_imp[,covar.time] - df_imp[id,covar.time]) <= epsilon)
  
  # The initial risk set could have zero participants (no subjects pass the euclidean distance test)
  # If that is the case, increase epsilon until we have at least ONE subject in the risk set. 
  # Note: it is necessary to check to have at least one subject in the risk set before moving on because
  # for proper imputation to occur we need to check that the largest observed time (for the endpoint being imputed) 
  # is an event. Hence, need to have at least one participant in the risk set to check for second condition.
  
  while(length(risk.set) == 0){
    epsilon = epsilon + 0.001
    risk.set <- which(df_imp[,var.time] > df_imp[id,var.time] & 
                        abs(df_imp[,covar.time] - df_imp[id,covar.time]) <= epsilon)}
  
  # Once we have at least one participant in the risk set, we recovered their information. 
  df_risk_set <- df_imp[risk.set,]
  
  # Iterate until we find a risk set that: i. has at least 5/cut off No. of participants, and 
  # ii. the largest observed time (for the endpoint being imputed) is an event
  # The loop is until max.epsilon because the maximum.epsilon = tau, will guaranteed that the second condition is met.
  
  while(epsilon <= max.epsilon){
    
    n.events = sum(df_risk_set[,var.delta])
    
    if(n.events >= cut.off & df_risk_set[which.max(df_risk_set[,var.time]),var.delta]){break}
    
    # If conditions are not met, increase epsilon. 
    epsilon = epsilon + 0.001
    risk.set <- which(df_imp[,var.time] > df_imp[id,var.time] & 
                        abs(df_imp[,covar.time] - df_imp[id,covar.time]) <= epsilon)
    
    df_risk_set <- df_imp[risk.set,]
  }
  
  # Once we have the appropriate risk set, we estimate the KM among the individuals in the risk set.
  km <- survfit(Surv(eval(parse(text = var.time)), eval(parse(text = var.delta))) ~ 1, data = df_risk_set)
  km <- summary(km)
  
  # Draw from an uniform distribution
  u <- runif(1)
  
  # Inverse transform
  impute.index = which.max(km$surv <= u)
  return(list(impute = km$time[impute.index], epsilon = epsilon)) 
  
}






