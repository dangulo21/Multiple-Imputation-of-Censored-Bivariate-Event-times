# Multiple-Imputation-of-Censored-Bivariate-Event-times
Multiple Imputation of Censored  Bivariate Event-times via Inverse Transform and Nonparametric Gibbs Sampling

The *data* folder has all the datasets created for the simulation analysis under two different censoring distributions, uniform and exponential. The file names indicate the sampling configuration used (regarding correlation and the marginal censoring rate). 

The *code* folder contains four R files: `00_helpers.R`, `IT-MI.R`, `benchmark.R`, and `final-results.R`. The former only has functions later called in the other R files. The `IT-MI.R` file recreates one replication result using our IT-MI method under a specific sampling distribution. The user can specify all parameters needed, such as correlation: $\rho \in ${0.5,0.8}, marginal censoring rate $\in ${0,40,50}, sample size: $n \in ${100,250,500}, and censoring distribution: uniform or exponential. Similarly, the `benchmark.R` file recreates one replication result using the Dabrowska bivariate survival estimator. Simulation results presented in the paper show an average of over 1000 replications. IT-MI and benchmark approaches take a simulation number from 1 to 1000 and create the replication result for that seed. Lastly, `final-results.R` file combines the results of the 1000 simulations for each sampling scenario.


Finally, the *application* folder contains two files: `IT-MI-app.R` and `benchmark-app.R`. Each file computes the application results, our methodology, and Dabrowska's comparison, respectively. 
