# Multiple-Imputation-of-Censored-Bivariate-Event-times
Multiple Imputation of Censored  Bivariate Event-times via Inverse Transform and Nonparametric Gibbs Sampling

The *data* folder has all the datasets created for the simulation analysis under two different censoring distributions, uniform and exponential. The file names indicate the sampling configuration used (regarding correlation and the marginal censoring rate). 

The *code* folder contains three R files: 00_helpers.R, IT-MI.R, and benchmark.R. The former only has functions later called in the other R files. The IT-MI.R file recreates one replication result using our IT-MI method under a specific sampling distribution. The user can specify all parameters needed, such as correlation $\rho \in \{0.5,0.8\}$, marginal censoring rate, sample size, and censoring distribution. Similarly, the benchmark.R file recreates one replication result using the Dabrowska bivariate survival estimator. Simulation results presented in the paper show an average of over 1000 replications. IT-MI and benchmark approaches take a simulation number from 1 to 1000 and create the replication result for that seed.
