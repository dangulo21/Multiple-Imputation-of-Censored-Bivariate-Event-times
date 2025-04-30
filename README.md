# Multiple Imputation of Censored  Bivariate Event-times via Inverse Transform and Nonparametric Gibbs Sampling

## Descrption 

The `Multiple Imputation of Censored Bivariate Event-times` houses all data and source files for the Multiple Imputation of Censored  Bivariate Event-times via Inverse Transform and Nonparametric Gibbs Sampling article by Daniela Angulo and Susan Murray. 

## Data sources
 
The *data* folder contains all the datasets created for the simulation analysis. The folder names indicate time-to-event distribution and censoring-times distribution, respectively. For example, `exponential-uniform` refers to the data configuration of bivariate exponential event-times with uniform censoring. Within each folder, file names indicate the sampling configuration used (regarding correlation and censoring rate). 

## Source files

The *code* folder contains five `R` files: `00_helpers.R`, `IT-MI.R`, `B-IT-MI.R`, `Darbowska.R`, and `final-results.R`. Please see a short description of each of these files below:

1.  `00_helpers.R`: This file contains the function `restriced_marginal(...)` that is recursevely call in the IT-MI algorithm. There, the user can define user-specific parameters such as $\epsilon$, $d_{\mathcal{R}}$ or $n_{\mathcal{R}}$. 

The former only has functions later called in the other R files. The `IT-MI.R` file recreates one replication result using our IT-MI method under a specific sampling distribution. The user can specify all parameters needed, such as correlation: $\rho \in ${0.5,0.8}, marginal censoring rate $\in ${0,40,50}, sample size: $n \in ${100,250,500}, and censoring distribution: uniform or exponential. Similarly, the `benchmark.R` file recreates one replication result using the Dabrowska bivariate survival estimator. Simulation results presented in the paper show an average of over 1000 replications. IT-MI and benchmark approaches take a simulation number from 1 to 1000 and create the replication result for that seed. Lastly, `final-results.R` file combines the results of the 1000 simulations for each sampling scenario.

## HIV/AIDS application 

Finally, the *application* folder contains two files: `IT-MI-app.R` and `benchmark-app.R`. Each file computes the application results, our methodology, and Dabrowska's comparison, respectively. 
