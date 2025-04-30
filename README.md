# Multiple Imputation of Censored  Bivariate Event-times via Inverse Transform and Nonparametric Gibbs Sampling

## Descrption 

The `Multiple Imputation of Censored Bivariate Event-times` houses all data and source files for the Multiple Imputation of Censored  Bivariate Event-times via Inverse Transform and Nonparametric Gibbs Sampling article by Daniela Angulo and Susan Murray. 

## Data sources
 
The *data* folder contains all the datasets created for the simulation analysis. The folder names indicate time-to-event distribution and censoring-times distribution, respectively. For example, `exponential-uniform` refers to the data configuration of bivariate exponential event-times with uniform censoring. Within each folder, file names indicate the sampling configuration used (regarding correlation and censoring rate). 

## Source files

The *code* folder contains five `R` files: `00_helpers.R`, `IT-MI.R`, `B-IT-MI.R`, `Darbowska.R`. Please see a short description of each of these files below:

1.  `00_helpers.R`: This file contains the function `restriced_marginal(...)` that is recursevely call in the IT-MI algorithm. There, the user can define user-specific parameters such as $\epsilon$, $d_{\mathcal{R}}$ or $n_{\mathcal{R}}$.
2.  `IT-MI.R`: This source files implements the inverse transform multiple imputation (IT-MI) method. The user can define parameters for the simulation configuration such as: distribution between event-times and censoring-times, correlation between event-times, \rho \in ${0.5,0.8}, and sample size: $n \in ${100,250,500}.
3. `B-IT-MI.R`: This source files implements the bootstrap verion of IT-MI, called the B-IT-MI method. The user can define parameters for the simulation configuration such as: distribution between event-times and censoring-times, correlation between event-times, \rho \in ${0.5,0.8}, and sample size: $n \in ${100,250,500}.
4.  `Darbowska.R`: This source files implements the Dabrowska estimator (using functions from the `mhazard` package).


## HIV/AIDS application 

Finally, the *application* folder contains two files: `IT-MI-app.R` and `benchmark-app.R`. Each file computes the application results, our methodology, and Dabrowska's comparison, respectively. 
