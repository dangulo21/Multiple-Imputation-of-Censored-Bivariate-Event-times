# Multiple Imputation of Censored  Bivariate Event-times via Inverse Transform and Nonparametric Gibbs Sampling

## Descrption 

This repository contains the data and source code associated with the article:  **_Multiple Imputation of Censored Bivariate Event-Times via Inverse Transform and Nonparametric Gibbs Sampling_**  by *Daniela Angulo* and *Susan Murray*.

## Data sources
 
The *data* folder contains all the datasets created for the simulation analysis. Folder names indicate the time-to-event distribution and censoring-time distribution used. For example, `exponential-uniform` refers to datasets with bivariate exponential event-times and uniform censoring. Within each folder, file names reflect the sampling configuration, including correlation and censoring rate.

## Source files

The *code* folder contains five `R` files: `00_helpers.R`, `IT-MI.R`, `B-IT-MI.R`, `Darbowska.R`. Below is a brief description of each:

1.  `00_helpers.R`:  Contains the function `restriced_marginal(...)` which is recursively called in the IT-MI algorithm. Users can specify parameters such as $\epsilon$, $d_{\mathcal{R}}$, and $n_{\mathcal{R}}$.
2.  `IT-MI.R`: Implements the inverse transform multiple imputation (IT-MI) method. Users can define simulation parameters such as the distributions of event-times and censoring-times, event-time correlation $\rho \in {0.5, 0.8}$, and sample sizes $n \in {100, 250, 500}$.
3. `B-IT-MI.R`: Implements the bootstrap version of IT-MI, known as the B-IT-MI method. Users can define simulation parameters such as the distributions of event-times and censoring-times, event-time correlation $\rho \in {0.5, 0.8}$, and sample sizes $n \in {100, 250, 500}$.
4.  `Darbowska.R`: Implements the Dabrowska estimator using functions from the `mhazard` package.


## HIV/AIDS application 

The `application.R` script applies the proposed methods to real-world data from 6,691 HIV/AIDS patients initiating antiretroviral therapy (ART) across five Latin American regions. The data is sourced from the `survSpearman` R package.

We estimate the bivariate survival probabilities from ART initiation to the following events:
1. **Virological failure** ($T_1$)
2. **Major regimen change** ($T_2$)

Both the IT-MI methodology and Dabrowska's estimator are used for analysis.
