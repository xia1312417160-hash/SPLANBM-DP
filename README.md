# SPLANBM-DP

Code for "Semi-parametric Bayesian inference to shape-restricted partially linear additive mixed-effects negative binomial models for longitudinal data" by Miaojie Xia, Ya Liu, and Jiang Du.

## Files

- `SPLANBM_DP.R`: model and Gibbs sampler.
- `simulation_example.R`: reproducible simulation example.

## Run

Install the R packages `BayesLogit` and `mvtnorm`, then run the following command from the RStudio Console or with the Source button:

```r
source("simulation_example.R")
```

The posterior means of $\beta_1$, $\beta_2$, and $r$ are printed, and the estimated functions are displayed in the RStudio Plots pane. Running the script with `Rscript` in the Terminal or as a Background Job will not display the plot in that pane.

The shape codes are `1` increasing, `2` decreasing, `3` convex, and `4` concave.

The example uses the paper's Case A setting with $(n,n_i)=(50,10)$, $K=3$, $G=200$, 4,000 MCMC iterations, and 2,000 burn-in iterations.
