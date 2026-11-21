# Monte Carlo Methods

My solutions to the lab sessions of the Monte Carlo methods course at Ensimag in Python. All examples are set in the Black-Scholes model. The course is ongoing, and new sessions are added as they are completed.

## Contents

| Notebook | Topic | What it does |
|---|---|---|
| `01_brownian_hitting_time.ipynb` | Brownian motion simulation | Estimates $\mathbb E[\tau_a \wedge T]$ for the first hitting time of a level $a$, compares it with the closed-form value, and shows the bias from monitoring the path on a discrete grid. |
| `02_variance_reduction.ipynb` | Variance reduction | Prices an option on the discrete maximum, $(\max_j S_{t_j} - K S_T)^+$, with crude Monte Carlo, antithetic variables and a control variate ($S_T$). Then learns an importance sampling drift for a call with a Robbins-Monro stochastic algorithm. |
| `03_importance_sampling.ipynb` | Importance sampling | Prices a discretely monitored down-and-out call. The Girsanov drift $\lambda$ is chosen by Newton's method on an empirical variance criterion, which narrows the 95% confidence interval by a factor of about 3 for the same number of paths. |
| `04_sensitivities.ipynb` | Greeks by Monte Carlo | Computes the delta of a capped call, $(S_T-K)^+\mathbf 1_{\{S_T\le 1.3K\}}$, by finite differences with common random numbers and by the likelihood ratio method. |


The exercise statements (énoncés) come from the course.

Course references:
- P. Glasserman, *Monte Carlo Methods in Financial Engineering*, Springer, 2004.
- E. Gobet, *Monte Carlo Methods and Stochastic Processes: From Linear to Non-Linear*, CRC Press, 2016.
- G. Pagès, *Numerical Probability: An Introduction with Applications to Finance*, Springer, 2018.
