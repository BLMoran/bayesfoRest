# Model Extraction Functions for bayesfoRest

S3 generics and methods for extracting posterior draws, pooled effect
summaries, and model type information from various Bayesian
meta-analysis backends. Supported backends: brms, MetaStan, bayesmeta,
rstanarm, rstan, cmdstanr, bayesma, RoBMA (brma), R2jags (JAGS).

All methods produce output in a canonical internal format so that
downstream plotting and table functions are backend-agnostic.
