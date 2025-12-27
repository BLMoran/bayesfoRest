# Introduction to Creating Forest Plots for Bayesian Meta Analyses

``` r
library(metafor)
library(brms)
library(bayesfoRest)

# Load data
binary_outcome <- binary_outcome

# Prepare with {metafor}
binary_outcome <- escalc(
  measure="OR", 
  ai=Event_Intervention, bi=Outcome_Intervention_No,
  ci=Event_Control, di=Outcome_Control_No,
  n1i=N_Intervention, n2i=N_Control,  
  data=binary_outcome) |> 
  dplyr::mutate(sei = sqrt(vi))

# Run the brms model
model_bin <- brms::brm(yi | se(sei) ~ 1 + (1 | Author),
             family = gaussian(),
             data = binary_outcome,
             prior = brms::prior(normal(0, 1), class = Intercept) +
                     brms::prior(cauchy(0, 0.5), class = sd),
             iter = 4000,
             cores = parallel::detectCores(),
             control = list(adapt_delta = 0.99),
             backend = "cmdstanr",
             silent = 2)

# Create forest plot
forest_plot_or <- bayes_forest(
  model = model_bin,
  data = binary_outcome,
  measure = "OR",
  studyvar = Author,
  year = Year,
  c_n = N_Control,
  i_n = N_Intervention,
  c_event = Event_Control,
  i_event = Event_Intervention,
  title = "Treatment Effect on Binary Outcome",
  subtitle = "Bayesian Random-Effects Meta-Analysis (Odds Ratio)",
  xlim = c(0.1, 3.5)
)

# Display Plot
forest_plot_or
```

![](reference/figures/forest_plot_binary.png)
