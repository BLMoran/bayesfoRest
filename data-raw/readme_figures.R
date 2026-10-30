# Regenerates the figures in man/figures/ that README.Rmd displays.
# Run from the package root: source("data-raw/readme_figures.R")

devtools::load_all()
library(metafor)
library(brms)

fig <- function(plot, file, width = 16, height = 9) {
  ggplot2::ggsave(file.path("man/figures", file), plot,
                  width = width, height = height, dpi = 300, bg = "white")
}

# Data and model --------------------------------------------------------------

binary_es <- escalc(
  measure = "OR",
  ai = Event_Intervention, bi = Outcome_Intervention_No,
  ci = Event_Control, di = Outcome_Control_No,
  n1i = N_Intervention, n2i = N_Control,
  data = binary_outcome) |>
  dplyr::mutate(sei = sqrt(vi))

model_bin <- brm(
  yi | se(sei) ~ 1 + (1 | Author),
  family = gaussian(),
  data = binary_es,
  prior = prior(normal(0, 1), class = Intercept) +
          prior(cauchy(0, 0.5), class = sd),
  iter = 4000,
  cores = 4,
  control = list(adapt_delta = 0.99),
  seed = 1234,
  silent = 2,
  refresh = 0,
  backend = "cmdstanr")

model_bin_ce <- brm(
  yi | se(sei) ~ 1,
  family = gaussian(),
  data = binary_es,
  prior = prior(normal(0, 1), class = Intercept),
  iter = 4000,
  cores = 4,
  control = list(adapt_delta = 0.99),
  seed = 1234,
  silent = 2,
  refresh = 0,
  backend = "cmdstanr")

forest_args <- list(
  model = model_bin,
  data = binary_es,
  measure = "OR",
  studyvar = quote(Author),
  yearvar = quote(Year),
  n_ctrl = quote(N_Control),
  n_int = quote(N_Intervention),
  event_ctrl = quote(Event_Control),
  event_int = quote(Event_Intervention),
  xlim = c(0.1, 3.5),
  add_pred = quote(TRUE))

forest <- function(...) do.call(bayes_forest, utils::modifyList(forest_args, list(...)))

# Forest plots ----------------------------------------------------------------

fig(forest(), "forest_plot_binary.png")
fig(forest(model = model_bin_ce, add_pred = FALSE), "forest_plot_common_effect.png")
fig(forest(plot_output = "boxplot"), "forest_plot_boxplot.png")
fig(forest(subgroup = TRUE, subgroupvar = quote(Subgroup)),
    "forest_plot_subgroup.png", height = 12)
fig(forest(add_rob = TRUE, rob_tool = "rob2"),
    "forest_plot_rob.png")
fig(forest(add_rob = TRUE,
           add_rob_legend = TRUE, rob_tool = "rob2"),
    "forest_plot_rob_legend.png", width = 22, height = 12)
fig(forest(add_rob = TRUE, rob_tool = "rob2", exclude_high_rob = TRUE),
    "forest_excl_high_rob.png")

# Risk of bias ----------------------------------------------------------------

fig(rob_plot(binary_outcome, rob_tool = "rob2", add_rob_legend = TRUE),
    "rob_viz.png", width = 12, height = 6)

fig(rob_summary(binary_outcome, rob_tool = "rob2", title = "Risk of Bias Summary"),
    "rob_summary.png", width = 10, height = 6)

# Overall and funnel -----------------------------------------------------------

fig(overall_plot(model_bin, binary_es, measure = "OR", add_tau = TRUE,
                 add_pred = TRUE, incl_mu_prior = TRUE, incl_tau_prior = TRUE,
                 add_null_range = TRUE, mu_xlim = c(0.3, 1.5), tau_xlim = c(0, 1),
                 plot_arrangement = "vertical"), "overall_plot.png")

fig(funnel_plot(model_bin, binary_es, measure = "OR",
                studyvar = Author, yearvar = Year,
                yi_col = "yi", sei_col = "sei", label_studies = TRUE),
    "funnel_plot.png", width = 6, height = 5)

# Sensitivity -----------------------------------------------------------------

priors <- list(
  vague = list(name = "Vague",
               mu_prior = bf_normal(0, 10), tau_prior = bf_half_cauchy(0, 1)),
  weakreg = list(name = "Weakly Regularising",
                 mu_prior = bf_normal(0, 1), tau_prior = bf_half_cauchy(0, 0.5)),
  informative = list(name = "Informative",
                     mu_prior = bf_normal(0, 0.5), tau_prior = bf_half_cauchy(0, 0.3)))

sens_args <- list(model = model_bin, data = binary_es, priors = priors,
                  measure = "OR", rob_var = quote(Overall),
                  exclude_high_rob = TRUE, incl_common_effect = TRUE,
                  xlim = c(0.25, 1.5),
                  add_null_range = TRUE)

fig(do.call(sensitivity_plot, c(sens_args, add_probs = FALSE)),
    "sens_plot_no_prob.png", width = 18, height = 10)
fig(do.call(sensitivity_plot, c(sens_args, add_probs = TRUE)),
    "sens_plot_prob.png", width = 22, height = 10)

fig(ecdf_model_plot(model_bin, binary_es, priors = priors["weakreg"],
                    measure = "OR", rob_var = Overall, exclude_high_rob = TRUE,
                    incl_common_effect = TRUE,
                    xlim = c(0.25, 1.5), add_null_range = TRUE),
    "ecdf_plot.png", width = 10, height = 6)

fig(ecdf_prior_plot(model_bin, binary_es, priors = priors, measure = "OR",
                    xlim = c(0.25, 1.5), add_null_range = TRUE),
    "ecdf_prior_plot.png", width = 10, height = 6)
