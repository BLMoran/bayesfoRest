#' Studies Comparing Medications on Post-Operative Nausea and Vomiting
#'
#' Outcome data from randomised controlled trials comparing anti-emetic medications (eg ondansetron, dexamethasone) against comparator anti-emetics or placebo in PONV outcomes for patients undergoing surgery.
#'
#' @format A data frame with 12 rows and 18 variables:
#' \describe{
#'   \item{Author}{Study author}
#'   \item{Year}{Year study published}
#'   \item{Subgroup}{Surgery subgroup}
#'   \item{Control}{Control intervention}
#'   \item{Intervention}{Intervention medication}
#'   \item{N_Total}{Total number of patients in the study}
#'   \item{N_Control}{Number of patients in the control group}
#'   \item{N_Intervention}{Number of patients in the intervention group}
#'   \item{Outcome_Control_Yes}{Number of patients in the control group with the outcome (PONV)}
#'   \item{Outcome_Control_No}{Number of patients in the control group without the outcome (PONV)}
#'   \item{Outcome_Intervention_Yes}{Number of patients in the intervention group with the outcome (PONV)}
#'   \item{Outcome_Intervention_No}{Number of patients in the intervetion group without the outcome (PONV)}
#'   \item{D1}{Domain 1 of the Risk of Bias 2.0 tool}
#'   \item{D2}{Domain 2 of the Risk of Bias 2.0 tool}
#'   \item{D3}{Domain 3 of the Risk of Bias 2.0 tool}
#'   \item{D4}{Domain 4 of the Risk of Bias 2.0 tool}
#'   \item{D5}{Domain 5 of the Risk of Bias 2.0 tool}
#'   \item{Overall}{Overall risk of bias}
#' }
"binary_outcome"