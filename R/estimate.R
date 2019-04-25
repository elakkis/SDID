#' Estimates the treatment effect using Synthetic Difference in Difference method
#'
#' @param df  Balanced panel data with the following columns: unit, time, treatment, Y (unit id, time id, treatment indicator and outcome variable)
#' @param covariates A vector of covariate names to be included. df must contain columns with names included in this variable
#' @param treatment_time  Time when the treatment occured
#' @param treated_units A vector of treated units' ids
#' @param weight_type : Type of weights to use. Allowed options are "SDID" for synthetic DID weights or "Kernel" for exponential kernel weighting
#' @param resampling_procedure : The algorithm used to obtain distribution of treatment effect. Allowed options are "Jackknife" and "Bootstrap"
#' @return A list that includes:
#' * est : Point estimate of Average Treatment effect on treated (ATT)
#' * sigma : Standard error of ATT estimate
#' * att_arr : N0 size array of resampled estimates of ATT using the resampling procedure
#' * CI_upper : 95 % confidence interval upper bound for the estimate
#' * CI_lower : 95 % confidence interval lower bound for the estimate
estimate <- function(df, covariates, treatment_time, treated_units, weight_type = "SDID", resampling_procedure = "Jackknife")
{
  #prepare the data
  data <- prepare_data(df, covariates, treatment_time, treated_units)
  
  #estimate treatment effect
  #prepare the weights
  dzeta <-  var(df$Y)
  rho <-  1
  tol <-  1e-10
  w_t <- time_weights(data$Y00, data$Y01, dzeta, rho, tol, weight_type)
  w_i <- unit_weights(data$Y00, data$Y10, dzeta, rho, tol, weight_type)
  
  #get point estimate
  est <- att(data$X, data$Y, -1, w_i, w_t)
  
  #inference:
  att_ests <- inference(data$X, data$Y, data$Y00, data$Y01, data$Y10, dzeta, rho, tol, weight_type, resampling_procedure)
  
  result <- NULL
  result$est <- est
  result$sigma <- sqrt(var(att_ests))
  result$att_arr <- att_ests
  result$CI_upper <- quantile(att_ests, 0.95)
  result$CI_lower <- quantile(att_ests, 0.05)
}
