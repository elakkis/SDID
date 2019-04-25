#' Creates a simulated balanced panel data to evaluate the method
#'
#' @param N       Number of units in the dataset to be created
#' @param T Number of units in the dataset to be created
#' @param treat_effect    the magnitude of treatment effect
#' @param treated_units a vector of treated units
#' @param treat_date    date at which the treatment happens
#' @return A dataframe with the following columns
#' * unit : Unit id
#' * time z : Time id
#' * treatment : Indicator variable for treatment
#' * Y : Outcome variable
simulate_panel <- function(N, T, treat_effect, treated_units, treat_date)
{
  beta <- treat_effect
  
  # First generate our data using the time constant effects
  data <- data.frame(unit = 1:N, fe = rnorm(N))
  
  # We expand our data by nobs
  data <- data[rep(1:N, each = T),]
  
  #add time fixed effects
  time_fe <- rep(rnorm(T), N)
  data$te <- time_fe
  
  
  #add time index
  data$time <- matrix(rep(1:T, times = N))
  
  data$treatment <- (data$time >= treat_date)*(data$unit  %in% treated_units)
  
  # Finally we are ready to simulate our y variables
  data$Y <- beta * data$treatment + data$fe + data$te + rnorm(N*T)
  
  return(data[,c("unit", "time", "treatment", "Y")])
}