#' Prepares the panel dataset and outputs object necessary for the estimation
#'
#' @param df  Balanced panel data with the following columns: unit, time, treatment, Y (unit id, time id, treatment indicator and outcome variable)
#' @param covariates A vector of covariate names to be included. df must contain columns with names included in this variable
#' @param treatment_time  Time when the treatment occured
#' @param treated_units A vector of treated units' ids
#' @return A list that includes:
#' * Y : N*T x 1 matrix of outcome variables
#' * X : Data matrix X used in WLS
#' * Y00 : T0 x N0 matrix of outcomes for control group before the treatment
#' * Y01 : T-T0 x N0 matrix of outcomes for control group after treatment
#' * Y10 : T0 x N - N0 matrix of outcomes for treated group before treatment
#' * n_covariates : number of covariates used in the estimation
prepare_data <- function(df, covariates, treatment_time, treated_units)
{
  #take in dataframe and construct matrices Y00, Y01, Y11
  #as well as the covariance matrix we will need for the estimation
  
  
  #check if all the necessary columns exist
  for (col in covariates) 
  {
    if(! col %in% colnames(df))
    {
      stop("The column does not exist, check the format of the dataframe")
    }
  }
  
  if (! "time" %in% colnames(df))
  {
    stop("time column does not exist in the dataframe")
  }
  if (! "unit" %in% colnames(df))
  {
    stop("unit column does not exist in the dataframe")
  }
  if (! "treatment" %in% colnames(df))
  {
    stop("treatment column does not exist in the dataframe")
  }
  if (! "Y" %in% colnames(df))
  {
    stop("Y column does not exist in the dataframe")
  }

  
  
  T <-  length(unique(df$time))
  T0 <-  length(unique(df[df$time < treatment_time,]$time))
  
  N <- length(unique(df$unit))
  N0 <- N - length(treated_units)
  
  
  #check if the panel is balanced
  if (dim(df)[1] != N*T)
  {
    stop("The panel is not balanced!")
  }

  
  #sort the dataframe, just in case:
  df <- df[order(df$unit, df$time),] 

  
  #first prepare the matrices
  Y00 <- df[(df$time < treatment_time) & (! df$unit %in% treated_units),]$Y;
  Y00 <- matrix(Y00, ncol = T0, byrow = TRUE)
  Y00 <- t(Y00)
  
  Y01 <- df[(df$time >= treatment_time) & (! df$unit %in% treated_units),]$Y;
  Y01 <- matrix(Y01, ncol = T - T0, byrow = TRUE)
  Y01 <- t(Y01)
  
  
  Y10 <- df[(df$time < treatment_time) & (df$unit %in% treated_units),]$Y;
  Y10 <- matrix(Y10, ncol = T0, byrow = TRUE)
  Y10 <- t(Y10)
  
  
  

  #now create time and unit dummies
  base_matr <- diag(1, N - 1, N - 1)
  unit_dummies <- matrix(rep(base_matr, each = T), nrow = (N - 1)*T, ncol = N - 1)
  unit_dummies <- rbind(unit_dummies, matrix(0, nrow = T, ncol = N - 1))

  base_matr <- rbind(diag(x = 1, nrow = T - 1, ncol = T - 1), matrix(0, 1, T-1))
  time_dummies <- do.call(rbind, replicate(N, base_matr, simplify = FALSE))

  
  #now create X matrix
  X <- as.matrix(cbind(df$treatment, df[,covariates], unit_dummies, time_dummies))
  dim(X)
  
  Y <- as.matrix(df$Y)
  
  
  result <- list()
  result$Y <- Y
  result$X <- X
  result$Y00 <- Y00
  result$Y01 <- Y01
  result$Y10 <- Y10
  result$n_covariates <- length(covariates)
  
  return(result)
}

