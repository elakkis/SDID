# SDID
Synthetic Difference in Difference

Implementation of the Synthetic Difference in Difference Method proposed by Arkhangelsky et al (2019)

Example code:


#simulate the data


#simulate the data
N <- 100 \n
T <- 40
treat_effect <- 1
treated_units <- 5
treat_date <- 20
covariates <- NULL
df <- simulate_panel(N, T, treat_effect, treated_units, treat_date)

#estimate
result <- estimate(df, covariates, treat_date, treated_units)


#output the results
print(result$est)
print(result$sigma)
print(result$CI_lower)
print(result$CI_upper)
