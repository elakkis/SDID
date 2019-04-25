# SDID
Synthetic Difference in Difference

Implementation of the Synthetic Difference in Difference Method proposed by Arkhangelsky et al (2019)

Example code:


#simulate the data


#simulate the data
N <- 100  <br />
T <- 40 <br />
treat_effect <- 1 <br />
treated_units <- 5<br />
treat_date <- 20<br />
covariates <- NULL<br />
df <- simulate_panel(N, T, treat_effect, treated_units, treat_date)<br />

#estimate<br />
result <- estimate(df, covariates, treat_date, treated_units)<br />


#output the results<br />
print(result$est)<br />
print(result$sigma)<br />
print(result$CI_lower)<br />
print(result$CI_upper)<br />
