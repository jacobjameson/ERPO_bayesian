# diagnostics.R
# Run MCMC diagnostics on all JAGS fits using the superdiag package

# 1) Load required libraries
library(superdiag)   # for MCMC diagnostics
library(coda)        # for as.mcmc() conversions

# 2) Collect your fitted JAGS objects into a named list
models <- list(
  Model1 = jags1,
  Model2 = jags2,
  Model3 = jags3,
  Model4 = jags4,
  Model5 = jags5
)

# 3) Run diagnostics for each model
for (nm in names(models)) {
  cat("\n==============================\n")
  cat("Diagnostics for", nm, "\n")
  cat("==============================\n")
  
  # Convert the R2jags object to an mcmc.list
  mcmc_fit <- as.mcmc(models[[nm]])
  
  # Superdiag runs a battery of convergence and mixing tests and produces plots
  print(superdiag(
    mcmc          = mcmc_fit,
    plot          = TRUE,                              # show diagnostic plots
  ))
  
}

mcmc_fit <- as.mcmc(models[['Model5']])
superdiag(
  mcmc          = mcmc_fit,
  plot          = TRUE,                              # show diagnostic plots
)


superdiag(
  mcmc = mcmc_fit,
  plot = TRUE,
  parms = c("beta1","beta2","delta1","delta2","r","IRR_5_year")
)

pars <- c("beta1","beta2","delta1","delta2","r","IRR_5_year")
library(coda)
# assume jags5.mcmc is your mcmc.list from running jags()
mcmc_sub <- lapply(mcmc_fit, function(chain) {
  as.mcmc(chain[, pars])
})
mcmc_sub <- as.mcmc.list(mcmc_sub)

superdiag(
  mcmcoutput   = mcmc_sub,
  burnin       = 2000,    # match what you used in JAGS
  plot         = TRUE
)
