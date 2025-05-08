################################################################################
# erpo_analysis_working_jags.R
# Complete workflow following a previously successful pattern
################################################################################

# 0) Load libraries
library(R2jags)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1) Load data
if (!exists("analysis_data")) {
  load("manuscript/src/erpo_analysis_data.RData")
}
set.seed(123)

# 2) Create state and year factors
state_factor <- as.numeric(factor(analysis_data$state))
year_factor  <- as.numeric(factor(analysis_data$year))

# 3) Build year dummy matrix
year_dummies <- model.matrix(~ factor(year) - 1, data = analysis_data)
colnames(year_dummies) <- paste0("year_", sort(unique(analysis_data$year)))
# drop first to avoid collinearity
year_dummies <- year_dummies[, -1]

# 4) Build lagged year dummies
n_obs <- nrow(analysis_data)
n_years <- ncol(year_dummies)
year_dummies_lag1 <- matrix(0, n_obs, n_years)
year_dummies_lag2 <- matrix(0, n_obs, n_years)

for (s in unique(state_factor)) {
  idx <- which(state_factor == s)
  if (length(idx) > 1) {
    year_dummies_lag1[idx[-1], ] <- year_dummies[idx[-length(idx)], ]
  }
  if (length(idx) > 2) {
    year_dummies_lag2[idx[-(1:2)], ] <- year_dummies[idx[1:(length(idx)-2)], ]
  }
}

# 5) Create lagged policy variables directly rather than matrices
erpo_lag1 <- analysis_data %>%
  group_by(state) %>%
  arrange(state, year) %>%
  mutate(erpo_lag1 = lag(erpo, 1, default = 0)) %>%
  ungroup() %>%
  pull(erpo_lag1)

erpo_phase_in_lag1 <- analysis_data %>%
  group_by(state) %>%
  arrange(state, year) %>%
  mutate(erpo_phase_in_lag1 = lag(erpo_phase_in, 1, default = 0)) %>%
  ungroup() %>%
  pull(erpo_phase_in_lag1)

erpo_lag2 <- analysis_data %>%
  group_by(state) %>%
  arrange(state, year) %>%
  mutate(erpo_lag2 = lag(erpo, 2, default = 0)) %>%
  ungroup() %>%
  pull(erpo_lag2)

erpo_phase_in_lag2 <- analysis_data %>%
  group_by(state) %>%
  arrange(state, year) %>%
  mutate(erpo_phase_in_lag2 = lag(erpo_phase_in, 2, default = 0)) %>%
  ungroup() %>%
  pull(erpo_phase_in_lag2)

# 6) Add covariates
# HFR as primary covariate
hfr <- analysis_data$HFR

# Secondary covariates (LASSO shrinkage)
covariates <- c("under_100_percent","x100_199_percent","x200_399_percent","x400_percent",
                "white","black","hispanic","asian",
                "adults_19_25","adults_26_34","adults_35_54","adults_55_64","x65")
cov_matrix <- as.matrix(analysis_data[, covariates])

# 7) Assemble JAGS data
jags_data <- list(
  n = nrow(analysis_data),
  n_states = length(unique(state_factor)),
  n_years = ncol(year_dummies),
  n_covariates = length(covariates),
  
  deaths = analysis_data$deaths,
  offset = log(analysis_data$population),
  log_rate_lag1 = analysis_data$log_rate_lag1,
  log_rate_lag2 = analysis_data$log_rate_lag2,
  
  state = state_factor,
  year_dummy = year_dummies,
  
  erpo = analysis_data$erpo,
  erpo_phase_in = analysis_data$erpo_phase_in,
  erpo_lag1 = erpo_lag1,
  erpo_phase_in_lag1 = erpo_phase_in_lag1,
  erpo_lag2 = erpo_lag2,
  erpo_phase_in_lag2 = erpo_phase_in_lag2,
  
  hfr = hfr,
  covariates = cov_matrix,
  
  prior_sd = 0.064
)

# 8) Define improved JAGS model
jags_model <- "model {
  # Likelihood
  for (i in 1:n) {
    log_mu[i] <- offset[i] +
                  delta1 * log_rate_lag1[i] +
                  delta2 * log_rate_lag2[i] +
                  beta1 * erpo[i] +
                  beta2 * erpo_phase_in[i] +
                  gamma_hfr * hfr[i] +
                  inprod(covariates[i,], gamma_cov) +
                  state_effect[state[i]] +
                  inprod(year_dummy[i,], year_effect) -
                  delta1 * (beta1 * erpo_lag1[i] + beta2 * erpo_phase_in_lag1[i]) -
                  delta2 * (beta1 * erpo_lag2[i] + beta2 * erpo_phase_in_lag2[i])
                 
    # Negative binomial likelihood with dispersion parameter r
    p[i] <- r / (r + exp(log_mu[i]))
    deaths[i] ~ dnegbin(p[i], r)
  }
  
  # Priors - improved from original model
  # Use normal priors for AR coefficients with appropriate constraints
  delta1_raw ~ dnorm(0.5, 4)  # More informative prior centered at 0.5
  delta1 <- max(0, min(delta1_raw, 1))  # Constrain between 0 and 1
  
  delta2_raw ~ dnorm(0, 4)    # Prior centered at 0
  delta2 <- max(-1, min(delta2_raw, 1))  # Constrain between -1 and 1
  
  # Policy effect priors
  beta1 ~ dnorm(0, tau_beta)  # Immediate ERPO effect
  beta2 ~ dnorm(0, tau_beta)  # Phase-in ERPO effect
  
  # More informative prior for HFR
  gamma_hfr ~ dnorm(0, 25)    # SD = 0.2, tighter than original
  
  # Better parameterization for dispersion
  log_r ~ dnorm(5.5, 0.5)     # Log-normal prior centered around ~250
  r <- exp(log_r)
  
  # LASSO priors for covariates - using improved hyperprior
  for (j in 1:n_covariates) {
    gamma_cov[j] ~ ddexp(0, tau_lasso)
  }
  
  # Half-normal prior for tau_lasso instead of half-Cauchy
  tau_lasso_raw ~ dnorm(0, 1) T(0,)
  tau_lasso <- tau_lasso_raw
  
  # Random effects for states with proper sum-to-zero constraint
  for (s in 1:n_states) {
    state_effect_raw[s] ~ dnorm(0, tau_state)
  }
  state_effect_mean <- mean(state_effect_raw[1:n_states])
  for (s in 1:n_states) {
    state_effect[s] <- state_effect_raw[s] - state_effect_mean
  }
  
  # Fixed effects for years
  for (t in 1:n_years) {
    year_effect[t] ~ dnorm(0, tau_year)
  }
  
  # Hyperpriors
  tau_beta <- pow(prior_sd, -2)
  tau_state <- pow(0.5, -2)
  tau_year <- pow(0.25, -2)
  
  # Calculate the total effect at 5 years
  effect_5_year <- beta1 + beta2
  IRR_5_year <- exp(effect_5_year)
  
  # Calculate effects over time
  for (t in 1:11) {
    year_since[t] <- t - 1
    
    # Phase-in over 5 years
    phase[t] <- min(year_since[t]/5, 1)
    effect[t] <- beta1 + beta2 * phase[t]
    IRR[t] <- exp(effect[t])
  }
}"

# 9) Parameters to monitor
params <- c("beta1", "beta2", "delta1", "delta2", "gamma_hfr", "gamma_cov", "r",
            "effect_5_year", "IRR_5_year", "effect", "IRR", "year_since", "phase")

# 10) Create improved initialization function for better mixing
inits <- function(chain) {
  # Create chain-specific initialization
  # This helps with mixing by starting chains in different regions
  list(
    delta1_raw = 0.5 + 0.05 * (chain - 2),     # Different starting points for each chain
    delta2_raw = 0.1 + 0.05 * (chain - 2),
    beta1 = -0.01 * chain,                     # Small negative values varying by chain
    beta2 = -0.01 * (4 - chain),               # Small negative values in reverse order
    gamma_hfr = 0.05 * (chain - 2),            # Centered around 0
    gamma_cov = rep(0, jags_data$n_covariates),
    log_r = 5.5 + 0.1 * (chain - 2),           # Centered on log(250)
    tau_lasso_raw = 0.8 + 0.1 * chain,
    state_effect_raw = rep(0, jags_data$n_states),
    year_effect = rep(0, jags_data$n_years)
  )
}

# 11) Run JAGS model with improved settings
# Run a test first with one chain to ensure it compiles
cat("Running JAGS test with 1 chain...\n")
jags_test <- jags(
  data = jags_data,
  parameters.to.save = c("beta1", "beta2", "delta1", "delta2", "r"),
  model.file = textConnection(jags_model),
  inits = function() { inits(1) },
  n.chains = 1,
  n.iter = 1000,
  n.burnin = 500,
  n.thin = 2
)

# If test successful, run full model with improved settings
if (exists("jags_test")) {
  cat("Test successful. Running full JAGS model with improved settings...\n")
  
  jags_full <- jags(
    data = jags_data,
    parameters.to.save = params,
    model.file = textConnection(jags_model),
    inits = inits,
    n.chains = 3,
    n.iter = 30000,        # Increased iterations
    n.burnin = 10000,      # Increased burn-in
    n.thin = 10,           # Increased thinning
    DIC = TRUE
  )
  
  # 12) Summarize results
  print(jags_full)
  
  # Add convergence diagnostics
  cat("\nGelman-Rubin statistics for key parameters:\n")
  gelman.diag(jags_full$BUGSoutput$sims.array[, , c("beta1", "beta2", "delta1", "delta2", "r")])
  
  # Continue with the rest of your analysis...
  # 13) Extract IRR results
  samps <- jags_full$BUGSoutput$sims.matrix
  irr5 <- samps[, "IRR_5_year"]
  
  cat("Median IRR @ 5 years:", median(irr5), "\n")
  cat("95% CI:", quantile(irr5, c(0.025, 0.975)), "\n")
  cat("Probability IRR < 1:", mean(irr5 < 1), "\n")
  
  # 14) Create data frame for IRR over time
  years <- 0:10
  n_years_plot <- length(years)
  
  irr_summary <- data.frame(
    year = years,
    median = numeric(n_years_plot),
    lower = numeric(n_years_plot),
    upper = numeric(n_years_plot),
    prob_less_than_1 = numeric(n_years_plot)
  )
  
  for (i in 1:n_years_plot) {
    irr_col <- paste0("IRR[", i, "]")
    irr_summary$median[i] <- median(samps[, irr_col])
    irr_summary$lower[i] <- quantile(samps[, irr_col], 0.025)
    irr_summary$upper[i] <- quantile(samps[, irr_col], 0.975)
    irr_summary$prob_less_than_1[i] <- mean(samps[, irr_col] < 1)
  }
  
  # Print IRR summary
  print(irr_summary)
}





################################################################################
# ERPO Effects Visualization
# Creates trace plot and density plots based on JAGS model results
################################################################################

# Load necessary packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) # For combining plots

# Define color palette
journal_colors <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00")

#---------------------------------------------------------
# 1. Extract results from jags_full model
#---------------------------------------------------------

# Extract MCMC samples and create summary stats
n_samples <- nrow(jags_full$BUGSoutput$sims.matrix)
samples <- jags_full$BUGSoutput$sims.matrix

# Get IRR samples for each time point
irr_samples <- data.frame(
  sample_id = rep(1:n_samples, 11),
  year = rep(0:10, each = n_samples),
  irr = NA
)

# Fill with IRR values
for (i in 1:11) {
  param_name <- paste0("IRR[", i, "]")
  irr_samples$irr[irr_samples$year == (i-1)] <- samples[, param_name]
}

# Calculate summary statistics
irr_summary <- irr_samples %>%
  group_by(year) %>%
  summarize(
    median = median(irr),
    mean = mean(irr),
    lower_95 = quantile(irr, 0.025),
    upper_95 = quantile(irr, 0.975),
    prob_less_than_1 = mean(irr < 1)
  )

#---------------------------------------------------------
# 2. Create trace plot with MCMC samples
#---------------------------------------------------------

# Create trace plot showing individual MCMC samples
trace_plot <- ggplot() +
  # Add individual traces (sample a subset to avoid overplotting)
  geom_line(
    data = irr_samples %>% filter(sample_id %% 20 == 0), # Sample every 20th MCMC sample
    aes(x = year, y = irr, group = sample_id),
    color = journal_colors[3], alpha = 0.09, size = 0.3
  ) +
  # Add median line
  geom_line(
    data = irr_summary,
    aes(x = year, y = median),
    color = 'black', size = 1.2
  ) +
  # Add reference line at IRR = 1
  geom_hline(
    yintercept = 1, 
    linetype = "longdash", 
    color = "gray30", 
    size = 0.8
  ) +
  annotate(
    "text", x = 4.2, y = 1.02, 
    label = "No effect (IRR = 1)",
    hjust = 0, vjust = 0, 
    size = 3, color = "black", 
    fontface = "italic"
  ) +
  # Labels and scales
  labs(
    title = "Effect of ERPO Laws Over Time",
    subtitle = "Individual MCMC traces with median effect highlighted",
    x = "Years Since Implementation",
    y = "Incidence Rate Ratio (IRR)"
  ) +
  scale_x_continuous(breaks = 0:10, minor_breaks = NULL) +
  scale_y_continuous(
    breaks = seq(0.9, 1.1, 0.05),
  ) +
  # Apply theme
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, color = "gray30", margin = margin(b = 15)),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.text = element_text(size = 10, color = "gray20"),
    axis.ticks = element_line(color = "gray70", size = 0.3),
    panel.grid.major = element_line(color = "gray92", size = 0.4),
    panel.grid.minor = element_line(color = "gray96", size = 0.2),
    plot.margin = margin(20, 25, 20, 20)
  )

trace_plot
#---------------------------------------------------------
# 3. Create density plot for selected years
#---------------------------------------------------------

# Select specific years to show distributions
selected_years <- c(1, 2, 3, 4, 5)
density_data <- irr_samples %>%
  filter(year %in% selected_years) %>%
  mutate(Year = factor(paste("Year", year)))

# Create the density plot
density_plot <- ggplot(
  density_data, 
  aes(x = irr, fill = Year, color = Year)
) +
  geom_density(alpha = 0.65, size = 0.7) +
  geom_rug(alpha = 0.7, show.legend = FALSE, size = 0.4) +
  geom_vline(xintercept = 1, linetype = "longdash", color = "gray30", size = 0.8) +
  annotate(
    "text", x = 1.01, y = 22, 
    label = "No effect (IRR = 1)",
    hjust = 0, vjust = 0, size = 3, 
    color = "black", 
    fontface = "italic"
  ) +
  labs(
    title = "Posterior Distributions of ERPO Effect by Year",
    subtitle = "Incidence Rate Ratio distributions for years 1-5 after implementation",
    x = "Incidence Rate Ratio (IRR)",
    y = "Density"
  ) +
  scale_x_continuous(
    breaks = seq(0.9, 1.1, by = 0.05),
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = journal_colors[1:5], name = "") +
  scale_color_manual(values = journal_colors[1:5], name = "") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, color = "gray30", margin = margin(b = 15)),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.text = element_text(size = 10, color = "gray20"),
    axis.ticks = element_line(color = "gray70", size = 0.3),
    panel.grid.major = element_line(color = "gray92", size = 0.4),
    panel.grid.minor = element_line(color = "gray96", size = 0.2),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1.2, "lines"),
    legend.margin = margin(t = 5, r = 0, b = 5, l = 5),
    legend.box.margin = margin(l = 10),
    plot.margin = margin(20, 25, 20, 20)
  ) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.8, size = 4), ncol = 1),
    color = FALSE  # Hide color legend (redundant with fill)
  )
density_plot
#---------------------------------------------------------
# 4. Combine plots and save
#---------------------------------------------------------

# Combine both plots vertically
combined_plot <- trace_plot / density_plot

# Display the combined plot
print(combined_plot)

# Save the combined plot
ggsave("erpo_effects_combined.png", combined_plot, width = 12, height = 10, dpi = 300)

#---------------------------------------------------------
# 5. Counterfactual analysis (deaths prevented)
#---------------------------------------------------------

# Build a lookup table of IRRs by years-since
irr_lookup <- irr_summary %>%
  select(years_since = year, IRR = median)

# For each state-year, compute years_since and join on IRR
prevented_df <- analysis_data %>%
  # Get each state's first ERPO year
  group_by(state) %>%
  mutate(
    impl_year = if_else(any(erpo > 0), min(year[erpo > 0]), NA_real_),
    years_since = if_else(!is.na(impl_year), year - impl_year, NA_real_)
  ) %>%
  ungroup() %>%
  # Cap years_since between 0 and 10
  mutate(
    years_since = pmin(pmax(years_since, 0), 10)
  ) %>%
  # Join to get the median IRR for that "years_since"
  left_join(irr_lookup, by = "years_since") %>%
  # Where there was no ERPO yet IRR=1 (no effect)
  mutate(IRR = if_else(is.na(IRR), 1, IRR)) %>%
  # Compute counterfactual and prevented deaths
  mutate(
    deaths_cf = deaths / IRR,
    deaths_saved = deaths_cf - deaths
  )

# Summarize over the full sample
prevention_summary <- prevented_df %>%
  summarize(
    total_observed = sum(deaths),
    total_counterfactual = sum(deaths_cf),
    total_prevented = sum(deaths_saved)
  )

# Print prevention summary
print(prevention_summary)



























# Load necessary libraries
library(R2jags)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(patchwork)

# Load the data if not already loaded
if(!exists("analysis_data")) {
  load("manuscript/src/erpo_analysis_data.RData")
}

# Set seed for reproducibility
set.seed(123)

# -----------------------------------------------------
# 1. Prepare data for the AR model
# -----------------------------------------------------

# Create factors for states and years
state_factor <- as.numeric(factor(analysis_data$state))
year_factor <- as.numeric(factor(analysis_data$year))

# Create time indicator matrices for years
year_dummies <- model.matrix(~ factor(year) - 1, data = analysis_data)
colnames(year_dummies) <- paste0("year_", unique(analysis_data$year))

# Remove the first column as reference
year_dummies <- year_dummies[, -1]

# Create lagged year indicator matrices
year_dummies_lag1 <- matrix(0, nrow = nrow(analysis_data), ncol = ncol(year_dummies))
year_dummies_lag2 <- matrix(0, nrow = nrow(analysis_data), ncol = ncol(year_dummies))

# Fill in lagged year indicators
for (s in unique(state_factor)) {
  state_idx <- which(state_factor == s)
  if (length(state_idx) > 1) {
    year_dummies_lag1[state_idx[-1], ] <- year_dummies[state_idx[-(length(state_idx))], ]
  }
  if (length(state_idx) > 2) {
    year_dummies_lag2[state_idx[-(1:2)], ] <- year_dummies[state_idx[1:(length(state_idx)-2)], ]
  }
}

# -----------------------------------------------------
# 2. Define and run the autoregressive JAGS model
# -----------------------------------------------------

# Define the JAGS model
# Define the JAGS model
jags_model5 <- "
model {
  # Likelihood
  for (i in 1:n) {
    log_mu[i] <- offset[i] + 
                 delta1 * log_rate_lag1[i] + 
                 delta2 * log_rate_lag2[i] + 
                 beta1 * erpo[i] + 
                 beta2 * erpo_phase_in[i] + 
                 state_effect[state[i]] + 
                 inprod(year_dummy[i,], year_effect) - 
                 delta1 * (beta1 * erpo_lag1[i] + beta2 * erpo_phase_in_lag1[i]) -
                 delta2 * (beta1 * erpo_lag2[i] + beta2 * erpo_phase_in_lag2[i])
                 
    # Negative binomial likelihood with dispersion parameter r
    p[i] <- r / (r + exp(log_mu[i]))
    deaths[i] ~ dnegbin(p[i], r)
  }
  
  # Priors
  delta1 ~ dunif(0, 1)          # First-order AR coefficient
  delta2 ~ dunif(-1, 1)         # Second-order AR coefficient
  beta1 ~ dnorm(0, tau_beta)    # Immediate ERPO effect
  beta2 ~ dnorm(0, tau_beta)    # Phase-in ERPO effect
  r ~ dgamma(0.1, 0.1)          # Dispersion parameter
  
  # Random effects for states with proper sum-to-zero constraint
  for (s in 1:n_states) {
    state_effect_raw[s] ~ dnorm(0, tau_state)
  }
  state_effect_mean <- sum(state_effect_raw) / n_states
  for (s in 1:n_states) {
    state_effect[s] <- state_effect_raw[s] - state_effect_mean
  }
  
  # Fixed effects for years
  for (t in 1:n_years) {
    year_effect[t] ~ dnorm(0, tau_year)
  }
  
  # Hyperpriors
  tau_beta <- pow(prior_sd, -2)
  tau_state <- pow(0.5, -2)
  tau_year <- pow(0.25, -2)
  
  # Calculate the total effect at 5 years
  effect_5_year <- beta1 + beta2
  IRR_5_year <- exp(effect_5_year)
  
  # Calculate effects over time
  for (t in 1:11) {
    year_since[t] <- t - 1
    
    # Phase-in over 5 years
    phase[t] <- min(year_since[t]/5, 1)
    effect[t] <- beta1 + beta2 * phase[t]
    IRR[t] <- exp(effect[t])
  }
}
"

# Prepare data for JAGS model 5
jags_data5 <- list(
  n = nrow(analysis_data),
  n_states = length(unique(state_factor)),
  n_years = ncol(year_dummies),
  deaths = analysis_data$deaths,
  offset = log(analysis_data$population),
  log_rate_lag1 = analysis_data$log_rate_lag1,
  log_rate_lag2 = analysis_data$log_rate_lag2,
  state = state_factor,
  year_dummy = year_dummies,
  erpo = analysis_data$erpo,
  erpo_phase_in = analysis_data$erpo_phase_in,
  erpo_lag1 = analysis_data %>% 
    group_by(state) %>% 
    arrange(state, year) %>% 
    mutate(erpo_lag1 = lag(erpo, 1, default = 0)) %>% 
    ungroup() %>% 
    pull(erpo_lag1),
  erpo_phase_in_lag1 = analysis_data %>% 
    group_by(state) %>% 
    arrange(state, year) %>% 
    mutate(erpo_phase_in_lag1 = lag(erpo_phase_in, 1, default = 0)) %>% 
    ungroup() %>% 
    pull(erpo_phase_in_lag1),
  erpo_lag2 = analysis_data %>% 
    group_by(state) %>% 
    arrange(state, year) %>% 
    mutate(erpo_lag2 = lag(erpo, 2, default = 0)) %>% 
    ungroup() %>% 
    pull(erpo_lag2),
  erpo_phase_in_lag2 = analysis_data %>% 
    group_by(state) %>% 
    arrange(state, year) %>% 
    mutate(erpo_phase_in_lag2 = lag(erpo_phase_in, 2, default = 0)) %>% 
    ungroup() %>% 
    pull(erpo_phase_in_lag2),
  prior_sd = 0.064  # Same as in Stan model
)

# Parameters to monitor
params5 <- c("beta1", "beta2", "delta1", "delta2", "r", "effect_5_year", 
             "IRR_5_year", "effect", "IRR", "year_since", "phase")

# Run JAGS model 5
cat("Running Model 5: Autoregressive Negative Binomial with debiasing terms...\n")
jags5 <- jags(
  data = jags_data5,
  parameters.to.save = params5,
  model.file = textConnection(jags_model5),
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 2000,
  n.thin = 5
)

# Extract results
results5 <- list(
  beta1 = jags5$BUGSoutput$sims.matrix[, "beta1"],
  beta2 = jags5$BUGSoutput$sims.matrix[, "beta2"],
  delta1 = jags5$BUGSoutput$sims.matrix[, "delta1"],
  delta2 = jags5$BUGSoutput$sims.matrix[, "delta2"],
  r = jags5$BUGSoutput$sims.matrix[, "r"],
  effect_5_year = jags5$BUGSoutput$sims.matrix[, "effect_5_year"],
  IRR_5_year = jags5$BUGSoutput$sims.matrix[, "IRR_5_year"]
)

# Calculate summary statistics
results5_summary <- list(
  median_irr = median(results5$IRR_5_year),
  mean_irr = mean(results5$IRR_5_year),
  sd_irr = sd(results5$IRR_5_year),
  lower_95 = quantile(results5$IRR_5_year, 0.025),
  upper_95 = quantile(results5$IRR_5_year, 0.975),
  lower_80 = quantile(results5$IRR_5_year, 0.1),
  upper_80 = quantile(results5$IRR_5_year, 0.9),
  prob_reduction = mean(results5$IRR_5_year < 1)
)

# Add to results table
model5_row <- data.frame(
  Model = "Autoregressive NB with Debiasing",
  Median_IRR = round(results5_summary$median_irr, 3),
  Mean_IRR = round(results5_summary$mean_irr, 3),
  SD = round(results5_summary$sd_irr, 3),
  Lower_95CI = round(results5_summary$lower_95, 3),
  Upper_95CI = round(results5_summary$upper_95, 3),
  Lower_80CI = round(results5_summary$lower_80, 3),
  Upper_80CI = round(results5_summary$upper_80, 3),
  Prob_Reduction = round(results5_summary$prob_reduction, 3)
)

# Print results for model 5
cat("\n=== Model 5: Autoregressive NB with Debiasing ===\n")
cat(sprintf("IRR at 5 years: %.3f (95%% CI: %.3f to %.3f)\n", 
            results5_summary$median_irr, 
            results5_summary$lower_95, 
            results5_summary$upper_95))
cat(sprintf("Probability of reduction: %.1f%%\n", results5_summary$prob_reduction * 100))
cat(sprintf("AR coefficients: delta1 = %.3f, delta2 = %.3f\n", 
            median(results5$delta1), 
            median(results5$delta2)))
cat(sprintf("DIC: %.2f\n\n", jags5$BUGSoutput$DIC))

# Extract time-varying effects
time_effects <- data.frame(
  Year = 0:10,
  IRR = NA,
  Lower_95CI = NA,
  Upper_95CI = NA
)

for (t in 1:11) {
  irr_samples <- jags5$BUGSoutput$sims.matrix[, paste0("IRR[", t, "]")]
  time_effects$IRR[t] <- median(irr_samples)
  time_effects$Lower_95CI[t] <- quantile(irr_samples, 0.025)
  time_effects$Upper_95CI[t] <- quantile(irr_samples, 0.975)
}

# Plot time-varying effects
ggplot(time_effects, aes(x = Year, y = IRR, ymin = Lower_95CI, ymax = Upper_95CI)) +
  geom_line(size = 1) +
  geom_ribbon(alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgray") +
  labs(
    title = "Effect of ERPO Laws Over Time",
    subtitle = "From Autoregressive Model with Debiasing",
    x = "Years Since Implementation",
    y = "Incidence Rate Ratio (IRR)"
  ) +
  theme_minimal()

# Save plot

# Save full model 5 results
results5_full <- list(
  jags_model = jags5,
  summary = results5_summary,
  time_effects = time_effects
)

# Return the results for potential further analysis
results5_full



library(dplyr)

# 1. Build a small lookup table of IRRs by years-since (1..11)
irr_lookup <- time_effects %>%
  mutate(years_since = 0:10) %>%            # 0 through 10
  select(years_since, IRR)

# 2. For each state-year, compute years_since and join on IRR
prevented_df <- analysis_data %>%
  # get each state’s first ERPO year
  group_by(state) %>%
  mutate(
    impl_year = if_else(any(erpo>0), min(year[erpo>0]), NA_real_),
    years_since = if_else(!is.na(impl_year), year - impl_year, NA_real_)
  ) %>%
  ungroup() %>%
  # cap years_since between 0 and 10, and drop pre‐ERPO years
  mutate(
    years_since = pmin(pmax(years_since, 0), 10)
  ) %>%
  # join to get the median IRR for that “years_since”
  left_join(irr_lookup, by = "years_since") %>%
  # where there was no ERPO yet IRR=1 (no effect)
  mutate(IRR = if_else(is.na(IRR), 1, IRR)) %>%
  # 3. compute counterfactual and prevented
  mutate(
    deaths_cf    = deaths / IRR,
    deaths_saved = deaths_cf - deaths
  )

# 4. Summarize over the full sample
prevented_df %>%
  summarize(
    total_observed     = sum(deaths),
    total_counterfactual = sum(deaths_cf),
    total_prevented    = sum(deaths_saved)
  ) -> prevention_summary

print(prevention_summary)



# Load necessary packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) # For combining plots

# Define color palette to match your style
journal_colors <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00")

# 1. TRACE PLOT WITH INDIVIDUAL MCMC SAMPLES
# ------------------------------------------
# Extract samples for each time point (Years 0-10)
# We'll create a data frame with n_samples rows for each year

# First, extract MCMC samples from your JAGS output
# Assuming jags5 is your model output from the code snippet
n_samples <- nrow(jags5$BUGSoutput$sims.matrix)
n_years <- 7

# Create empty data frame to store all samples
trace_data <- data.frame(
  sample_id = rep(1:n_samples, n_years),
  year = rep(0:6, each = n_samples),
  irr = NA
)

# Fill with IRR values from JAGS output
for(i in 1:n_years) {
  param_name <- paste0("IRR[", i, "]")
  trace_data$irr[trace_data$year == (i-1)] <- jags5$BUGSoutput$sims.matrix[, param_name]
}

# Calculate summary statistics for each year
summary_data <- trace_data %>%
  group_by(year) %>%
  summarize(
    median_irr = median(irr),
    mean_irr = mean(irr),
    lower_95 = quantile(irr, 0.025),
    upper_95 = quantile(irr, 0.975)
  )


# Create the trace plot
trace_plot <- ggplot() +
  # Add individual traces (sample a subset to avoid overplotting)
  geom_line(data = trace_data %>% 
              filter(sample_id %% 20 == 0), # Sample every 20th MCMC sample
            aes(x = year, y = irr, group = sample_id),
            color = journal_colors[3], alpha = 0.09, size = 0.3) +
  # Add median line
  geom_line(data = summary_data,
            aes(x = year, y = median_irr),
            color = 'black', size = 1.2) +
  # Add reference line at IRR = 1
  geom_hline(yintercept = 1, linetype = "longdash", color = "gray30", size = 0.8) +
  annotate("text", x = 4.2, y = 1.02, label = "No effect (IRR = 1)", 
           hjust = 0, vjust = 0, size = 3, color = "black", fontface = "italic") +
  # Labels and scales
  labs(
    title = "Effect of ERPO Laws Over Time",
    subtitle = "Individual MCMC traces with median effect highlighted",
    x = "Years Since Implementation",
    y = "Incidence Rate Ratio (IRR)",
  ) +
  scale_x_continuous(breaks = -1:10, minor_breaks = NULL) +
  scale_y_continuous(breaks = seq(0.7, 1.3, 0.1), 
                     minor_breaks = seq(0.65, 1.35, 0.05),
                     limits = c(0.65, 1.15),
                     expand = c(0.02, 0)) +
  # Apply your custom theme
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, color = "gray30", margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = "gray40", hjust = 0, margin = margin(t = 10)),
    # Axis elements
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.text = element_text(size = 10, color = "gray20"),
    axis.ticks = element_line(color = "gray70", size = 0.3),
    # Grid lines
    panel.grid.major = element_line(color = "gray92", size = 0.4),
    panel.grid.minor = element_line(color = "gray96", size = 0.2),
    # Overall plot margins
    plot.margin = margin(20, 25, 20, 20)
  )

# 2. DENSITY PLOT FOR SELECTED YEARS
# ----------------------------------
# Select specific years to show distributions (years 1, 2, 3, 4, 5)
selected_years <- c(1, 2, 3, 4, 5)
density_data <- trace_data %>%
  filter(year %in% selected_years) %>%
  mutate(Year = factor(paste("Year", year)))

# Create the density plot
density_plot <- ggplot(density_data, aes(x = irr, fill = Year, color = Year)) +
  geom_rug(alpha = 0.7, show.legend = FALSE, size = 0.4) +
  geom_density(alpha = 0.65, size = 0.7) +
  geom_vline(xintercept = 1, linetype = "longdash", color = "gray30", size = 0.8) +
  annotate("text", x = 1.02, y = 3, label = "No effect (IRR = 1)", 
           hjust = 0, vjust = 0, size = 3, color = "black", angle = 90, fontface = "italic") +
  labs(
    title = "Posterior Distributions of ERPO Effect by Year",
    subtitle = "Incidence Rate Ratio distributions for years 1-5 after implementation",
    x = "Incidence Rate Ratio (IRR)",
    y = "Density",
    caption = "Note: Values < 1 indicate a reduction in incidents; values > 1 indicate an increase."
  ) +
  scale_x_continuous(breaks = seq(0.7, 1.3, 0.1), 
                     minor_breaks = seq(0.65, 1.35, 0.05),
                     limits = c(0.65, 1.1),
                     expand = c(0.02, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = journal_colors[1:5], name = "") +
  scale_color_manual(values = journal_colors[1:5], name = "") +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, color = "gray30", margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = "gray40", hjust = 0, margin = margin(t = 10)),
    # Axis elements
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.text = element_text(size = 10, color = "gray20"),
    axis.ticks = element_line(color = "gray70", size = 0.3),
    # Grid lines
    panel.grid.major = element_line(color = "gray92", size = 0.4),
    panel.grid.minor = element_line(color = "gray96", size = 0.2),
    # Change legend position to right and adjust spacing
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1.2, "lines"),
    legend.margin = margin(t = 5, r = 0, b = 5, l = 5),
    legend.box.margin = margin(l = 10),
    # Overall plot margins
    plot.margin = margin(20, 25, 20, 20)
  ) +
  # Configure legend to display vertically
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.8, size = 4), ncol = 1),
    color = FALSE  # Hide color legend (redundant with fill)
  )

# Print individual plots or combine them
# To view individual plots:
# print(trace_plot)
# print(density_plot)

# To combine both plots vertically
combined_plot <- trace_plot / density_plot
print(combined_plot)

# Save the combined plot
ggsave("erpo_effects_combined.png", combined_plot, width = 12, height = 10, dpi = 300)