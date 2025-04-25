# Load necessary libraries
library(tidyverse)
library(R2jags)
library(ggplot2)
library(gridExtra)
library(patchwork)  # For combining plots

# Load data
load("outputs/data/erpo_analysis_data.RData")

set.seed(123)

# -----------------------------------------------------
# 1. Data exploration and visualization
# -----------------------------------------------------

trend_data <- analysis_data %>%
  mutate(has_erpo = state %in% unique(analysis_data$state[analysis_data$erpo > 0]),
         erpo_status = if_else(has_erpo, "ERPO States", "Non-ERPO States")) %>%
  group_by(year, erpo_status) %>%
  summarize(
    total_deaths = sum(deaths),
    total_population = sum(population),
    avg_rate = total_deaths / total_population * 100000,
    .groups = "drop"
  )

ggplot(trend_data, aes(x = year, y = avg_rate, color = erpo_status)) +
  geom_line(size = 1) +
  geom_point() +
  labs(
    title = "Trends in Firearm Suicide Rates by ERPO Status",
    subtitle = "Average rates per 100,000 population",
    x = "Year",
    y = "Firearm Suicide Rate per 100,000",
    color = "State Group"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Create before-after plot for ERPO states
erpo_states <- analysis_data %>%
  group_by(state) %>%
  summarize(
    implemented_erpo = any(erpo > 0),
    min_erpo_year = if_else(implemented_erpo, min(year[erpo > 0]), NA_real_)
  ) %>%
  filter(implemented_erpo & !is.na(min_erpo_year))

# Create before-after data
before_after_data <- analysis_data %>%
  inner_join(erpo_states, by = "state") %>%
  mutate(
    period = case_when(
      year < min_erpo_year ~ "Before",
      year == min_erpo_year ~ "Implementation Year",
      year > min_erpo_year ~ "After"
    ),
    years_since = year - min_erpo_year
  )

# Create before-after plot
ggplot(before_after_data, 
       aes(x = years_since, y = crude_rate, group = state)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_line(alpha = 0.7) +
  facet_wrap(~ state, ncol = 4) +
  geom_point(data = before_after_data %>% filter(period == "Implementation Year"), 
             size = 2, shape = 18) +
  labs(
    title = "Firearm Suicide Rates Before and After ERPO Implementation",
    subtitle = "By state, centered on implementation year (0)",
    x = "Years Since ERPO Implementation",
    y = "Firearm Suicide Rate per 100,000"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# -----------------------------------------------------
# 2. Prepare data for JAGS models
# -----------------------------------------------------

# Create factors for states and years to use in models
state_factor <- as.numeric(factor(analysis_data$state))
year_factor <- as.numeric(factor(analysis_data$year))

# Create dummy variables for years (to use as fixed effects)
year_dummies <- model.matrix(~ factor(year) - 1, data = analysis_data)
colnames(year_dummies) <- paste0("year_", unique(analysis_data$year))

# -----------------------------------------------------
# 3. Define and run simple JAGS models
# -----------------------------------------------------

# Model 1: Simple model with binary ERPO indicator
jags_model1 <- "
model {
  # Likelihood
  for (i in 1:n) {
    log_rate[i] <- log(population[i]) + beta0 + beta_erpo * erpo[i] + year_effect[year[i]]
    deaths[i] ~ dpois(exp(log_rate[i]))
  }
  
  # Priors
  beta0 ~ dnorm(0, 0.001)
  beta_erpo ~ dnorm(0, 10)  # Weakly informative prior centered at 0
  
  # Year effects
  for (j in 1:n_years) {
    year_effect[j] ~ dnorm(0, 0.1)
  }
  
  # Derived quantities
  IRR_erpo <- exp(beta_erpo)
}
"

# Prepare data for model 1
jags_data1 <- list(
  n = nrow(analysis_data),
  n_years = length(unique(year_factor)),
  deaths = analysis_data$deaths,
  population = analysis_data$population,
  erpo = analysis_data$erpo,
  year = year_factor
)

# Parameters to monitor
params1 <- c("beta0", "beta_erpo", "IRR_erpo")

# Run JAGS model 1
cat("Running Model 1: Simple Poisson with ERPO binary indicator...\n")
jags1 <- jags(
  data = jags_data1,
  parameters.to.save = params1,
  model.file = textConnection(jags_model1),
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 1000,
  n.thin = 5
)

# Model 2: Model with phase-in effect
jags_model2 <- "
model {
  # Likelihood
  for (i in 1:n) {
    log_rate[i] <- log(population[i]) + beta0 + beta_phase * erpo_phase_in[i] + year_effect[year[i]]
    deaths[i] ~ dpois(exp(log_rate[i]))
  }
  
  # Priors
  beta0 ~ dnorm(0, 0.001)
  beta_phase ~ dnorm(0, 10)  # Weakly informative prior centered at 0
  
  # Year effects
  for (j in 1:n_years) {
    year_effect[j] ~ dnorm(0, 0.1)
  }
  
  # Derived quantities
  IRR_phase <- exp(beta_phase)
}
"

# Prepare data for model 2
jags_data2 <- list(
  n = nrow(analysis_data),
  n_years = length(unique(year_factor)),
  deaths = analysis_data$deaths,
  population = analysis_data$population,
  erpo_phase_in = analysis_data$erpo_phase_in,
  year = year_factor
)

# Parameters to monitor
params2 <- c("beta0", "beta_phase", "IRR_phase")

# Run JAGS model 2
cat("Running Model 2: Simple Poisson with ERPO phase-in...\n")
jags2 <- jags(
  data = jags_data2,
  parameters.to.save = params2,
  model.file = textConnection(jags_model2),
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 1000,
  n.thin = 5
)

# Model 3: Negative binomial model with state random effects and binary ERPO
jags_model3 <- "
model {
  # Likelihood
  for (i in 1:n) {
    log_mu[i] <- log(population[i]) + beta0 + beta_erpo * erpo[i] + year_effect[year[i]] + state_effect[state[i]]
    p[i] <- r / (r + exp(log_mu[i]))
    deaths[i] ~ dnegbin(p[i], r)
  }
  
  # Priors
  beta0 ~ dnorm(0, 0.001)
  beta_erpo ~ dnorm(0, 10)  # Weakly informative prior centered at 0
  r ~ dgamma(0.1, 0.1)      # Prior for the dispersion parameter
  
  # Random effects for states
  for (s in 1:n_states) {
    state_effect[s] ~ dnorm(0, tau_state)
  }
  
  # Year effects
  for (j in 1:n_years) {
    year_effect[j] ~ dnorm(0, 0.1)
  }
  
  # Hyperpriors
  tau_state ~ dgamma(0.1, 0.1)
  
  # Derived quantities
  IRR_erpo <- exp(beta_erpo)
}
"

# Prepare data for model 3
jags_data3 <- list(
  n = nrow(analysis_data),
  n_states = length(unique(state_factor)),
  n_years = length(unique(year_factor)),
  deaths = analysis_data$deaths,
  population = analysis_data$population,
  erpo = analysis_data$erpo,
  state = state_factor,
  year = year_factor
)

# Parameters to monitor
params3 <- c("beta0", "beta_erpo", "IRR_erpo", "r", "tau_state")

# Run JAGS model 3
cat("Running Model 3: Negative Binomial with state random effects and binary ERPO...\n")
jags3 <- jags(
  data = jags_data3,
  parameters.to.save = params3,
  model.file = textConnection(jags_model3),
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 1000,
  n.thin = 5
)

# Model 4: Negative binomial model with state random effects and phase-in effect
jags_model4 <- "
model {
  # Likelihood
  for (i in 1:n) {
    log_mu[i] <- log(population[i]) + beta0 + beta_phase * erpo_phase_in[i] + year_effect[year[i]] + state_effect[state[i]]
    p[i] <- r / (r + exp(log_mu[i]))
    deaths[i] ~ dnegbin(p[i], r)
  }
  
  # Priors
  beta0 ~ dnorm(0, 0.001)
  beta_phase ~ dnorm(0, 10)  # Weakly informative prior centered at 0
  r ~ dgamma(0.1, 0.1)       # Prior for the dispersion parameter
  
  # Random effects for states
  for (s in 1:n_states) {
    state_effect[s] ~ dnorm(0, tau_state)
  }
  
  # Year effects
  for (j in 1:n_years) {
    year_effect[j] ~ dnorm(0, 0.1)
  }
  
  # Hyperpriors
  tau_state ~ dgamma(0.1, 0.1)
  
  # Derived quantities
  IRR_phase <- exp(beta_phase)
}
"

# Prepare data for model 4
jags_data4 <- list(
  n = nrow(analysis_data),
  n_states = length(unique(state_factor)),
  n_years = length(unique(year_factor)),
  deaths = analysis_data$deaths,
  population = analysis_data$population,
  erpo_phase_in = analysis_data$erpo_phase_in,
  state = state_factor,
  year = year_factor
)

# Parameters to monitor
params4 <- c("beta0", "beta_phase", "IRR_phase", "r", "tau_state")

# Run JAGS model 4
cat("Running Model 4: Negative Binomial with state random effects and phase-in effect...\n")
jags4 <- jags(
  data = jags_data4,
  parameters.to.save = params4,
  model.file = textConnection(jags_model4),
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 1000,
  n.thin = 5
)

# 4. Extract and visualize results
# -----------------------------------------------------

# Function to extract posterior samples and calculate probabilities
extract_jags_results <- function(jags_fit, param_name) {
  # Extract posterior samples
  posterior <- jags_fit$BUGSoutput$sims.matrix
  
  # Get IRR samples
  irr_samples <- posterior[, param_name]
  
  # Calculate summary statistics
  median_irr <- median(irr_samples)
  mean_irr <- mean(irr_samples)
  sd_irr <- sd(irr_samples)
  lower_95 <- quantile(irr_samples, 0.025)
  upper_95 <- quantile(irr_samples, 0.975)
  lower_80 <- quantile(irr_samples, 0.1)
  upper_80 <- quantile(irr_samples, 0.9)
  
  # Calculate probability of reduction
  prob_reduction <- mean(irr_samples < 1)
  
  return(list(
    irr_samples = irr_samples,
    median = median_irr,
    mean = mean_irr,
    sd = sd_irr,
    lower_95 = lower_95,
    upper_95 = upper_95,
    lower_80 = lower_80,
    upper_80 = upper_80,
    prob_reduction = prob_reduction
  ))
}

# Extract results
results1 <- extract_jags_results(jags1, "IRR_erpo")
results2 <- extract_jags_results(jags2, "IRR_phase")
results3 <- extract_jags_results(jags3, "IRR_erpo")
results4 <- extract_jags_results(jags4, "IRR_phase")

# Create results table
results_table <- data.frame(
  Model = c("Simple Poisson (Binary)", "Simple Poisson (Phase-in)", 
            "Negative Binomial (Binary)", "Negative Binomial (Phase-in)"),
  Median_IRR = c(results1$median, results2$median, results3$median, results4$median),
  Mean_IRR = c(results1$mean, results2$mean, results3$mean, results4$mean),
  SD = c(results1$sd, results2$sd, results3$sd, results4$sd),
  Lower_95CI = c(results1$lower_95, results2$lower_95, results3$lower_95, results4$lower_95),
  Upper_95CI = c(results1$upper_95, results2$upper_95, results3$upper_95, results4$upper_95),
  Lower_80CI = c(results1$lower_80, results2$lower_80, results3$lower_80, results4$lower_80),
  Upper_80CI = c(results1$upper_80, results2$upper_80, results3$upper_80, results4$upper_80),
  Prob_Reduction = c(results1$prob_reduction, results2$prob_reduction, 
                     results3$prob_reduction, results4$prob_reduction)
)

# Round numeric columns
results_table <- results_table %>%
  mutate(across(where(is.numeric), round, 3))

# Create a forest plot of model results
results_plot_data <- results_table %>%
  mutate(CI_text = paste0(Median_IRR, " (", Lower_95CI, "-", Upper_95CI, ")"))

ggplot(results_plot_data, 
                      aes(y = Model, x = Median_IRR, xmin = Lower_95CI, xmax = Upper_95CI)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray") +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2) +
  geom_text(aes(x = Upper_95CI + 0.03, label = paste0(round(Prob_Reduction * 100, 1), "%")), 
            hjust = 0, size = 3) +
  labs(
    title = "Effect of ERPO Laws on Firearm Suicide Rates",
    subtitle = "Incidence Rate Ratios with 95% Credible Intervals",
    x = "Incidence Rate Ratio (IRR)",
    y = NULL,
    caption = "Percentages indicate posterior probability of a reduction in suicide rates"
  ) +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0))

# Create density plots of the posterior distributions
density_plot_data <- data.frame(
  Model = rep(c("Simple Poisson (Binary)", "Simple Poisson (Phase-in)", 
                "Negative Binomial (Binary)", "Negative Binomial (Phase-in)"), 
              c(length(results1$irr_samples), length(results2$irr_samples), 
                length(results3$irr_samples), length(results4$irr_samples))),
  IRR = c(results1$irr_samples, results2$irr_samples, 
          results3$irr_samples, results4$irr_samples)
)

journal_colors <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00")

# Define model order
model_order <- c("Simple Poisson (Binary)", 
                 "Simple Poisson (Phase-in)", 
                 "Negative Binomial (Binary)", 
                 "Negative Binomial (Phase-in)")

# Convert Model to factor with specific order
density_plot_data$Model <- factor(density_plot_data$Model, levels = model_order)

ggplot(density_plot_data, aes(x = IRR, fill = Model, color = Model)) +
  geom_rug(alpha = 0.7, show.legend = FALSE, size = 0.4) +
  geom_density(alpha = 0.65, size = 0.7) +
  geom_vline(xintercept = 1, linetype = "longdash", color = "gray30", size = 0.8) +
  annotate("text", x = 1.02, y = 25, label = "No effect (IRR = 1)", 
           hjust = 0, vjust = 0, size = 3, color = "black", angle = 90, fontface = "italic") +
  labs(
    title = "Posterior Distributions of ERPO Effect",
    subtitle = "Incidence Rate Ratios across different models",
    x = "Incidence Rate Ratio (IRR)",
    y = "Density",
    caption = "Note: Values < 1 indicate a reduction in incidents; values > 1 indicate an increase."
  ) +
  scale_x_continuous(breaks = seq(0.7, 1.3, 0.1), 
                     minor_breaks = seq(0.65, 1.35, 0.05),
                     limits = c(0.65, 1.1),
                     expand = c(0.02, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = journal_colors, name = "Model Type") +
  scale_color_manual(values = journal_colors, name = "Model Type") +
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

# Create trace plots for key parameters
trace_plot1 <- traceplot(jags3, parameters = c("beta_erpo"))
trace_plot2 <- traceplot(jags4, parameters = c("beta_phase"))

# -----------------------------------------------------
# 5. Create an evidence plot with existing literature

# -----------------------------------------------------
# 7. Print results and summary
# -----------------------------------------------------

cat("\n=== Effects of ERPO Laws on Firearm Suicide Rates ===\n")
cat("------------------------------------------------------\n")
cat("Results from different Bayesian models (JAGS):\n\n")
print(results_table)

cat("\nDIC (Deviance Information Criterion) for model comparison:\n")
dic_comparison <- data.frame(
  Model = c("Simple Poisson (Binary)", "Simple Poisson (Phase-in)", 
            "Negative Binomial (Binary)", "Negative Binomial (Phase-in)"),
  DIC = c(jags1$BUGSoutput$DIC, jags2$BUGSoutput$DIC, 
          jags3$BUGSoutput$DIC, jags4$BUGSoutput$DIC)
)
print(dic_comparison)


# Save all results
results <- list(
  jags_models = list(jags1 = jags1, jags2 = jags2, jags3 = jags3, jags4 = jags4),
  results_table = results_table,
  dic_comparison = dic_comparison,
  posterior_samples = list(
    model1 = results1,
    model2 = results2,
    model3 = results3,
    model4 = results4
  ),
  trend_data = trend_data,
  before_after_data = before_after_data
)

save(results, file = "outputs/erpo_jags_results.RData")















# Fifth JAGS Model with Autoregressive Structure
# This script adds a more complex autoregressive model similar to the Stan version

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
time_effects_plot <- ggplot(time_effects, aes(x = Year, y = IRR, ymin = Lower_95CI, ymax = Upper_95CI)) +
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
ggsave("manuscript/figs/time_effects_plot.pdf", time_effects_plot, width = 8, height = 6)

# Save full model 5 results
results5_full <- list(
  jags_model = jags5,
  summary = results5_summary,
  time_effects = time_effects
)

# Return the results for potential further analysis
results5_full











# Difference-in-Differences and Event Study Analysis for ERPO Laws
# This script adds DiD and event study analyses to evaluate ERPO laws

# Load necessary libraries
library(tidyverse)
library(R2jags)
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
# 1. Prepare data for DiD and Event Study Analysis
# -----------------------------------------------------

# Create relative time variable
event_data <- analysis_data %>%
  group_by(state) %>%
  mutate(
    implementation_year = if_else(any(erpo > 0), min(year[erpo > 0]), NA_real_),
    rel_year = year - implementation_year,
    treated = !is.na(implementation_year)
  ) %>%
  ungroup()

# Keep relative years within a reasonable range (-5 to +5)
event_data <- event_data %>%
  mutate(
    rel_year_capped = case_when(
      is.na(rel_year) ~ NA_real_,
      rel_year < -5 ~ -5,
      rel_year > 5 ~ 5,
      TRUE ~ rel_year
    )
  )

# -----------------------------------------------------
# 2. Bayesian Difference-in-Differences Analysis
# -----------------------------------------------------

# Define the JAGS model for DiD
jags_did_model <- "
model {
  # Likelihood
  for (i in 1:n) {
    log_mu[i] <- log(population[i]) + beta0 + 
                 alpha_state[state[i]] + 
                 alpha_year[year[i]] + 
                 beta_did * treated[i] * post[i]
    
    p[i] <- r / (r + exp(log_mu[i]))
    deaths[i] ~ dnegbin(p[i], r)
  }
  
  # Priors
  beta0 ~ dnorm(0, 0.001)
  beta_did ~ dnorm(0, 10)     # DiD effect
  r ~ dgamma(0.1, 0.1)        # Dispersion parameter
  
  # State fixed effects
  for (s in 1:n_states) {
    alpha_state[s] ~ dnorm(0, tau_state)
  }
  
  # Year fixed effects
  for (t in 1:n_years) {
    alpha_year[t] ~ dnorm(0, tau_year)
  }
  
  # Hyperpriors
  tau_state ~ dgamma(0.1, 0.1)
  tau_year ~ dgamma(0.1, 0.1)
  
  # Derived quantities
  IRR_did <- exp(beta_did)
}
"

# Create state and year factors
state_factor <- as.numeric(factor(analysis_data$state))
year_factor <- as.numeric(factor(analysis_data$year))

# Create treatment and post indicators for DiD
did_data <- analysis_data %>%
  group_by(state) %>%
  mutate(
    treated = as.numeric(any(erpo > 0)),
    implementation_year = if_else(treated == 1, min(year[erpo > 0]), NA_real_),
    post = as.numeric(treated == 1 & year >= implementation_year)
  ) %>%
  ungroup()

# Prepare data for JAGS DiD model
jags_did_data <- list(
  n = nrow(did_data),
  n_states = length(unique(state_factor)),
  n_years = length(unique(year_factor)),
  deaths = did_data$deaths,
  population = did_data$population,
  state = as.numeric(factor(did_data$state)),
  year = as.numeric(factor(did_data$year)),
  treated = did_data$treated,
  post = did_data$post
)

# Parameters to monitor for DiD
did_params <- c("beta0", "beta_did", "IRR_did", "r")

# Run JAGS DiD model
cat("Running Bayesian Difference-in-Differences model...\n")
jags_did <- jags(
  data = jags_did_data,
  parameters.to.save = did_params,
  model.file = textConnection(jags_did_model),
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 2000,
  n.thin = 5
)

# Extract DiD results
did_results <- list(
  beta_did = jags_did$BUGSoutput$sims.matrix[, "beta_did"],
  IRR_did = jags_did$BUGSoutput$sims.matrix[, "IRR_did"],
  r = jags_did$BUGSoutput$sims.matrix[, "r"]
)

did_summary <- list(
  median_irr = median(did_results$IRR_did),
  mean_irr = mean(did_results$IRR_did),
  sd_irr = sd(did_results$IRR_did),
  lower_95 = quantile(did_results$IRR_did, 0.025),
  upper_95 = quantile(did_results$IRR_did, 0.975),
  prob_reduction = mean(did_results$IRR_did < 1)
)

# Print DiD results
cat("\n=== Bayesian Difference-in-Differences Results ===\n")
cat(sprintf("IRR: %.3f (95%% CI: %.3f to %.3f)\n", 
            did_summary$median_irr, 
            did_summary$lower_95, 
            did_summary$upper_95))
cat(sprintf("Probability of reduction: %.1f%%\n", did_summary$prob_reduction * 100))
cat(sprintf("DIC: %.2f\n\n", jags_did$BUGSoutput$DIC))

# -----------------------------------------------------
# 3. Bayesian Event Study Analysis
# -----------------------------------------------------

# Define the JAGS model for Event Study
jags_event_model <- "
model {
  # Likelihood
  for (i in 1:n) {
    # Only include effect for states that implemented ERPO
    effect_term[i] <- equals(treated[i], 1) * 
                      inprod(rel_year_dummy[i,], beta_rel)
    
    log_mu[i] <- log(population[i]) + beta0 + 
                 alpha_state[state[i]] + 
                 alpha_year[year[i]] + 
                 effect_term[i]
    
    p[i] <- r / (r + exp(log_mu[i]))
    deaths[i] ~ dnegbin(p[i], r)
  }
  
  # Priors
  beta0 ~ dnorm(0, 0.001)
  r ~ dgamma(0.1, 0.1)      # Dispersion parameter
  
  # Relative year effects (event study coefficients)
  # Year 0 (implementation year) is the reference
  beta_rel[1] ~ dnorm(0, 10)  # Year -5 or earlier
  for (t in 2:5) {            # Years -4 to -1
    beta_rel[t] ~ dnorm(0, 10)
  }
  # beta_rel[6] = 0            # Year 0 (implementation) is reference
  for (t in 6:10) {           # Years 1 to 5
    beta_rel[t] ~ dnorm(0, 10)
  }
  beta_rel[6] <- 0            # Year 0 as reference
  
  # State fixed effects
  for (s in 1:n_states) {
    alpha_state[s] ~ dnorm(0, tau_state)
  }
  
  # Year fixed effects
  for (t in 1:n_years) {
    alpha_year[t] ~ dnorm(0, tau_year)
  }
  
  # Hyperpriors
  tau_state ~ dgamma(0.1, 0.1)
  tau_year ~ dgamma(0.1, 0.1)
  
  # Derived quantities - IRRs for each relative year
  for (t in 1:10) {
    if (t != 6) {  # Skip reference year
      IRR_rel[t] <- exp(beta_rel[t])
    } else {
      IRR_rel[t] <- 1  # Reference year has IRR = 1
    }
  }
}
"

# Create relative year dummies for event study
# We'll create dummies for years -5 to 5 (with 0 as reference)
rel_years <- -5:5
rel_year_dummies <- matrix(0, nrow = nrow(event_data), ncol = length(rel_years))

for (i in 1:length(rel_years)) {
  if (rel_years[i] != 0) {  # Skip reference year
    rel_year_dummies[, i] <- as.numeric(event_data$rel_year_capped == rel_years[i])
  }
}

# Prepare data for JAGS Event Study model
jags_event_data <- list(
  n = nrow(event_data),
  n_states = length(unique(state_factor)),
  n_years = length(unique(year_factor)),
  deaths = event_data$deaths,
  population = event_data$population,
  state = as.numeric(factor(event_data$state)),
  year = as.numeric(factor(event_data$year)),
  treated = as.numeric(!is.na(event_data$implementation_year)),
  rel_year_dummy = rel_year_dummies
)

# Parameters to monitor for Event Study
event_params <- c("beta_rel", "IRR_rel", "r")

# Run JAGS Event Study model
cat("Running Bayesian Event Study model...\n")
jags_event <- jags(
  data = jags_event_data,
  parameters.to.save = event_params,
  model.file = textConnection(jags_event_model),
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 2000,
  n.thin = 5
)

# Extract event study results
event_results <- list()
rel_irr <- matrix(NA, nrow = jags_event$BUGSoutput$n.sims, ncol = length(rel_years))

for (i in 1:length(rel_years)) {
  if (rel_years[i] != 0) {  # Skip reference year (0)
    param_name <- paste0("IRR_rel[", i, "]")
    rel_irr[, i] <- jags_event$BUGSoutput$sims.matrix[, param_name]
  } else {
    # Reference year has IRR = 1
    rel_irr[, i] <- 1
  }
}

# Create summary data frame for event study results
event_summary <- data.frame(
  rel_year = rel_years,
  irr = apply(rel_irr, 2, median),
  irr_mean = apply(rel_irr, 2, mean),
  lower_95 = apply(rel_irr, 2, quantile, 0.025),
  upper_95 = apply(rel_irr, 2, quantile, 0.975)
)

# Set IRR and CI to exactly 1 for reference year
event_summary$irr[event_summary$rel_year == 0] <- 1
event_summary$lower_95[event_summary$rel_year == 0] <- 1
event_summary$upper_95[event_summary$rel_year == 0] <- 1

# Create event study plot
event_plot <- ggplot(event_summary, 
                     aes(x = rel_year, y = irr, ymin = lower_95, ymax = upper_95)) +
  geom_point(size = 3) +
  geom_errorbar(width = 0.2) +
  geom_line(linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "solid", color = "darkgray") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "darkgray") +
  labs(
    title = "Event Study: Effect of ERPO Laws on Firearm Suicide Rates",
    subtitle = "Relative to implementation year (year 0)",
    x = "Years Relative to ERPO Implementation",
    y = "Incidence Rate Ratio (IRR)"
  ) +
  scale_x_continuous(breaks = -5:5) +
  theme_minimal()

# Save event study plot
ggsave("manuscript/figs/event_study_plot.pdf", event_plot, width = 10, height = 6)

# Print event study results
cat("\n=== Bayesian Event Study Results ===\n")
cat("Relative Year IRRs (compared to implementation year):\n")
print(event_summary)
cat(sprintf("DIC: %.2f\n\n", jags_event$BUGSoutput$DIC))

# -----------------------------------------------------
# 4. Combined visualization of all models
# -----------------------------------------------------

# Create a combined data frame of all model results
# Assuming results_table exists with the first 4 models
# and model5_row from the previous script

combined_results <- rbind(
  results_table,  # Original 4 models
  model5_row,     # Model 5 (Autoregressive)
  data.frame(     # DiD model
    Model = "Difference-in-Differences",
    Median_IRR = round(did_summary$median_irr, 3),
    Mean_IRR = round(did_summary$mean_irr, 3),
    SD = round(did_summary$sd_irr, 3),
    Lower_95CI = round(did_summary$lower_95, 3),
    Upper_95CI = round(did_summary$upper_95, 3),
    Lower_80CI = NA,  # Not calculated for DiD
    Upper_80CI = NA,  # Not calculated for DiD
    Prob_Reduction = round(did_summary$prob_reduction, 3)
  )
)

# Create an updated forest plot with all models
forest_plot_combined <- ggplot(combined_results, 
                               aes(y = reorder(Model, Median_IRR), 
                                   x = Median_IRR, 
                                   xmin = Lower_95CI, 
                                   xmax = Upper_95CI)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray") +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2) +
  geom_text(aes(x = Upper_95CI + 0.03, 
                label = paste0(round(Prob_Reduction * 100, 1), "%")), 
            hjust = 0, size = 3) +
  labs(
    title = "Effect of ERPO Laws on Firearm Suicide Rates",
    subtitle = "Comparing Different Modeling Approaches",
    x = "Incidence Rate Ratio (IRR)",
    y = NULL,
    caption = "Percentages indicate posterior probability of a reduction in suicide rates"
  ) +
  scale_x_continuous(breaks = seq(0.6, 1.2, 0.1)) +
  coord_cartesian(xlim = c(0.6, 1.2)) +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0))

# Save the combined forest plot
ggsave("manuscript/figs/forest_plot_combined.pdf", forest_plot_combined, width = 10, height = 8)

# Combine DiD, Event Study, and Model 5 results for saving
did_event_results <- list(
  jags_did = jags_did,
  did_summary = did_summary,
  jags_event = jags_event,
  event_summary = event_summary,
  forest_plot_combined = forest_plot_combined
)

# Save results
save(did_event_results, file = "manuscript/src/erpo_did_event_results.RData")

# Return results for potential further use
did_event_results