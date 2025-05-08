library(readxl)
library(tidyverse)
library(dplyr)
library(tidyr)
library(lubridate)

NFA_1999_2020 <- read_excel("raw/placebo/NFA_1999_2020.xlsx") %>%
  select(-Notes) %>%
  filter(!is.na(State)) %>%
  filter(Year < 2018)

NFA_1999_2020$Deaths <- as.numeric(NFA_1999_2020$Deaths)
NFA_1999_2020$Deaths[is.na(NFA_1999_2020$Deaths)] <- 5

NFA_1999_2020$Population <- as.numeric(NFA_1999_2020$Population)
NFA_1999_2020$`Crude Rate` <- NFA_1999_2020$Deaths/NFA_1999_2020$Population * 100000

NFA_2018_2023 <- read_excel("raw/placebo/NFA_2018_2023.xlsx") %>%
  select(-Notes)  %>%
  filter(!is.na(State))

NFA_2018_2023$Deaths <- as.numeric(NFA_2018_2023$Deaths)
NFA_2018_2023$Deaths[is.na(NFA_2018_2023$Deaths)] <- 5

NFA_2018_2023$Population <- as.numeric(NFA_2018_2023$Population)
NFA_2018_2023$`Crude Rate` <- NFA_2018_2023$Deaths/NFA_2018_2023$Population * 100000

NFA_1999_2023 <- rbind(NFA_1999_2020, NFA_2018_2023)

NFA_1999_2023 <- NFA_1999_2023 %>%
  rename_with(~ gsub(" ", "_", .x)) %>%
  rename_with(~ tolower(.x)) %>%
  select(-state_code, -year_code) 


#-------------------------------------------------------------------------------

# Create a data frame for ERPO law implementation
erpo_laws <- data.frame(
  state = c("California", "Colorado", "Connecticut", "Delaware", "District of Columbia", 
            "Florida", "Hawaii", "Illinois", "Indiana", "Maryland", 
            "Massachusetts", "Michigan", "Minnesota", "Nevada", "New Jersey", 
            "New Mexico", "New York", "Oregon", "Rhode Island", "Vermont", 
            "Virginia", "Washington"),
  implementation_date = c("2016-01-01", "2019-04-12", "1999-10-01", "2018-12-27", "2019-01-30",
                          "2018-03-09", "2020-01-01", "2019-01-01", "2005-07-01", "2018-10-01",
                          "2018-08-17", "2024-01-01", "2024-01-01", "2020-01-01", "2019-09-01",
                          "2020-05-20", "2019-08-24", "2018-01-01", "2018-06-01", "2018-04-11",
                          "2020-07-01", "2016-12-08"),
  law_enforcement = c(TRUE, TRUE, TRUE, TRUE, TRUE, 
                      TRUE, TRUE, TRUE, TRUE, TRUE, 
                      TRUE, TRUE, TRUE, TRUE, TRUE, 
                      TRUE, TRUE, TRUE, TRUE, TRUE, 
                      TRUE, TRUE),
  family_member = c(TRUE, TRUE, TRUE, TRUE, TRUE, 
                    FALSE, TRUE, TRUE, FALSE, TRUE, 
                    TRUE, TRUE, TRUE, TRUE, TRUE, 
                    FALSE, TRUE, TRUE, FALSE, TRUE, 
                    TRUE, TRUE),
  health_professional = c(FALSE, TRUE, TRUE, FALSE, TRUE, 
                          FALSE, TRUE, FALSE, FALSE, TRUE, 
                          TRUE, TRUE, FALSE, FALSE, FALSE, 
                          FALSE, TRUE, FALSE, FALSE, FALSE, 
                          FALSE, FALSE),
  educator = c(TRUE, TRUE, FALSE, FALSE, FALSE, 
               FALSE, TRUE, FALSE, FALSE, FALSE, 
               TRUE, FALSE, FALSE, FALSE, FALSE, 
               FALSE, TRUE, FALSE, FALSE, FALSE, 
               FALSE, FALSE),
  max_duration = c(5, 1, NA, 1, 1, 
                   1, 1, 1, NA, 1, 
                   1, 1, 1, 1, NA, 
                   1, 1, 1, 1, 0.5, 
                   0.5, 1)  # Duration in years, NA for "until terminated"
)

# Convert implementation_date to Date format
erpo_laws$implementation_date <- as.Date(erpo_laws$implementation_date)

generate_erpo_panel <- function(
    erpo_data,
    start_year = 1999,
    end_year   = 2023,
    all_states = unique(NFA_1999_2023$state)
) {
  # 1) build all state-year combos
  panel <- expand.grid(
    state = all_states,
    year  = seq(start_year, end_year),
    stringsAsFactors = FALSE
  ) %>%
    arrange(state, year) %>%
    # initialize all ERPO flags to 0
    mutate(
      erpo                   = 0L,
      erpo_law_enforcement   = 0L,
      erpo_family            = 0L,
      erpo_health            = 0L,
      erpo_educator          = 0L
    )
  
  # 2) loop through each law, flip to 1 for year >= impl_year
  for (i in seq_len(nrow(erpo_data))) {
    st   <- erpo_data$state[i]
    impl_year <- year(erpo_data$implementation_date[i])
    
    # simple binary turn-on
    panel <- panel %>%
      mutate(
        erpo = if_else(state == st & year >= impl_year, 1L, erpo),
        erpo_law_enforcement = if_else(state == st & year >= impl_year, 
                                       as.integer(erpo_data$law_enforcement[i]), 
                                       erpo_law_enforcement),
        erpo_family          = if_else(state == st & year >= impl_year, 
                                       as.integer(erpo_data$family_member[i]), 
                                       erpo_family),
        erpo_health          = if_else(state == st & year >= impl_year, 
                                       as.integer(erpo_data$health_professional[i]), 
                                       erpo_health),
        erpo_educator        = if_else(state == st & year >= impl_year, 
                                       as.integer(erpo_data$educator[i]), 
                                       erpo_educator)
      )
  }
  
  return(panel)
}


erpo_panel <- generate_erpo_panel(erpo_laws)

# Merge with firearm suicide data
# First aggregate the suicide data to total (both sexes)
NFA_total <- NFA_1999_2023 
# Merge the data
merged_data <- left_join(erpo_panel, NFA_total, by = c("state", "year"))

# Create lagged variables for autoregressive terms
merged_data <- merged_data %>%
  arrange(state, year) %>%
  group_by(state) %>%
  mutate(
    deaths_lag1 = lag(deaths, 1),
    deaths_lag2 = lag(deaths, 2),
    population_lag1 = lag(population, 1),
    population_lag2 = lag(population, 2),
    log_rate_lag1 = log(deaths_lag1 / (population_lag1/1e5)),
    log_rate_lag2 = log(deaths_lag2 / (population_lag2/1e5))
  )

# Create the phase-in indicator (increasing from 0 to 1 over 5 years after implementation)
merged_data <- merged_data %>%
  group_by(state) %>%
  arrange(state, year) %>%
  mutate(
    # Find the first year with any ERPO implementation
    first_erpo_year = ifelse(any(erpo > 0), min(year[erpo > 0]), Inf),
    # Create phase-in variable that increases from 0 to 1 over 5 years
    erpo_phase_in = case_when(
      year < first_erpo_year ~ 0,
      year >= first_erpo_year + 5 ~ 1,
      TRUE ~ (year - first_erpo_year + first(erpo[year == first_erpo_year])) / 5
    )
  ) %>%
  select(-first_erpo_year) %>%
  ungroup()

# Handle missing values in log rates (for zeros)
merged_data$log_rate_lag1[is.infinite(merged_data$log_rate_lag1)] <- log(0.5/(merged_data$population_lag1[is.infinite(merged_data$log_rate_lag1)]/1e5))
merged_data$log_rate_lag2[is.infinite(merged_data$log_rate_lag2)] <- log(0.5/(merged_data$population_lag2[is.infinite(merged_data$log_rate_lag2)]/1e5))

# Drop the first two years which have missing lagged variables
analysis_data <- merged_data %>%
  filter(!is.na(log_rate_lag2),
         state != 'Connecticut', state != 'District of Columbia')

source("src/stateyeardata.R")
data
analysis_data
analysis_data <- merge(
  analysis_data,
  data,
  by = c("state", "year")
)

# Save the prepared data



# 0) Load libraries
library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bayesplot)
library(loo)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


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

# 5) Build policy matrices P, P1, P2
policy_vars <- c("erpo","erpo_phase_in")
P  <- as.matrix(analysis_data[, policy_vars])
P1 <- as.matrix(analysis_data %>%
                  group_by(state) %>%
                  arrange(state, year) %>%
                  mutate(across(all_of(policy_vars), ~lag(.,1,0))) %>%
                  ungroup() %>%
                  select(all_of(policy_vars)))
P2 <- as.matrix(analysis_data %>%
                  group_by(state) %>%
                  arrange(state, year) %>%
                  mutate(across(all_of(policy_vars), ~lag(.,2,0))) %>%
                  ungroup() %>%
                  select(all_of(policy_vars)))

# 6) Prepare covariates (without LASSO shrinkage)
# Combine all covariates into a single matrix
all_covariates <- c("HFR", "under_100_percent", "x100_199_percent", "x200_399_percent", 
                    "x400_percent", "white", "black", "hispanic", "asian",
                    "adults_19_25", "adults_26_34", "adults_35_54", "adults_55_64", "x65")
X <- as.matrix(analysis_data[, all_covariates])

# 7) Assemble Stan data
placebo_stan_data <- list(
  N = n_obs,
  n_states = length(unique(state_factor)),
  n_years = n_years,
  KP = ncol(P),
  Kx = ncol(X),
  
  y = as.integer(analysis_data$deaths),
  l1 = as.numeric(analysis_data$log_rate_lag1),
  l2 = as.numeric(analysis_data$log_rate_lag2),
  offset = as.numeric(log(analysis_data$population)),
  offset1 = as.numeric(log(analysis_data$population_lag1)),
  offset2 = as.numeric(log(analysis_data$population_lag2)),
  
  P = P,
  P1 = P1,
  P2 = P2,
  
  X = X,
  
  T = year_dummies,
  T1 = year_dummies_lag1,
  T2 = year_dummies_lag2,
  
  S = state_factor,
  prior_sd = 0.064
)

# 8) Write Stan model to file
placebo_stan_code <- '
// ERPO_model.stan
// Based on Schell et al.\'s approach but without LASSO
data {
  int<lower=1> N;                  // number of observations
  int<lower=1> n_states;           // number of states
  int<lower=1> n_years;            // number of year effects
  int<lower=1> KP;                 // number of policy variables
  int<lower=1> Kx;                 // number of covariates
  
  int<lower=0> y[N];               // outcome count (deaths)
  vector[N] l1;                    // lagged outcome, lag 1 (log scale)
  vector[N] l2;                    // lagged outcome, lag 2 (log scale)
  
  matrix[N, KP] P;                 // policy matrix
  matrix[N, KP] P1;                // lagged policy matrix, lag 1
  matrix[N, KP] P2;                // lagged policy matrix, lag 2
  
  vector[N] offset;                // offset (log population)
  vector[N] offset1;               // lagged offset, lag 1
  vector[N] offset2;               // lagged offset, lag 2
  
  matrix[N, Kx] X;                 // all covariates without shrinkage
  
  matrix[N, n_years] T;            // time effects matrix
  matrix[N, n_years] T1;           // time effects matrix, lag 1
  matrix[N, n_years] T2;           // time effects matrix, lag 2
  
  int<lower=1, upper=n_states> S[N];  // state index for each observation
  
  real<lower=0> prior_sd;          // prior SD for policy effects
}

parameters {
  real alpha;                      // intercept
  real<lower=0, upper=1> delta1;   // AR(1) coefficient (constrained 0-1)
  real<lower=-1, upper=1> delta2;  // AR(2) coefficient (constrained -1 to 1)
  
  vector[KP] beta;                 // policy effects
  vector[Kx] gamma_X;              // covariate effects
  
  vector[n_years] zeta;            // year effects
  
  real<lower=0> invphi;            // negative binomial dispersion
  
  vector[n_states] state_raw;      // raw state effects
}

transformed parameters {
  vector[N] log_mu;                // log mean of negative binomial
  real phi = 1.0 / invphi;         // dispersion parameter
  vector[n_states] state_effect;   // centered state effects
  
  // Center state effects to sum to zero
  state_effect = state_raw - mean(state_raw);
  
  // Calculate log_mu following Schell\'s approach
  log_mu = alpha + 
           delta1 * (l1 - offset1 - P1 * beta) + 
           delta2 * (l2 - offset2 - P2 * beta) + 
           P * beta + 
           X * gamma_X + 
           T * zeta;
           
  // Add state effects and offset
  for (i in 1:N) {
    log_mu[i] = log_mu[i] + state_effect[S[i]] + offset[i];
  }
}

model {
  // Priors - using the same as Schell et al.
  alpha ~ normal(0, sqrt(10));
  delta1 ~ normal(0.5, 1);         // Already constrained 0-1 in parameters block
  delta2 ~ normal(0, 1);           // Already constrained -1 to 1 in parameters block
  
  beta ~ normal(0, prior_sd);
  gamma_X ~ normal(0, 0.1);
  zeta ~ normal(0, 0.2);
  
  // Dispersion and state effects
  invphi ~ normal(0, 0.1) T[0,];
  state_raw ~ normal(0, 0.5);
  
  // Likelihood
  y ~ neg_binomial_2_log(log_mu, phi);
}

generated quantities {
  vector[N] log_lik;
  
  // Calculate IRRs for time-varying effects
  real effect_0 = beta[1];
  real effect_1 = beta[1] + 0.2 * beta[2];
  real effect_2 = beta[1] + 0.4 * beta[2];
  real effect_3 = beta[1] + 0.6 * beta[2];
  real effect_4 = beta[1] + 0.8 * beta[2];
  real effect_5 = beta[1] + beta[2];
  
  real IRR_0 = exp(effect_0);
  real IRR_1 = exp(effect_1);
  real IRR_2 = exp(effect_2);
  real IRR_3 = exp(effect_3);
  real IRR_4 = exp(effect_4);
  real IRR_5 = exp(effect_5);
  
  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_log_lpmf(y[i] | log_mu[i], phi);
  }
}
'

# Write the Stan model to a file
writeLines(placebo_stan_code, "erpo_model.stan")

# 9) Compile the Stan model
placebo_stan_model <- stan_model(file = "erpo_model.stan")

# 10) Find optimization-based initialization values
init_opt <- optimizing(placebo_stan_model, data = placebo_stan_data, 
                       algorithm = "LBFGS", 
                       as_vector = FALSE)

# Use these values for all chains (similar to Schell et al.)
init_values <- list(
  alpha = init_opt$par$alpha,
  delta1 = init_opt$par$delta1,
  delta2 = init_opt$par$delta2,
  beta = init_opt$par$beta,
  gamma_X = init_opt$par$gamma_X,
  zeta = init_opt$par$zeta,
  invphi = init_opt$par$invphi,
  state_raw = init_opt$par$state_raw
)

# Make a list of identical initialization values for each chain
init_list <- list(init_values, init_values, init_values)

# 11) Run Stan sampling
placebo_stan_fit <- sampling(
  placebo_stan_model,
  data = placebo_stan_data,
  init = init_list,
  chains = 3,
  iter = 10000,
  warmup = 2000,
  thin = 5,
  control = list(adapt_delta = 0.9, max_treedepth = 14),
  seed = 123
)

# 12) Extract and summarize results
print(placebo_stan_fit, pars = c("beta", "delta1", "delta2", "invphi", 
                         "IRR_0", "IRR_1", "IRR_2", "IRR_3", "IRR_4", "IRR_5"))

# 13) Check convergence with trace plots
mcmc_trace(placebo_stan_fit, pars = c("beta[1]", "beta[2]", "delta1", "delta2", "invphi"))

# 14) Extract posterior samples and create time-varying IRR data
placebo_samples <- extract(placebo_stan_fit)

# Create a data frame for IRR over time
irr_summary <- data.frame(
  year = 0:5,
  median = c(
    median(placebo_samples$IRR_0),
    median(placebo_samples$IRR_1),
    median(placebo_samples$IRR_2),
    median(placebo_samples$IRR_3),
    median(placebo_samples$IRR_4),
    median(placebo_samples$IRR_5)
  ),
  lower = c(
    quantile(placebo_samples$IRR_0, 0.025),
    quantile(placebo_samples$IRR_1, 0.025),
    quantile(placebo_samples$IRR_2, 0.025),
    quantile(placebo_samples$IRR_3, 0.025),
    quantile(placebo_samples$IRR_4, 0.025),
    quantile(placebo_samples$IRR_5, 0.025)
  ),
  upper = c(
    quantile(placebo_samples$IRR_0, 0.975),
    quantile(placebo_samples$IRR_1, 0.975),
    quantile(placebo_samples$IRR_2, 0.975),
    quantile(placebo_samples$IRR_3, 0.975),
    quantile(placebo_samples$IRR_4, 0.975),
    quantile(placebo_samples$IRR_5, 0.975)
  ),
  prob_less_than_1 = c(
    mean(placebo_samples$IRR_0 < 1),
    mean(placebo_samples$IRR_1 < 1),
    mean(placebo_samples$IRR_2 < 1),
    mean(placebo_samples$IRR_3 < 1),
    mean(placebo_samples$IRR_4 < 1),
    mean(placebo_samples$IRR_5 < 1)
  )
)

# Print IRR summary
print(irr_summary)

# 15) Plot the time-varying IRR
irr_plot <- ggplot(irr_summary, aes(x = year, y = median)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(name = "Incident Rate Ratio (IRR)") +
  scale_x_continuous(name = "Years Since Implementation", breaks = 0:5) +
  labs(title = "Effect of ERPO Laws on Firearm Deaths",
       subtitle = "Immediate and phase-in effects over time") +
  theme_minimal()

# Display the plot
print(irr_plot)

# 16) Create trace plot with MCMC samples (your original visualization)
placebo_n_samples <- length(placebo_samples$IRR_0)

# Create a data frame for all posterior samples
irr_samples <- data.frame(
  sample_id = rep(1:placebo_n_samples, 6),
  year = rep(0:5, each = placebo_n_samples),
  irr = c(
    placebo_samples$IRR_0,
    placebo_samples$IRR_1,
    placebo_samples$IRR_2,
    placebo_samples$IRR_3,
    placebo_samples$IRR_4,
    placebo_samples$IRR_5
  )
)

