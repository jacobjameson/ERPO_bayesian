
###############################################################################
#  erpo_analysis_stan.R
#  Bayesian NB‐AR(2) model of ERPO effects on firearm suicides, following
#  Schell 2024 – three policy splines (immediate, 0–5 yr phase-in, 5–10 yr).
###############################################################################

## 0. Packages -----------------------------------------------------------------
library(rstan)        ; rstan_options(auto_write = TRUE)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bayesplot)
library(patchwork)
library(loo)

options(mc.cores = parallel::detectCores())
set.seed(123)

## 1.  Load analysis_data ------------------------------------------------------
if (!exists("analysis_data"))
  load("manuscript/src/erpo_analysis_data.RData")   # <- created earlier

## 2.  Factors & time dummies --------------------------------------------------
state_factor <- as.numeric(factor(analysis_data$state))
year_factor  <- as.numeric(factor(analysis_data$year))

year_dummies <- model.matrix(~ 0 + factor(year), analysis_data)
colnames(year_dummies) <- paste0("year_", sort(unique(analysis_data$year)))
year_dummies <- year_dummies[, -1]                       # drop reference

## 3.  Lagged time dummies (shift within state) -------------------------------
n_obs   <- nrow(analysis_data)
n_years <- ncol(year_dummies)
year_dummies_lag1 <- year_dummies_lag2 <- matrix(0, n_obs, n_years)

for (s in unique(state_factor)) {
  idx <- which(state_factor == s)
  if (length(idx) > 1)
    year_dummies_lag1[idx[-1], ]       <- year_dummies[idx[-length(idx)], ]
  if (length(idx) > 2)
    year_dummies_lag2[idx[-(1:2)], ]   <- year_dummies[idx[1:(length(idx)-2)], ]
}

## 4.  Policy matrices  --------------------------------------------------------
policy_vars <- c("erpo", "erpo_phase_in", "erpo_later_phase")
P  <- as.matrix(analysis_data[ , policy_vars])

make_policy_lag <- function(df, k)
  df %>%
  group_by(state) %>% arrange(state, year) %>%
  mutate(across(all_of(policy_vars), ~lag(.x, k, 0))) %>%
  ungroup() %>% select(all_of(policy_vars)) %>% as.matrix()

P1 <- make_policy_lag(analysis_data, 1)
P2 <- make_policy_lag(analysis_data, 2)

## 5.  Covariates (no shrinkage here) -----------------------------------------
xvars <- c("HFR","under_100_percent","x100_199_percent","x200_399_percent",
           "x400_percent","white","black","hispanic","asian",
           "adults_19_25","adults_26_34","adults_35_54","adults_55_64","x65")
X <- as.matrix(analysis_data[ , xvars])

## 6.  Assemble list for Stan --------------------------------------------------
stan_data <- list(
  N        = n_obs,
  n_states = length(unique(state_factor)),
  n_years  = n_years,
  KP       = ncol(P),
  Kx       = ncol(X),
  
  y        = as.integer(analysis_data$deaths),
  l1       = analysis_data$log_rate_lag1,
  l2       = analysis_data$log_rate_lag2,
  offset   = log(analysis_data$population),
  offset1  = log(analysis_data$population_lag1),
  offset2  = log(analysis_data$population_lag2),
  
  P  = P,  P1 = P1,  P2 = P2,
  X  = X,
  T  = year_dummies,  T1 = year_dummies_lag1,  T2 = year_dummies_lag2,
  S  = state_factor,
  prior_sd = 0.064                                     # Schell prior for suicide
)

## 7.  Write Stan code ---------------------------------------------------------
stan_code <- '
data {
  int<lower=1> N; int<lower=1> n_states; int<lower=1> n_years;
  int<lower=1> KP; int<lower=1> Kx;

  int<lower=0> y[N];
  vector[N] l1; vector[N] l2;
  matrix[N,KP]  P;  matrix[N,KP]  P1;  matrix[N,KP]  P2;
  vector[N] offset; vector[N] offset1; vector[N] offset2;
  matrix[N,Kx] X;
  matrix[N,n_years] T; matrix[N,n_years] T1; matrix[N,n_years] T2;
  int<lower=1,upper=n_states> S[N];
  real<lower=0> prior_sd;
}
parameters {
  real alpha;
  real<lower=0,upper=1>  delta1;
  real<lower=-1,upper=1> delta2;
  vector[KP] beta;
  vector[Kx] gamma_X;
  vector[n_years] zeta;
  real<lower=0> invphi;
  vector[n_states] state_raw;
}
transformed parameters {
  vector[n_states] state_eff = state_raw - mean(state_raw);
  real phi = 1.0 / invphi;
  vector[N] log_mu =
        alpha
      + delta1 .* (l1 - offset1 - P1*beta)
      + delta2 .* (l2 - offset2 - P2*beta)
      + P*beta + X*gamma_X + T*zeta + offset;

  for (i in 1:N) log_mu[i] += state_eff[S[i]];
}
model {
  alpha  ~ normal(0, sqrt(10));
  delta1 ~ normal(0.5, 1);
  delta2 ~ normal(0,   1);

  beta     ~ normal(0, prior_sd);
  gamma_X  ~ normal(0, 0.1);
  zeta     ~ normal(0, 0.2);
  invphi   ~ normal(0, 0.1) T[0,];
  state_raw~ normal(0, 0.5);

  y ~ neg_binomial_2_log(log_mu, phi);
}
generated quantities {
  vector[6] IRR;
  IRR[1] = exp(beta[1]);                     //  year 0
  IRR[2] = exp(beta[1] + 0.2*beta[2]);       //  year 1
  IRR[3] = exp(beta[1] + 0.4*beta[2]);       //  year 2
  IRR[4] = exp(beta[1] + 0.6*beta[2]);       //  year 3
  IRR[5] = exp(beta[1] + 0.8*beta[2]);       //  year 4
  IRR[6] = exp(beta[1] +       beta[2]);     //  >=5 yr
}
'

writeLines(stan_code, "erpo_model.stan")
stan_mod <- stan_model("erpo_model.stan")

## 8.  Initialise at posterior mode (LBFGS) -----------------------------------
opt <- optimizing(stan_mod, data = stan_data, algorithm = "LBFGS",
                  as_vector = FALSE)

init_fun <- function() list(
  alpha     = opt$par$alpha,
  delta1    = opt$par$delta1,
  delta2    = opt$par$delta2,
  beta      = opt$par$beta,
  gamma_X   = opt$par$gamma_X,
  zeta      = opt$par$zeta,
  invphi    = opt$par$invphi,
  state_raw = opt$par$state_raw
)

## 9.  Sample ------------------------------------------------------------------
fit <- sampling(stan_mod, data = stan_data,
                init   = init_fun,
                chains = 3, iter = 10000, warmup = 2000, thin = 5,
                control = list(adapt_delta = 0.9, max_treedepth = 14),
                seed = 123)

print(fit, pars = c("beta","delta1","delta2","invphi","IRR"))

###############################################################################
# 10.  Post-processing, plots, counter-factuals
###############################################################################

## A)  Extract IRR matrix ------------------------------------------------------
draws <- rstan::extract(fit)$IRR             #  iter × 6 matrix
irr_summary <- data.frame(
  year = 0:5,
  median = apply(draws, 2, median),
  lower  = apply(draws, 2, quantile, .025),
  upper  = apply(draws, 2, quantile, .975),
  prob_less_than_1 = apply(draws, 2, function(v) mean(v < 1))
)

## B)  Publication-style plots -------------------------------------------------
journal_colors <- colorRampPalette(c("#F5BB67", "#0072B2"))(6)


###  Trace / spaghetti --------------------------------------------------------
ns   <- nrow(draws)          # iterations (rows)
years <- 0:5
irr_df <- as_tibble(draws) |>
  set_names(paste0("y", years)) |>
  mutate(draw = row_number()) |>
  pivot_longer(starts_with("y"),
               names_to  = "year",
               values_to = "irr",
               names_pattern = "y(\\d)") |>
  mutate(year = as.integer(year))


trace_plot <- ggplot() +
  geom_line(data = irr_df %>% filter(draw %% 20 == 0),
            aes(year, irr, group = draw),
            colour = journal_colors[3], alpha = .09, size = .3) +
  geom_line(data = irr_summary, aes(year, median),
            colour = "black", size = 1.2) +
  # Add median line
  geom_line(
    data = irr_summary,
    aes(x = year, y = median),
    color = 'black', size = 1.2
  ) +
  geom_segment(
    aes(x = -0.5, y = 1, xend = 0, yend = 1),
    color = "black", size = 1.2
  ) +
  geom_segment(
    aes(x = 0, y = 1, xend = 0, yend = irr_summary$median[1]),
    color = "black", size = 1.2
  ) +
  # Add reference line at IRR = 1
  geom_hline(
    yintercept = 1, 
    linetype = "longdash", 
    color = "gray30", 
    size = 0.8
  ) +
  annotate(
    "text", x = 2.2, y = 1.01, 
    label = "No effect (IRR = 1)",
    hjust = 0, vjust = 0, 
    size = 3, color = "black", 
    fontface = "italic"
  ) +
  # Labels and scales
  labs(
    subtitle = "Individual MCMC traces with median effect highlighted",
    x = "Years Since Implementation",
    y = "Incidence Rate Ratio (IRR)"
  ) +
  scale_x_continuous(breaks = 0:5, minor_breaks = NULL) +
  scale_y_continuous(
    breaks = seq(0.8, 1.3, 0.05),
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

###  Density facets -----------------------------------------------------------
dens_plot <- ggplot(irr_df,
                    aes(irr, fill = factor(year), colour = factor(year))) +
  geom_density(alpha = 0.65, size = 0.7) +
  geom_rug(alpha = 0.7, show.legend = FALSE, size = 0.4) +
  geom_vline(xintercept = 1, linetype = "longdash", color = "gray30", size = 0.8) +
  annotate(
    "text", x = 1.01, y = 25, 
    label = "No effect (IRR = 1)",
    hjust = 0, vjust = 0, size = 3, 
    color = "black", angle = 90, 
    fontface = "italic"
  ) +
  labs(
    subtitle = "Incidence Rate Ratio distributions for years 0-5 after implementation",
    x = "Incidence Rate Ratio (IRR)",
    y = "Density"
  ) +
  scale_x_continuous(
    breaks = seq(0.8, 1.3, 0.05),
    minor_breaks = seq(0.65, 1.35, 0.05)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = journal_colors[1:6], name = "") +
  scale_color_manual(values = journal_colors[1:6], name = "") +
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


###  Combine & save -----------------------------------------------------------
combined_plot <- trace_plot + dens_plot
combined_plot
ggsave("erpo_effects_combined2.png", combined_plot,
       width = 12, height = 5, dpi = 300)

## C)  Counter-factual prevented deaths ---------------------------------------
cf_df <- analysis_data %>%
  filter(year >= 2001) %>%
  group_by(state) %>%
  mutate(first_year = min(year[erpo == 1], na.rm = TRUE),
         yrs_since  = pmin(pmax(year - first_year, 0), 5),
         yrs_since  = ifelse(is.finite(first_year), yrs_since, NA)) %>%
  ungroup() %>%
  left_join(irr_summary %>% select(yrs_since = year, IRR = median),
            by = "yrs_since") %>%
  mutate(IRR = ifelse(is.na(IRR), 1, IRR),
         deaths_cf   = deaths / IRR,
         deaths_saved = deaths_cf - deaths)

prevent_summary <- summarise(cf_df,
                             observed      = sum(deaths),
                             counterfactual= sum(deaths_cf),
                             prevented     = sum(deaths_saved))

print(prevent_summary, digits = 0)

# 1.  Identify *universal* adoption year for every state
universal_df <- analysis_data %>% 
  mutate(first_erpo_actual = if_else(erpo == 1, year, NA_real_)) %>% 
  group_by(state) %>% 
  summarise(first_erpo_actual = min(first_erpo_actual, na.rm = TRUE)) %>% 
  mutate(first_erpo_univ = pmin(first_erpo_actual, 2018, na.rm = TRUE),
         first_erpo_univ = if_else(is.infinite(first_erpo_univ), 2018, first_erpo_univ))  # never-adopters get 2018

# 2.  Merge back and create “yrs since” under *universal* adoption
univ_cf <- analysis_data %>% 
  left_join(universal_df, by = "state") %>% 
  mutate(yrs_since_univ = pmax(pmin(year - first_erpo_univ, 5), 0),
         # if adoption hasn’t happened yet, no effect
         yrs_since_univ = if_else(year < first_erpo_univ, NA_real_, yrs_since_univ)) %>% 
  left_join(irr_summary %>% select(yrs_since_univ = year, IRR_univ = median),
            by = "yrs_since_univ") %>% 
  mutate(IRR_univ = if_else(is.na(IRR_univ), 1, IRR_univ),
         deaths_cf_univ = deaths / IRR_univ,
         deaths_saved_univ = deaths_cf_univ - deaths)

# 3.  Summarise
totals_univ <- univ_cf %>% 
  summarise(observed      = sum(deaths),
            counterfactual= sum(deaths_cf_univ),
            prevented     = sum(deaths_saved_univ))

# 4.  Compare to the “actual-world” preventions you already computed
additional_saved <- totals_univ$prevented - prevent_summary$prevented

cat(sprintf(
  "Observed deaths:         %6.0f\nCounter-factual (universal 2018 ERPO): %6.0f\n",
  totals_univ$observed, totals_univ$counterfactual))
cat(sprintf(
  "Prevented deaths under universal adoption: %6.0f\n",
  totals_univ$prevented))
cat(sprintf(
  "Additional lives saved vs. actual world:   %6.0f\n",
  additional_saved))

## D)  Save everything ---------------------------------------------------------
save(fit, irr_summary, trace_plot, dens_plot, combined_plot,
     prevent_summary,
     file = "erpo_analysis_stan_results.RData")


library(rstan)
library(dplyr)

## 1.  pull posterior draws of the linear predictor
log_mu_draws <- rstan::extract(fit, pars = "log_mu")$log_mu  # matrix: iterations × N

## 2.  posterior median prediction for each observation
log_mu_med <- apply(log_mu_draws, 2, median)      # length N vector

## 3.  convert to predicted *rates* (deaths per pop)
mu_med   <- exp(log_mu_med)                       # predicted mean counts
rate_pred <- mu_med / analysis_data$population    # same denominator as observed

## 4.  observed rate
rate_obs  <- analysis_data$deaths / analysis_data$population

## 5.  Bayesian R²  (squared correlation)
R2 <- cor(rate_pred, rate_obs, use = "complete.obs")^2
print(R2)


#------------------
###############################################################################
# 11.  Publication lollipop – observed vs counter-factual (treated states)
###############################################################################
library(scales)   # for comma()

## A)  add 2-letter abbreviations ---------------------------------------------
state_abbs <- tibble(state = state.name, abbr = state.abb) |>
  bind_rows(tibble(state = "District of Columbia", abbr = "DC"))

## B)  totals by state  --------------------------------------------------------
state_totals <- cf_df |>
  group_by(state) |>
  summarise(observed       = sum(deaths),
            counterfactual = sum(deaths_cf),
            prevented      = sum(deaths_saved),
            .groups = "drop") |>
  left_join(state_abbs, by = "state") |>
  # keep only states that *ever* had an ERPO (i.e. prevented > 0)
  filter(prevented > 0) |>
  mutate(prevented = round(prevented),
         state_f   = fct_reorder(abbr, prevented))        # order for plotting

## C)  long form for two lollipop end-points ----------------------------------
plot_df <- state_totals |>
  pivot_longer(c(observed, counterfactual),
               names_to  = "type",
               values_to = "deaths")

cols <- c(observed       = "#0072B2",
          counterfactual = "red")

## D)  beautiful lollipop plot -------------------------------------------------
p_lolli <- ggplot(plot_df,
                  aes(x = deaths, y = state_f, colour = type)) +
  ## stick
  geom_line(aes(group = state_f), linewidth = 1.1, colour = "grey70") +
  ## dots
  geom_point(size = 3) +
  ## right-end label = averted
  geom_text(data = filter(plot_df, type == "counterfactual"),
            aes(label = comma(prevented)),
            hjust  = -0.40, size = 3.3, colour = "grey30") +
  ## scales / labels
  scale_x_log10(labels = comma) +
  scale_colour_manual(values = cols,
                      breaks = c("observed", "counterfactual"),
                      labels = c("Observed deaths",
                                 "Counter-factual deaths")) +
  labs(x = NULL, y = NULL, colour = NULL) +
  theme_classic(base_size = 12) +
  # add annotation with an arrow point to a number
  annotate("segment", x = 30000, xend = 5000,
           y = 15, yend = 15,
           arrow = arrow(length = unit(0.2, "cm")),
           colour = "grey30") +
  annotate("text", x = 100000, y = 15,
           label = "Prevented deaths",
           size = 4, colour = "grey30") +
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
  ) 

print(p_lolli)

## E)  save – PDF (vector) & PNG ----------------------------------------------
ggsave("erpo_lollipop_averted.pdf", p_lolli,
       width = 7.5, height = 6, device = cairo_pdf)
ggsave("erpo_lollipop_averted.png", p_lolli,
       width = 7.5, height = 6, dpi = 300)