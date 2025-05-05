# Bayesian DiD & Event‐Study with brms
# ------------------------------------

# 1) Load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(ggplot2)

# 2) Prepare data
did_data <- analysis_data %>%
  group_by(state) %>%
  mutate(
    treated = as.numeric(any(erpo > 0)),
    impl_year = if_else(treated == 1, min(year[erpo > 0]), NA_real_),
    post = if_else(treated == 1 & year >= impl_year, 1, 0),
    rel_year = year - impl_year
  ) %>%
  ungroup() %>%
  mutate(
    rel_year_cat = factor(pmin(pmax(rel_year, -5), 5))  # clamp to [-5,5]
  )

# 3) Bayesian Difference‐in‐Differences
did_fit <- brm(
  bf(deaths ~ offset(log(population)) + treated:post +
       (1 | state) + (1 | year)),
  data    = did_data,
  family  = negbinomial(),
  prior   = c(
    prior(normal(0, 1), class = "b"),      # for treatment effect
    prior(gamma(0.1, 0.1), class = "shape") # NB dispersion
  ),
  chains  = 4, iter = 4000, warmup = 1000,
  seed    = 123
)

# 4) Extract DiD IRR posterior
did_effect <- did_fit %>%
  spread_draws(b_treatedpost) %>%
  median_hdci(exp(b_treatedpost), .width = c(.95)) %>%
  rename(IRR = `exp(b_treatedpost)`)

print(did_effect)

# 5) Plot DiD posterior
did_fit %>%
  conditional_effects("treated:post") %>%
  plot(points = FALSE)[[1]] +
  scale_y_log10() +
  labs(
    title = "Bayesian DiD: IRR for ERPO (treated × post)",
    y = "Incidence Rate Ratio (log scale)"
  )

# 6) Bayesian Event‐Study
es_fit <- brm(
  bf(deaths ~ offset(log(population)) +
       treated * rel_year_cat +
       (1 | state) + (1 | year)),
  data    = did_data,
  family  = negbinomial(),
  prior   = c(
    prior(normal(0, 1), class = "b"),
    prior(gamma(0.1, 0.1), class = "shape")
  ),
  chains  = 4, iter = 4000, warmup = 1000,
  seed    = 123
)

# 7) Extract event‐study IRRs correctly
library(tidybayes)
library(stringr)

es_coefs <- es_fit %>%
  # pull out each treated:rel_year_cat interaction coefficient
  spread_draws(
    `b_treated:rel_year_catM4`, `b_treated:rel_year_catM3`,
    `b_treated:rel_year_catM2`, `b_treated:rel_year_catM1`,
    `b_treated:rel_year_cat0`,  `b_treated:rel_year_cat1`,
    `b_treated:rel_year_cat2`,  `b_treated:rel_year_cat3`,
    `b_treated:rel_year_cat4`,  `b_treated:rel_year_cat5`
  ) %>%
  # pivot into long form
  pivot_longer(
    cols      = starts_with("b_treated:rel_year_cat"),
    names_to  = "param",
    values_to = "beta"
  ) %>%
  # extract the relative‐year from the name and compute IRR
  mutate(
    rel = str_remove(param, "b_treated:rel_year_cat"),
    IRR = exp(beta)
  ) %>%
  # summarize with 80% and 95% credible intervals
  group_by(rel) %>%
  median_hdci(IRR, .width = c(.80, .95)) %>%
  ungroup()

es_coefs


# 8) Plot event‐study IRRs
ggplot(es_coefs, aes(x = rel, y = IRR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper, fill = .width),
                  alpha = 0.4) +
  geom_point() +
  labs(
    title = "Bayesian Event‐Study: ERPO Effect Over Time",
    x = "Years Since Implementation",
    y = "Incidence Rate Ratio",
    fill = "Credible\nInterval"
  ) +
  theme_minimal()







# -------------------------------------------------------------------
# Event‐Study with Callaway & Sant’Anna’s did package
# for ERPO implementation and firearm suicides
# -------------------------------------------------------------------

# 1) Load libraries
library(tidyverse)
library(did)


data <- analysis_data %>%
  select(state, year, erpo, crude_rate, population, population_lag1, population_lag2, log_rate_lag1,log_rate_lag2 ) 

# if a state ever has an ERPO, make a var called treated = 1
data <- data %>%
  group_by(state) %>%
  mutate(treated = as.numeric(any(erpo > 0))) %>%
  mutate(post = if_else(treated == 1 & year >= min(year[erpo > 0]), 1, 0)) %>%
  ungroup()

# create a variable called gname that is the year for a state where treatment erpo starts
data <- data %>%
  group_by(state) %>%
  mutate(impl_year = min(year[erpo > 0])) %>%
  ungroup()

# make impl_year 0 if inf
data <- data %>%
  mutate(impl_year = if_else(impl_year == Inf, 0, impl_year))


data$state <- as.numeric(factor(data$state))


atts_1 <- att_gt(yname = "crude_rate",
                   tname = "year",
                   idname = "state",
                   gname = "impl_year",
                   control_group = "nevertreated",
                   clustervars   = "state",
                   panel = F,
                   xformla = ~log_rate_lag1 + log_rate_lag2 + population,
                   data = data,
                   allow_unbalanced_panel = F,
                   base_period = "universal",
                   cband = FALSE,
                   pl = TRUE,
                   cores = 4)
  
agg_event_1 <- aggte(atts_1, type = "dynamic", min_e = -5, max_e = 4, na.rm = TRUE)

# 5) Summary
summary(agg_event_1)

# 6) Plot
ggdid(agg_event_1) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "black") +
  geom_hline(yintercept =  0, color = "black") +
  geom_point() +
  geom_line() +
  geom_errorbar(width = 0.2) +
  theme_test() +
  theme(
    legend.position    = "none",
    text               = element_text(size = 10, face = "bold"),
    axis.ticks         = element_blank(),
    panel.grid.major.y = element_line(color = "gray", size = 0.25)
  ) +
  ggtitle("Never‐treated states as controls") +
  xlab("Years Since ERPO Implementation") +
  ylab("Firearm Suicides (per 100k)")

