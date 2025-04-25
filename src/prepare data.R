library(readxl)
library(tidyverse)
library(dplyr)
library(tidyr)
library(lubridate)


FA_1999_2020 <- read_excel("raw/FA_1999_2020.xlsx") %>%
  select(-Notes) %>%
  filter(!is.na(State)) %>%
  filter(Year < 2018)

FA_1999_2020$Deaths <- as.numeric(FA_1999_2020$Deaths)
FA_1999_2020$Deaths[is.na(FA_1999_2020$Deaths)] <- 5

FA_1999_2020$Population <- as.numeric(FA_1999_2020$Population)
FA_1999_2020$`Crude Rate` <- FA_1999_2020$Deaths/FA_1999_2020$Population * 100000

FA_2018_2023 <- read_excel("raw/FA_2018_2023.xlsx") %>%
  select(-Notes)  %>%
  filter(!is.na(State))

FA_2018_2023$Deaths <- as.numeric(FA_2018_2023$Deaths)
FA_2018_2023$Deaths[is.na(FA_2018_2023$Deaths)] <- 5

FA_2018_2023$Population <- as.numeric(FA_2018_2023$Population)
FA_2018_2023$`Crude Rate` <- FA_2018_2023$Deaths/FA_2018_2023$Population * 100000

FA_1999_2023 <- rbind(FA_1999_2020, FA_2018_2023)


FA_1999_2023 <- FA_1999_2023 %>%
  rename_with(~ gsub(" ", "_", .x)) %>%
  rename_with(~ tolower(.x)) %>%
  select(-state_code, -sex_code, -year_code) 

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
    all_states = unique(FA_1999_2023$state)
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
FA_total <- FA_1999_2023 %>%
  group_by(state, year) %>%
  summarize(
    deaths = sum(deaths),
    population = sum(population)/2,  
    crude_rate = sum(deaths) / (sum(population)/1e5)
  )

# Merge the data
merged_data <- left_join(erpo_panel, FA_total, by = c("state", "year"))

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
         state != 'Connecticut')

# Save the prepared data
save(analysis_data, file = "outputs/data/erpo_analysis_data.RData")
