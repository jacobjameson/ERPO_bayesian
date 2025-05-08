library(tidyverse)

# Read in raw files
KFF_econ <- read_csv("raw/state year/KFF econ.csv")
KFF_race <- read_csv("raw/state year/KFF race.csv")
KFF_pop <- read_csv("raw/state year/KFF pop.csv")


# -------- ECONOMIC DATA --------
econ_data <- KFF_econ %>%
  pivot_longer(
    cols = -Location,
    names_to = "year_group",
    values_to = "value"
  ) %>%
  separate(year_group, into = c("year", "income_group"), sep = "__") %>%
  mutate(
    year = as.integer(year),
    state = Location
  ) %>%
  select(state, year, income_group, value) %>%
  pivot_wider(names_from = income_group, values_from = value) %>%
  select(-`NA`, -Total) %>%
  clean_names()

# -------- RACE DATA --------
race_data <- KFF_race %>%
  mutate(across(-Location, as.character)) %>%
  pivot_longer(
    cols = -Location,
    names_to = "year_group",
    values_to = "value_raw"
  ) %>%
  separate(year_group, into = c("year", "race_group"), sep = "__") %>%
  mutate(
    year = as.integer(year),
    state = Location,
    value = as.numeric(value_raw),
    value = if_else(is.na(value), 0.01, value)
  ) %>%
  select(state, year, race_group, value) %>%
  pivot_wider(names_from = race_group, values_from = value) %>%
  select(-`NA`) %>%
  clean_names()

# -------- POP DATA --------

pop_data <- KFF_pop %>%
  mutate(across(-Location, as.character)) %>%
  pivot_longer(
    cols = -Location,
    names_to = "year_group",
    values_to = "value_raw"
  ) %>%
  separate(year_group, into = c("year", "age_group"), sep = "__") %>%
  mutate(
    year = as.integer(year),
    state = Location,
    value = as.numeric(value_raw),
    value = if_else(is.na(value), 0.01, value)
  ) %>%
  select(state, year, age_group, value) %>%
  pivot_wider(names_from = age_group, values_from = value) %>%
  select(-`NA`) %>%
  clean_names()


# -------- MERGE AND FINAL DATA --------
kff_combined <- econ_data %>%
  full_join(race_data, by = c("state", "year")) %>%
  full_join(pop_data, by = c("state", "year")) %>%
  arrange(state, year)

# View the result
print(kff_combined, n = 10)


library(tidyverse)
library(zoo)

# Step 1: Get full panel of all state-year combos
all_years <- tibble(year = 2000:2023)
all_states <- distinct(kff_combined, state)
full_panel <- crossing(all_states, all_years)

# Step 2: Merge with your existing data
kff_full <- full_panel %>%
  left_join(kff_combined, by = c("state", "year")) %>%
  arrange(state, year)

# Step 3: Impute using linear interpolation + extrapolation (rule = 2)
# Select variables to interpolate
vars_to_impute <- names(kff_combined)[!(names(kff_combined) %in% c("state", "year"))]

# Apply na.approx to each state group
kff_imputed <- kff_full %>%
  group_by(state) %>%
  mutate(across(all_of(vars_to_impute), ~ na.approx(., x = year, na.rm = FALSE, rule = 2))) %>%
  ungroup()

# View imputed rows (optional)
kff_imputed %>%  print(n = 20)

#--------------------

RAND <- read_excel("raw/state year/RAND.xlsx", 
                   sheet = "State-Level Data & Factor Score") %>%
  select(year=Year, state=STATE, HFR)


RAND_lagged <- RAND %>%
  mutate(year = year + 5) %>%
  filter(year >= 2000) 

all_years <- tibble(year = 2000:2023)
all_states <- distinct(RAND_lagged, state)
full_panel <- crossing(all_states, all_years)

# Step 2: Merge with RAND_lagged
RAND_full <- full_panel %>%
  left_join(RAND_lagged, by = c("state", "year")) %>%
  arrange(state, year)

# Step 3: Extrapolate HFR values to 2022 and 2023
RAND_imputed <- RAND_full %>%
  group_by(state) %>%
  mutate(HFR = na.approx(HFR, x = year, na.rm = FALSE, rule = 2)) %>%
  ungroup()



# Step 4: Merge RAND_imputed with kff_imputed

data <- kff_imputed %>%
  left_join(RAND_imputed, by = c("state", "year")) %>%
  arrange(state, year) 


data
