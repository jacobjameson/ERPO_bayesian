# sketch of an IPUMS-CPS workflow in R:
library(ipumsr)
ddi   <- read_ipums_ddi("cps_00001.xml")
data  <- read_ipums_micro(ddi) %>% 
  filter(YEAR >= 2000 & YEAR <= 2023) %>% 
  group_by(STATEFIP, YEAR) %>% 
  summarize(
    pov_rate    = mean(INCTOT < poverty_threshold),  # your poverty flag
    med_hh_inc  = median(HHINCOME),
    pct_unemp   = mean(EMPSTAT == 3),
    pct_bach    = mean(EDUC >= 16),
    pct_children_single = mean(KIDS > 0 & RELATE == 6), # RELATE==6: own child, HH head not spouse
    # …etc…
  )
