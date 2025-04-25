library(ggplot2)
library(dplyr)
library(patchwork)
library(lubridate)

# Your original ERPO laws data
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
                          "2020-07-01", "2016-12-08")
)

# Add implementation year and state abbreviations
erpo_laws$implementation_year <- year(as.Date(erpo_laws$implementation_date))
state_abbreviations <- data.frame(
  state = state.name,
  abbr = state.abb
)
# Add DC and correct Connecticut
state_abbreviations <- rbind(
  state_abbreviations,
  data.frame(state = "District of Columbia", abbr = "DC")
)

# Merge to get state abbreviations
erpo_laws <- erpo_laws %>%
  left_join(state_abbreviations, by = "state") %>%
  # Make sure Connecticut is included
  mutate(implementation_year = ifelse(state == "Connecticut", 1999, implementation_year))

# Create a grid layout matching the NYT image
states_grid <- data.frame(
  abbr = c(
    "AK", "", "", "", "", "", "", "", "", "", "ME", "",
    "", "", "", "", "", "", "", "", "VT", "NH", "MA", "",
    "WA", "MT", "ND", "SD", "MN", "WI", "MI", "", "NY", "CT", "RI", "",
    "OR", "ID", "WY", "NE", "IA", "IL", "IN", "OH", "PA", "NJ", "", "",
    "CA", "NV", "UT", "CO", "KS", "MO", "KY", "WV", "DC", "MD", "DE", "",
    "", "AZ", "NM", "OK", "AR", "TN", "VA", "NC", "", "", "", "",
    "", "", "TX", "LA", "MS", "AL", "GA", "SC", "", "", "", "",
    "HI", "", "", "", "", "", "FL", "", "", "", "", ""
  ),
  col = rep(0:11, 8),
  row = rep(0:7, each = 12)
)

# Remove empty cells from the grid
states_grid <- states_grid %>% filter(abbr != "")

# Define the facet years, now including 1999 for Connecticut
facet_years <- c(1999, 2005, 2016, 2018, 2019, 2020, 2024)

# Create a function to check if a state had implemented ERPO laws by a given year
is_implemented_by_year <- function(state_abbr, year) {
  state_data <- erpo_laws %>% filter(abbr == state_abbr)
  if(nrow(state_data) == 0) return(FALSE)
  return(state_data$implementation_year <= year)
}

# Create a function to check if a state was newly implemented in a specific year
is_newly_implemented_in_year <- function(state_abbr, year) {
  state_data <- erpo_laws %>% filter(abbr == state_abbr)
  if(nrow(state_data) == 0) return(FALSE)
  return(state_data$implementation_year == year)
}

# Create plots for each year
plot_list <- list()

# Create a base plot to ensure consistent sizing
base_plot <- ggplot() +
  coord_equal(xlim = c(-0.5, 11.5), ylim = c(-7.5, 0.5)) +
  theme_void()

for(year in facet_years) {
  # Add implementation status to the grid data
  year_data <- states_grid %>%
    mutate(
      implemented = sapply(abbr, function(s) is_implemented_by_year(s, year)),
      newly_implemented = sapply(abbr, function(s) is_newly_implemented_in_year(s, year))
    )
  
  # Count states with ERPO laws
  implemented_count <- sum(year_data$implemented, na.rm = TRUE)
  newly_implemented_count <- sum(year_data$newly_implemented, na.rm = TRUE)
  
  # Create a plot for this year
  p <- base_plot +
    geom_tile(data = year_data, aes(x = col, y = -row, fill = case_when(
      newly_implemented ~ "Newly Implemented",
      implemented ~ "Previously Implemented",
      TRUE ~ "Not Implemented"
    )), color = "white", width = 0.95, height = 0.95) +
    geom_text(data = year_data, aes(x = col, y = -row, label = abbr, 
                                    color = ifelse(implemented | newly_implemented, "black", "black")), 
              size = 3) +
    scale_fill_manual(values = c(
      "Not Implemented" = "#f1f1f1", 
      "Previously Implemented" = "#0072B2",
      "Newly Implemented" = "#f5bb67"
    )) +
    scale_color_manual(values = c("black", "white"), guide = "none") +
    labs(
      title = as.character(year),
      subtitle = paste0(implemented_count, " states with ERPO laws"),
      caption = ifelse(newly_implemented_count > 0, 
                       paste0(newly_implemented_count, " new implementation", 
                              ifelse(newly_implemented_count > 1, "s", "")), "")
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      plot.caption = element_text(hjust = 0.5, size = 8, color = "gray50"),
      plot.margin = margin(10, 10, 10, 10),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = "gray90")
    )
  
  plot_list[[as.character(year)]] <- p
}

# Create an improved legend plot that better fills the space
legend_plot <- base_plot +
  # Add a title at the same position as the year labels in the other plots
  # Create larger color squares in the center of the plot
  annotate("rect", xmin = 0, xmax = 3, ymin = 0, ymax = 1, fill = "#f1f1f1", color = "white") +
  annotate("rect", xmin = 0, xmax = 3, ymin = -2, ymax = -1, fill = "#0072B2", color = "white") +
  annotate("rect", xmin = 0, xmax = 3, ymin = -4, ymax = -3, fill = "#f5bb67", color = "white") +
  # Add centered, larger text to the right of each square
  annotate("text", x = 4, y = 0.5, label = "No ERPO law", hjust = 0, size = 4) +
  annotate("text", x = 4, y = -1.5, label = "Previously implemented", hjust = 0, size = 4) +
  annotate("text", x = 4, y = -3.5, label = "Newly implemented", hjust = 0, size = 4) +
  # Add subtitle explaining laws
  annotate("text", x = 5.5, y = -7.5, 
           label = "Extreme Risk Protection Order (ERPO) laws\nallow temporary removal of firearms from\nindividuals deemed at risk", 
           hjust = 0.5, vjust = 0, size = 3, lineheight = 1.2) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.background = element_rect(fill = "white", color = "gray90"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Arrange plots in a 2×4 grid
rows <- list()
rows[[1]] <- wrap_plots(plot_list[1:4], ncol = 4, widths = rep(1, 4))
rows[[2]] <- wrap_plots(c(plot_list[5:7], list(legend = legend_plot)), ncol = 4, widths = rep(1, 4))

# Combine the rows
combined_plot <- rows[[1]] / rows[[2]] +
  plot_layout(heights = rep(1, 2)) +
  plot_annotation(
    title = "Extreme Risk Protection Order Laws",
    subtitle = "States with enacted \"red flag\" laws, 1999-2024",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "serif"),
      plot.subtitle = element_text(size = 14, hjust = 0.5, family = "serif"),
      plot.caption = element_text(size = 10, color = "grey40", hjust = 0.5, family = "serif"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )


combined_plot

ggsave("erpo_laws_spread.png", width = 13, height = 9, dpi = 300)


