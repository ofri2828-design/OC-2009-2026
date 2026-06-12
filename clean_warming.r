library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(zoo)

# 1. Define the file path and read the data
file_path <- "OC_allyears_combined_new.csv"
df <- read_csv(file_path, show_col_types = FALSE)

# 2. Median filter within each profile
df <- df %>%
  arrange(Year, Cruise, Station, Depth) %>%
  group_by(Year, Cruise, Station) %>%
  mutate(
    Temperature_raw = Temperature,
    Salinity_raw = Salinity,
    Temperature_filt = zoo::rollmedian(Temperature, k = 5, fill = NA, align = "center"),
    Salinity_filt = zoo::rollmedian(Salinity, k = 5, fill = NA, align = "center")
  ) %>%
  ungroup()

# ==========================================
# 1. Surface Water Analysis (Depth 0 - 10m)
# ==========================================
surface_data <- df %>%
  filter(Depth >= 0 & Depth <= 10) %>%
  filter(!is.na(Temperature_filt) & !is.na(Month) & !is.na(`Bot. Depth`))

# Multiple Linear Regression (Controlling for Month and Bottom Depth)
lm_surface <- lm(Temperature_filt ~ Year + as.factor(Month) + `Bot. Depth`, data = surface_data)

cat("\n=== Surface Water Model Summary (0-10m) ===\n")
print(summary(lm_surface))

# Surface Water Plot
p_surface <- ggplot(surface_data, aes(x = Year, y = Temperature_filt)) +
  geom_point(aes(color = as.factor(Month)), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = TRUE, level = 0.95) +
  scale_color_viridis_d(option = "turbo") +
  labs(
    title = "A",
    x = "Year",
    y = "Temperature (°C)",
    color = "Month"
  ) +
  theme_minimal()

# ==========================================
# 2. LIW Analysis (Max Salinity per Profile)
# ==========================================
liw_data <- df %>%
  filter(!is.na(Temperature_filt) & !is.na(Salinity_filt)) %>%
  filter(`Bot. Depth` > 150) %>%
  filter(Depth > 100) %>%
  group_by(Year, Cruise, Station) %>%
  slice_max(order_by = Salinity_filt, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(!is.na(Month) & !is.na(`Bot. Depth`) & !is.na(Depth))

# Multiple Linear Regression (Controlling for Month and Bottom Depth)
lm_liw <- lm(Temperature_filt ~ Year + as.factor(Month) + `Bot. Depth` , data = liw_data)

cat("\n=== LIW Model Summary (Max Salinity) ===\n")
print(summary(lm_liw))

# LIW Plot
p_liw <- ggplot(liw_data, aes(x = Year, y = Temperature_filt)) +
  geom_point(aes(color = as.factor(Month)), alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = TRUE, fill = "#ffb3b3", level = 0.95) +
  scale_color_viridis_d(option = "turbo") +
  labs(
    title = "B",
    x = "Year",
    y = "Temperature (°C)",
    color = "Month"
  ) +
  theme_minimal()

# Combine both plots side-by-side using patchwork
combined_plots <- p_surface / p_liw + plot_layout(guides = "collect")
print(combined_plots)

# Save the combined plots to a high-resolution PNG file
ggsave(
  "warming_trends_plot2.png",
  plot = combined_plots,
  width = 10,
  height = 10,
  dpi = 300
)
