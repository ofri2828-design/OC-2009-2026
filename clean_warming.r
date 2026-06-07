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
  filter(!is.na(Temperature_filt))

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
# 2. Deep Water Analysis (Depth 200m-400m)
# ==========================================
deep_data <- df %>%
  filter(Depth >= 200 & Depth <= 400) %>%
  filter(!is.na(Temperature_filt))

# Multiple Linear Regression (Controlling for Depth, Bottom Depth, and Month)
lm_deep <- lm(Temperature_filt ~ Year + Depth + `Bot. Depth` + as.factor(Month), data = deep_data)

cat("\n=== Deep Water Model Summary (200-400m) ===\n")
print(summary(lm_deep))

# Deep Water Plot
p_deep <- ggplot(deep_data, aes(x = Year, y = Temperature_filt)) +
  geom_point(aes(color = Depth), alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", se = TRUE, fill = "#ffb3b3", level = 0.95) +
  scale_color_viridis_c(
    option = "turbo",
    direction = -1,
    guide = guide_colorbar(reverse = TRUE)
  ) +
  labs(
    title = "B",
    x = "Year",
    y = "Temperature (°C)",
    color = "Depth (m)"
  ) +
  theme_minimal()

# Combine both plots side-by-side using patchwork
combined_plots <- p_surface / p_deep + plot_layout(guides = "collect")
print(combined_plots)

# Save the combined plots to a high-resolution PNG file
ggsave(
  "warming_trends_plot1.png",
  plot = combined_plots,
  width = 10,
  height = 10,
  dpi = 300
)