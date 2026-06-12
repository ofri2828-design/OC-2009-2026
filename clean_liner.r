required_packages <- c("readr", "dplyr", "ggplot2", "oce", "patchwork", "segmented", "zoo")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  install.packages(new_packages)
}

library(readr)
library(dplyr)
library(ggplot2)
library(oce)
library(patchwork)
library(segmented)
library(zoo)

# 1. Define the file path
file_path <- "OC_allyears_combined_new.csv"

# 2. Read the data
df <- read_csv(file_path)

# 3. Initial filtering and median filter despiking
df <- df %>%
  filter(Depth >= 0 & Depth <= 10) %>%
  filter(Salinity < 41.5) %>%
  arrange(Year, `Day of Year`, Depth) %>%
  group_by(Year, Cruise, Station) %>%
  mutate(
    Temperature_raw = Temperature,
    Salinity_raw = Salinity,
    Temperature = zoo::rollmedian(Temperature, k = 5, fill = NA, align = "center"),
    Salinity = zoo::rollmedian(Salinity, k = 5, fill = NA, align = "center")
  ) %>%
  ungroup() %>%
  filter(!is.na(Temperature), !is.na(Salinity))

# 4. Calculate Brunt-Väisälä Frequency
df <- df %>%
  filter(Depth >= 0 & Depth <= 10) %>%
  group_by(Year, Cruise, Station) %>%
  mutate(
    SigmaTheta = swSigmaTheta(
      salinity = Salinity,
      temperature = Temperature,
      pressure = Depth,
      eos = "unesco"
    ),
    BV_Frequency = if(n() > 1) swN2(pressure = Depth, sigmaTheta = SigmaTheta) else NA
  ) %>%
  ungroup() %>%
  filter(is.finite(BV_Frequency))

# 5. Create the Graphs

# Picture 1: Graphs with letters and color by year
p_sal_letters <- ggplot(df, aes(x = `Day of Year`, y = Salinity, color = factor(Year), shape = factor(Year))) +
  geom_point(alpha = 0.8) +
  scale_shape_manual(values = LETTERS) +
  scale_color_viridis_d(option = "turbo") +
  labs(title = "A", x = "Day of Year", y = "Salinity (PSU)", color = "Year", shape = "Year") +
  theme_minimal() +
  coord_cartesian(xlim = c(50, 180))

p_temp_letters <- ggplot(df, aes(x = `Day of Year`, y = Temperature, color = factor(Year), shape = factor(Year))) +
  geom_point(alpha = 0.8) +
  scale_shape_manual(values = LETTERS) +
  scale_color_viridis_d(option = "turbo") +
  labs(title = "B", x = "Day of Year", y = "Temperature (°C)", color = "Year", shape = "Year") +
  theme_minimal() +
  coord_cartesian(xlim = c(50, 180))

p_bv_letters <- ggplot(df, aes(x = `Day of Year`, y = BV_Frequency, color = factor(Year), shape = factor(Year))) +
  geom_point(alpha = 0.8) +
  scale_shape_manual(values = LETTERS) +
  scale_color_viridis_d(option = "turbo") +
  scale_y_log10() +
  labs(title = "C", x = "Day of Year", y = expression("Brunt-Väisälä Freq "*(s^-2)), color = "Year", shape = "Year") +
  theme_minimal() +
  coord_cartesian(xlim = c(50, 180))

# Picture 2: Graphs with linear regression lines
p_sal_lm <- ggplot(df, aes(x = `Day of Year`, y = Salinity)) +
  geom_point(color = "deepskyblue", alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", fill = "lightblue", se = TRUE, level = 0.95) +
  labs(title = "A", x = "Day of Year", y = "Salinity (PSU)") +
  theme_minimal() +
  coord_cartesian(xlim = c(50, 180))

p_temp_lm <- ggplot(df, aes(x = `Day of Year`, y = Temperature)) +
  geom_point(color = "salmon", alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", fill = "mistyrose", se = TRUE, level = 0.95) +
  labs(title = "B", x = "Day of Year", y = "Temperature (°C)") +
  theme_minimal() +
  coord_cartesian(xlim = c(50, 180))

p_bv_lm <- ggplot(df, aes(x = `Day of Year`, y = BV_Frequency)) +
  geom_point(color = "palegreen3", alpha = 0.6) +
  geom_smooth(method = "lm", color = "darkgreen", fill = "palegreen1", se = TRUE, level = 0.95) +
  scale_y_log10() +
  labs(title = "C", x = "Day of Year", y = expression("Brunt-Väisälä Freq "*(s^-2))) +
  theme_minimal() +
  coord_cartesian(xlim = c(50, 180))

# Picture 3: Graphs by Bottom Depth
p_sal_depth <- ggplot(df, aes(x = `Day of Year`, y = Salinity)) +
  geom_point(aes(color = `Bot. Depth`), alpha = 0.8) +
  scale_color_viridis_c(option = "turbo", direction = -1, guide = guide_colorbar(reverse = TRUE)) +
  labs(title = "A", x = "Day of Year", y = "Salinity (PSU)", color = "Bottom Depth") +
  theme_minimal() +
  coord_cartesian(xlim = c(50, 180))

p_temp_depth <- ggplot(df, aes(x = `Day of Year`, y = Temperature)) +
  geom_point(aes(color = `Bot. Depth`), alpha = 0.8) +
  scale_color_viridis_c(option = "turbo", direction = -1, guide = guide_colorbar(reverse = TRUE)) +
  labs(title = "B", x = "Day of Year", y = "Temperature (°C)", color = "Bottom Depth") +
  theme_minimal() +
  coord_cartesian(xlim = c(50, 180))

p_bv_depth <- ggplot(df, aes(x = `Day of Year`, y = BV_Frequency)) +
  geom_point(aes(color = `Bot. Depth`), alpha = 0.8) +
  scale_color_viridis_c(option = "turbo", direction = -1, guide = guide_colorbar(reverse = TRUE)) +
  scale_y_log10() +
  labs(title = "C", x = "Day of Year", y = expression("Brunt-Väisälä Freq "*(s^-2)), color = "Bottom Depth") +
  theme_minimal() +
  coord_cartesian(xlim = c(50, 180))

# 6. Display and save the plots together using the patchwork package
combined_letters <- (p_sal_letters / p_temp_letters / p_bv_letters) + plot_layout(guides = "collect")
combined_lm <- (p_sal_lm / p_temp_lm / p_bv_lm) + plot_layout(guides = "collect")
combined_depth <- (p_sal_depth / p_temp_depth / p_bv_depth) + plot_layout(guides = "collect")

print(combined_letters)
print(combined_lm)
print(combined_depth)

ggsave("combined_letters_plot1.png", plot = combined_letters, width = 8, height = 10, dpi = 300)
ggsave("combined_lm_plot1.png", plot = combined_lm, width = 8, height = 10, dpi = 300)
ggsave("combined_depth_plot1.png", plot = combined_depth, width = 8, height = 10, dpi = 300)

# 7. Predict the transition using Segmented Regression
df$Day_of_Year <- df$`Day of Year`

lm_salinity <- lm(Salinity ~ Day_of_Year, data = df)
lm_temp <- lm(Temperature ~ Day_of_Year, data = df)
lm_bv <- lm(BV_Frequency ~ Day_of_Year, data = df)

set.seed(123)
seg_sal <- segmented(lm_salinity, seg.Z = ~ Day_of_Year)
seg_temp <- segmented(lm_temp, seg.Z = ~ Day_of_Year)
seg_bv <- segmented(lm_bv, seg.Z = ~ Day_of_Year)

# Extract Breakpoints and Confidence Intervals
ci_sal <- confint(seg_sal)
ci_temp <- confint(seg_temp)
ci_bv <- confint(seg_bv)

cat("\n--- Segmented Regression: Winter to Spring Transition ---\n")
cat(sprintf("Salinity Breakpoint: Day %.1f [95%% CI: %.1f - %.1f]\n", ci_sal[1, 1], ci_sal[1, 2], ci_sal[1, 3]))
cat(sprintf("Temperature Breakpoint: Day %.1f [95%% CI: %.1f - %.1f]\n", ci_temp[1, 1], ci_temp[1, 2], ci_temp[1, 3]))
cat(sprintf("BV Frequency Breakpoint: Day %.1f [95%% CI: %.1f - %.1f]\n", ci_bv[1, 1], ci_bv[1, 2], ci_bv[1, 3]))
cat("---------------------------------------------------------\n")

# Add fitted values to dataframe for plotting
df$seg_sal_fit <- fitted(seg_sal)
df$seg_temp_fit <- fitted(seg_temp)
df$seg_bv_fit <- fitted(seg_bv)

# Plot for Salinity (Segmented)
p_seg_sal <- ggplot(df, aes(x = `Day of Year`, y = Salinity)) +
  geom_point(alpha = 0.4, color = "deepskyblue") +
  geom_line(aes(y = seg_sal_fit), color = "blue", linewidth = 1) +
  geom_vline(xintercept = ci_sal[1, 1], linetype = "dashed", color = "black") +
  labs(title = "A", x = "Day of Year", y = "Salinity (PSU)") +
  theme_minimal() +
  coord_cartesian(xlim = c(50, 180))

# Plot for Temperature (Segmented)
p_seg_temp <- ggplot(df, aes(x = `Day of Year`, y = Temperature)) +
  geom_point(alpha = 0.4, color = "salmon") +
  geom_line(aes(y = seg_temp_fit), color = "red", linewidth = 1) +
  geom_vline(xintercept = ci_temp[1, 1], linetype = "dashed", color = "black") +
  labs(title = "B", x = "Day of Year", y = "Temperature (°C)") +
  theme_minimal() +
  coord_cartesian(xlim = c(50, 180))

# Plot for BV Frequency (Segmented)
p_seg_bv <- ggplot(df, aes(x = `Day of Year`, y = BV_Frequency)) +
  geom_point(alpha = 0.4, color = "palegreen3") +
  geom_line(aes(y = seg_bv_fit), color = "darkgreen", linewidth = 1) +
  geom_vline(xintercept = ci_bv[1, 1], linetype = "dashed", color = "black") +
  scale_y_log10() +
  labs(title = "C", x = "Day of Year", y = expression("Brunt-Väisälä Freq "*(s^-2))) +
  theme_minimal() +
  coord_cartesian(xlim = c(50, 180))

combined_seg <- (p_seg_sal / p_seg_temp / p_seg_bv) + plot_layout(guides = "collect")
print(combined_seg)

ggsave("combined_segmented_plot1.png", plot = combined_seg, width = 8, height = 10, dpi = 300)

# ==========================================
# Statistical Tests to Support the Text (No Graphs)
# ==========================================

# 1. Create a classification for stations (Shallow/Coastal vs Deep/Pelagic)
# Note: The 200m depth threshold is an example; adjust it to match your definition of "shallow stations"
df_stat <- df %>%
  mutate(Zone = ifelse(`Bot. Depth` <= 200, "Shallow", "Deep")) %>%
  mutate(Zone = factor(Zone, levels = c("Deep", "Shallow"))) # 'Deep' is set as the reference group

cat("\n======================================================\n")
cat("1. Thermal Dynamics - Is the spring warming rate faster in shallow waters?\n")
cat("======================================================\n")
# Interaction model: Tests if the temperature slope over time (Day of Year) differs between zones
lm_temp_stat <- lm(Temperature ~ `Day of Year` * Zone, data = df_stat)
print(summary(lm_temp_stat))


cat("\n======================================================\n")
cat("2. Haline Dynamics - Salinity\n")
cat("======================================================\n")
# Interaction model for salinity dynamics
lm_sal_stat <- lm(Salinity ~ `Day of Year` * Zone, data = df_stat)
print(summary(lm_sal_stat))


cat("\n======================================================\n")
cat("3. Stratification Timing (Stability Index - BV Frequency)\n")
cat("======================================================\n")
# Interaction model for Brunt-Väisälä Frequency to test stratification timing
lm_bv_stat <- lm(BV_Frequency ~ `Day of Year` * Zone, data = df_stat)
print(summary(lm_bv_stat))