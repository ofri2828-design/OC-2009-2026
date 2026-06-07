#!/usr/bin/env Rscript

# Temperature-Salinity diagram for the OC2026 cast database.
# This script uses only base R so it can run without installing packages.

input_file <- if (file.exists("output/cast_data.csv")) {
  "output/cast_data.csv"
} else {
  "cast_data.csv"
}

output_file <- "output/ts_diagram.png"
cruise_output_file <- "output/ts_diagram_by_cruise.png"
station_output_file <- "output/ts_diagram_by_station.png"
combined_station_output_file <- "output/ts_by_station/ts_diagrams_all_stations.png"

if (!file.exists(input_file)) {
  stop("Could not find ", input_file, ". Export the Cast_Data sheet to CSV first.")
}

cast <- read.csv(input_file, check.names = FALSE, fileEncoding = "UTF-8")

temp_col <- "Temp [°C]"
sal_col <- "Salinity [PSU]"
depth_col <- "Depth [m]"
bottom_depth_col <- "Bot. Depth [m]"
station_col <- "Station"
cruise_col <- "Cruise"

required <- c(temp_col, sal_col, depth_col, bottom_depth_col, station_col, cruise_col)
missing_cols <- setdiff(required, names(cast))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

ts <- cast[, required]
names(ts) <- c("temperature", "salinity", "depth", "bottom_depth", "station", "cruise")
ts$temperature <- as.numeric(ts$temperature)
ts$salinity <- as.numeric(ts$salinity)
ts$depth <- as.numeric(ts$depth)
ts$bottom_depth <- as.numeric(ts$bottom_depth)
ts <- ts[complete.cases(ts[, c("temperature", "salinity", "depth")]), ]

if (nrow(ts) == 0) {
  stop("No complete temperature, salinity, and depth records were found.")
}

# UNESCO 1983 equation of state at atmospheric pressure.
# Returns sigma-t, density minus 1000 kg m^-3.
sigma_t <- function(S, T) {
  rho_w <- 999.842594 +
    6.793952e-2 * T -
    9.095290e-3 * T^2 +
    1.001685e-4 * T^3 -
    1.120083e-6 * T^4 +
    6.536332e-9 * T^5

  A <- 0.824493 -
    4.0899e-3 * T +
    7.6438e-5 * T^2 -
    8.2467e-7 * T^3 +
    5.3875e-9 * T^4

  B <- -5.72466e-3 +
    1.0227e-4 * T -
    1.6546e-6 * T^2

  C <- 4.8314e-4

  rho_w + A * S + B * S^(3 / 2) + C * S^2 - 1000
}

pad_range <- function(x, fraction = 0.04) {
  r <- range(x, na.rm = TRUE)
  pad <- diff(r) * fraction
  if (pad == 0) pad <- 0.1
  c(r[1] - pad, r[2] + pad)
}

sal_range <- pad_range(ts$salinity)
temp_range <- pad_range(ts$temperature)
sal_grid <- seq(sal_range[1], sal_range[2], length.out = 160)
temp_grid <- seq(temp_range[1], temp_range[2], length.out = 160)
sigma_grid <- outer(sal_grid, temp_grid, Vectorize(function(S, T) sigma_t(S, T)))
sigma_levels <- pretty(range(sigma_grid, na.rm = TRUE), n = 10)

line_length <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) return(0)
  sum(sqrt(diff(x)^2 + diff(y)^2), na.rm = TRUE)
}

line_point <- function(x, y, fraction = 0.55) {
  if (length(x) < 2 || length(y) < 2) {
    return(c(x = x[1], y = y[1]))
  }
  segment_lengths <- sqrt(diff(x)^2 + diff(y)^2)
  cumulative_lengths <- c(0, cumsum(segment_lengths))
  target <- max(cumulative_lengths, na.rm = TRUE) * fraction
  point_index <- which.min(abs(cumulative_lengths - target))
  c(x = x[point_index], y = y[point_index])
}

draw_isopycnals <- function(x, y, z, levels, line_col = "black", label_col = "black", lwd = 0.8, label_cex = 0.9) {
  contour(
    x,
    y,
    z,
    add = TRUE,
    drawlabels = FALSE,
    col = line_col,
    lwd = lwd,
    levels = levels
  )

  isopycnal_lines <- contourLines(x, y, z, levels = levels)
  for (level in unique(vapply(isopycnal_lines, function(item) item$level, numeric(1)))) {
    level_lines <- isopycnal_lines[vapply(isopycnal_lines, function(item) item$level == level, logical(1))]
    longest <- level_lines[[which.max(vapply(level_lines, function(item) line_length(item$x, item$y), numeric(1)))]]
    label_point <- line_point(longest$x, longest$y)
    label_y_offset <- 0.018 * diff(par("usr")[3:4])
    text(
      label_point[["x"]],
      label_point[["y"]] + label_y_offset,
      labels = format(round(level, 2), trim = TRUE),
      col = label_col,
      cex = label_cex
    )
  }
}

depth_breaks <- pretty(ts$depth, n = 8)
depth_breaks <- depth_breaks[depth_breaks >= min(ts$depth) & depth_breaks <= max(ts$depth)]
if (length(depth_breaks) < 2) {
  depth_breaks <- seq(min(ts$depth), max(ts$depth), length.out = 6)
}

depth_pal <- hcl.colors(100, palette = "Viridis", rev = TRUE)
depth_scaled <- (ts$depth - min(ts$depth)) / diff(range(ts$depth))
depth_scaled[!is.finite(depth_scaled)] <- 0
point_cols <- depth_pal[pmax(1, pmin(100, floor(depth_scaled * 99) + 1))]

png(output_file, width = 2200, height = 1700, res = 220)
op <- par(
  mar = c(5.2, 5.4, 4.5, 6.8),
  xaxs = "i",
  yaxs = "i",
  family = "sans",
  cex.lab = 1.3,
  cex.axis = 1.18,
  cex.main = 1.35
)

plot(
  ts$salinity,
  ts$temperature,
  type = "n",
  xlim = sal_range,
  ylim = temp_range,
  xlab = "Salinity [PSU]",
  ylab = "Temperature [°C]",
  main = "Temperature-Salinity Diagram, OC2026 Cast Data"
)

draw_isopycnals(
  sal_grid,
  temp_grid,
  sigma_grid,
  levels = sigma_levels,
  line_col = "black",
  label_col = "black",
  lwd = 0.8,
  label_cex = 0.92
)

points(
  ts$salinity,
  ts$temperature,
  pch = 16,
  cex = 1.25,
  col = adjustcolor(point_cols, alpha.f = 0.78)
)

box()
mtext(expression(sigma[t] ~ "isopycnals [kg " * m^-3 * "]"), side = 3, line = 0.35, cex = 1.02)

legend_y <- seq(0, 1, length.out = 100)
usr <- par("usr")
x0 <- usr[2] + 0.018 * diff(usr[1:2])
x1 <- usr[2] + 0.048 * diff(usr[1:2])
y0 <- usr[3]
y1 <- usr[4]
for (i in seq_len(99)) {
  rect(
    x0,
    y1 - legend_y[i + 1] * diff(usr[3:4]),
    x1,
    y1 - legend_y[i] * diff(usr[3:4]),
    col = depth_pal[i],
    border = NA,
    xpd = NA
  )
}
axis(
  side = 4,
  at = y1 - (depth_breaks - min(ts$depth)) / diff(range(ts$depth)) * diff(usr[3:4]),
  labels = depth_breaks,
  las = 1,
  pos = x1,
  xpd = NA
)
mtext("Depth [m]", side = 4, line = 5.2)

legend(
  "bottomright",
  legend = sprintf("%s records; %s stations; cruises: %s",
                   format(nrow(ts), big.mark = ","),
                   length(unique(ts$station)),
                   paste(unique(ts$cruise), collapse = ", ")),
  bty = "n",
  cex = 0.78
)

par(op)
dev.off()

message("Wrote ", output_file)

cruises <- sort(unique(ts$cruise))
cruise_cols <- setNames(hcl.colors(length(cruises), palette = "Dark 3"), cruises)
cruise_pch <- setNames(rep(c(16, 17, 15, 18, 8, 3, 4), length.out = length(cruises)), cruises)

png(cruise_output_file, width = 2200, height = 1700, res = 220)
op <- par(
  mar = c(5.2, 5.4, 4.5, 2),
  xaxs = "i",
  yaxs = "i",
  family = "sans",
  cex.lab = 1.3,
  cex.axis = 1.18,
  cex.main = 1.35
)

plot(
  ts$salinity,
  ts$temperature,
  type = "n",
  xlim = sal_range,
  ylim = temp_range,
  xlab = "Salinity [PSU]",
  ylab = "Temperature [°C]",
  main = "Temperature-Salinity Diagram by Cruise"
)

draw_isopycnals(
  sal_grid,
  temp_grid,
  sigma_grid,
  levels = sigma_levels,
  line_col = "black",
  label_col = "black",
  lwd = 0.8,
  label_cex = 0.92
)

for (cruise in cruises) {
  cruise_data <- ts[ts$cruise == cruise, ]
  points(
    cruise_data$salinity,
    cruise_data$temperature,
    pch = cruise_pch[[cruise]],
    cex = 1.25,
    col = adjustcolor(cruise_cols[[cruise]], alpha.f = 0.74)
  )
}

box()
mtext(expression(sigma[t] ~ "isopycnals [kg " * m^-3 * "]"), side = 3, line = 0.35, cex = 1.02)

legend(
  "topright",
  legend = sprintf("%s (n=%s)", cruises, vapply(cruises, function(x) sum(ts$cruise == x), integer(1))),
  col = cruise_cols,
  pch = cruise_pch,
  pt.cex = 1,
  bty = "n",
  title = "Cruise",
  cex = 0.86
)

par(op)
dev.off()

message("Wrote ", cruise_output_file)

preferred_station_order <- c("H15", "H50", "H100", "H400", "H800")
stations <- c(
  preferred_station_order[preferred_station_order %in% unique(ts$station)],
  sort(setdiff(unique(ts$station), preferred_station_order))
)
station_cols <- setNames(hcl.colors(length(stations), palette = "Set 2"), stations)
station_pch <- setNames(rep(c(16, 17, 15, 18, 8, 3, 4), length.out = length(stations)), stations)
station_depth_scaled <- sqrt((ts$depth - min(ts$depth)) / diff(range(ts$depth)))
station_depth_scaled[!is.finite(station_depth_scaled)] <- 0
station_point_cex <- 0.9 + station_depth_scaled * 1.65

png(station_output_file, width = 2200, height = 1700, res = 220)
op <- par(
  mar = c(5.2, 5.4, 4.5, 2),
  xaxs = "i",
  yaxs = "i",
  family = "sans",
  cex.lab = 1.3,
  cex.axis = 1.18,
  cex.main = 1.35
)

plot(
  ts$salinity,
  ts$temperature,
  type = "n",
  xlim = sal_range,
  ylim = temp_range,
  xlab = "Salinity [PSU]",
  ylab = "Temperature [°C]",
  main = "Temperature-Salinity Diagram by Station"
)

draw_isopycnals(
  sal_grid,
  temp_grid,
  sigma_grid,
  levels = sigma_levels,
  line_col = "black",
  label_col = "black",
  lwd = 0.8,
  label_cex = 0.92
)

for (station in stations) {
  station_data <- ts[ts$station == station, ]
  points(
    station_data$salinity,
    station_data$temperature,
    pch = station_pch[[station]],
    cex = station_point_cex[ts$station == station],
    col = adjustcolor(station_cols[[station]], alpha.f = 0.74)
  )
}

box()
mtext(expression(sigma[t] ~ "isopycnals [kg " * m^-3 * "]"), side = 3, line = 0.35, cex = 1.02)

legend(
  "topright",
  legend = sprintf("%s (n=%s)", stations, vapply(stations, function(x) sum(ts$station == x), integer(1))),
  col = station_cols,
  pch = station_pch,
  pt.cex = 1,
  bty = "n",
  title = "Station",
  cex = 0.86
)

depth_legend <- pretty(ts$depth, n = 4)
depth_legend <- depth_legend[depth_legend >= min(ts$depth) & depth_legend <= max(ts$depth)]
depth_legend_cex <- 0.45 + sqrt((depth_legend - min(ts$depth)) / diff(range(ts$depth))) * 1.15
legend(
  "bottomright",
  legend = paste(depth_legend, "m"),
  pt.cex = depth_legend_cex,
  pch = 16,
  col = adjustcolor("grey20", alpha.f = 0.65),
  bty = "n",
  title = "Depth",
  cex = 0.82
)

par(op)
dev.off()

message("Wrote ", station_output_file)

station_dir <- "output/ts_by_station"
if (!dir.exists(station_dir)) {
  dir.create(station_dir, recursive = TRUE)
}

local_depth_scale_stations <- c("H15", "H50", "H100", "H400")

for (station in stations) {
  station_data <- ts[ts$station == station, ]
  station_file <- file.path(station_dir, paste0("ts_diagram_", station, ".png"))
  if (station %in% local_depth_scale_stations) {
    station_bottom <- max(station_data$bottom_depth, station_data$depth, na.rm = TRUE)
    station_depth_limits <- c(0, station_bottom)
    station_depth_breaks <- pretty(station_depth_limits, n = 6)
    station_depth_breaks <- station_depth_breaks[
      station_depth_breaks >= station_depth_limits[1] &
        station_depth_breaks <= station_depth_limits[2]
    ]
  } else {
    station_depth_limits <- range(ts$depth, na.rm = TRUE)
    station_depth_breaks <- depth_breaks
  }
  station_depth_scaled <- (station_data$depth - station_depth_limits[1]) / diff(station_depth_limits)
  station_depth_scaled[!is.finite(station_depth_scaled)] <- 0
  station_depth_cols <- depth_pal[pmax(1, pmin(100, floor(station_depth_scaled * 99) + 1))]

  png(station_file, width = 2200, height = 1700, res = 220)
  op <- par(
    mar = c(5.2, 5.4, 4.5, 6.8),
    xaxs = "i",
    yaxs = "i",
    family = "sans",
    cex.lab = 1.3,
    cex.axis = 1.18,
    cex.main = 1.35
  )

  plot(
    ts$salinity,
    ts$temperature,
    type = "n",
    xlim = sal_range,
    ylim = temp_range,
    xlab = "Salinity [PSU]",
    ylab = "Temperature [°C]",
    main = station
  )

  draw_isopycnals(
    sal_grid,
    temp_grid,
    sigma_grid,
    levels = sigma_levels,
    line_col = "black",
    label_col = "black",
    lwd = 0.8,
    label_cex = 0.92
  )

  points(
    station_data$salinity,
    station_data$temperature,
    pch = 16,
    cex = 1.25,
    col = adjustcolor(station_depth_cols, alpha.f = 0.78)
  )

  box()
  mtext(expression(sigma[t] ~ "isopycnals [kg " * m^-3 * "]"), side = 3, line = 0.35, cex = 1.02)

  legend_y <- seq(0, 1, length.out = 100)
  usr <- par("usr")
  x0 <- usr[2] + 0.018 * diff(usr[1:2])
  x1 <- usr[2] + 0.048 * diff(usr[1:2])
  y0 <- usr[3]
  y1 <- usr[4]
  for (i in seq_len(99)) {
    rect(
      x0,
      y1 - legend_y[i + 1] * diff(usr[3:4]),
      x1,
      y1 - legend_y[i] * diff(usr[3:4]),
      col = depth_pal[i],
      border = NA,
      xpd = NA
    )
  }
  axis(
    side = 4,
    at = y1 - (station_depth_breaks - station_depth_limits[1]) / diff(station_depth_limits) * diff(usr[3:4]),
    labels = station_depth_breaks,
    las = 1,
    pos = x1,
    xpd = NA
  )
  mtext("Depth [m]", side = 4, line = 5.2)

  par(op)
  dev.off()

  message("Wrote ", station_file)
}

png(combined_station_output_file, width = 3300, height = 2200, res = 220)
op <- par(
  mfrow = c(2, 3),
  mar = c(4.8, 5.0, 3.0, 5.0),
  oma = c(0, 0, 0, 0),
  xaxs = "i",
  yaxs = "i",
  family = "sans",
  cex.lab = 1.24,
  cex.axis = 1.12,
  cex.main = 1.27
)

for (i_station in seq_along(stations)) {
  station <- stations[[i_station]]
  station_data <- ts[ts$station == station, ]
  if (station %in% local_depth_scale_stations) {
    station_bottom <- max(station_data$bottom_depth, station_data$depth, na.rm = TRUE)
    station_depth_limits <- c(0, station_bottom)
    station_depth_breaks <- pretty(station_depth_limits, n = 5)
    station_depth_breaks <- station_depth_breaks[
      station_depth_breaks >= station_depth_limits[1] &
        station_depth_breaks <= station_depth_limits[2]
    ]
  } else {
    station_depth_limits <- range(ts$depth, na.rm = TRUE)
    station_depth_breaks <- depth_breaks
  }
  station_depth_scaled <- (station_data$depth - station_depth_limits[1]) / diff(station_depth_limits)
  station_depth_scaled[!is.finite(station_depth_scaled)] <- 0
  station_depth_cols <- depth_pal[pmax(1, pmin(100, floor(station_depth_scaled * 99) + 1))]

  plot(
    ts$salinity,
    ts$temperature,
    type = "n",
    xlim = sal_range,
    ylim = temp_range,
    xlab = "Salinity [PSU]",
    ylab = "Temperature [°C]",
    main = station
  )

  draw_isopycnals(
    sal_grid,
    temp_grid,
    sigma_grid,
    levels = sigma_levels,
    line_col = "black",
    label_col = "black",
    lwd = 0.7,
    label_cex = 0.78
  )

  points(
    station_data$salinity,
    station_data$temperature,
    pch = 16,
    cex = 1.05,
    col = adjustcolor(station_depth_cols, alpha.f = 0.78)
  )

  box()
  mtext(expression(sigma[t] ~ "isopycnals [kg " * m^-3 * "]"), side = 3, line = -0.2, cex = 0.82)
  mtext(LETTERS[[i_station]], side = 3, adj = 0, line = 0.35, font = 2, cex = 1.4)

  legend_y <- seq(0, 1, length.out = 100)
  usr <- par("usr")
  x0 <- usr[2] + 0.022 * diff(usr[1:2])
  x1 <- usr[2] + 0.052 * diff(usr[1:2])
  y0 <- usr[3]
  y1 <- usr[4]
  for (i in seq_len(99)) {
    rect(
      x0,
      y1 - legend_y[i + 1] * diff(usr[3:4]),
      x1,
      y1 - legend_y[i] * diff(usr[3:4]),
      col = depth_pal[i],
      border = NA,
      xpd = NA
    )
  }
  axis(
    side = 4,
    at = y1 - (station_depth_breaks - station_depth_limits[1]) / diff(station_depth_limits) * diff(usr[3:4]),
    labels = station_depth_breaks,
    las = 1,
    pos = x1,
    cex.axis = 0.94,
    xpd = NA
  )
  mtext("Depth [m]", side = 4, line = 4.0, cex = 0.98)
}

if (length(stations) < 6) {
  for (unused_panel in seq_len(6 - length(stations))) {
    plot.new()
  }
}

par(op)
dev.off()

message("Wrote ", combined_station_output_file)
