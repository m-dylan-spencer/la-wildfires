# ============================================================
# SUPPLEMENT (Space Time Cube ANALYSIS-ONLY): Wildfire x Crime (LA County), 2015–2025
# ============================================================

# ----------------------------------------
# Step 1: Load Packages
# ----------------------------------------
library(sf)
library(dplyr)
library(tmap)

# ----------------------------------------
# Step 2: Load and Reproject LA County Boundary
# ----------------------------------------
la_county <- st_read("") %>%
  st_transform(5070)  # Equal-area projection for spatial accuracy

# Break multi-polygon into individual parts
la_parts <- st_cast(la_county, "POLYGON")

# Retain only the largest polygon (assumed mainland)
mainland_la <- la_parts %>%
  mutate(area = st_area(.)) %>%
  slice_max(order_by = area, n = 1) %>%
  select(-area)

# ----------------------------------------
# Step 3. Create 1km x 1km grid
# ----------------------------------------
la_grid <- st_make_grid(mainland_la, cellsize = 1000, square = TRUE) %>%
  st_sf(grid_id = 1:length(.)) %>%
  st_intersection(mainland_la)

tmap_mode("view")

tm_shape(la_grid) +
  tm_borders(col = "gray") +
  tm_shape(mainland_la) +
  tm_borders(col = "black", lwd = 2)

# ----------------------------------------
# Step 4: Load and Combine Crime Data
# ----------------------------------------
crime_city <- st_read("") %>%
  st_transform(5070) %>%
  st_zm(drop = TRUE)

crime_county <- st_read("") %>%
  st_transform(5070) %>%
  st_zm(drop = TRUE)

all_crime <- bind_rows(crime_city, crime_county)

# ----------------------------------------
# Step 5: Aggregate Crimes to Grid Cells by Time
# ----------------------------------------
library(lubridate)
library(tidyr)
library(tidyverse)

all_crime <- all_crime %>%
  mutate(
    date_parsed = parse_date_time(date_cc, orders = c("mdy HM", "mdy HMS", "mdy")),
    year = year(date_parsed),
    month = month(date_parsed),
    year_month = floor_date(date_parsed, unit = "month")  # for monthly aggregation
  )

# Spatial join: assign each crime to a grid cell
crime_joined <- st_join(all_crime, la_grid, left = FALSE)

# Count crimes per grid cell per month
crime_counts <- crime_joined %>%
  st_drop_geometry() %>%
  group_by(grid_id, year_month) %>%
  summarise(n_crimes = n(), .groups = "drop")

write_csv(
  crime_counts,
  ""
)

print(head(crime_counts))

# ----------------------------------------
# Step 6: Create Space-Time Cube Format (Wide Matrix)
# ----------------------------------------
# Pivot wider: one column per time step (e.g., month)
cube_matrix <- crime_counts %>%
  pivot_wider(
    id_cols = grid_id,
    names_from = year_month,
    values_from = n_crimes,
    values_fill = 0  # Fill missing grid-months with zero crimes
  )

print(dim(cube_matrix))      
print(names(cube_matrix)[1:6])  

write_csv(
  cube_matrix,
  ""
)

# ----------------------------------------
# Step 8: Plot Sample Time Series for Grid Cells
# ----------------------------------------
# Reshape matrix back to long format for plotting
cube_long <- cube_matrix %>%
  pivot_longer(
    cols = -grid_id,
    names_to = "date",
    values_to = "n_crimes"
  ) %>%
  mutate(date = as.Date(date))

# Sample a few grid cells (or pick specific IDs)
sample_ids <- sample(unique(cube_long$grid_id), 6)

# ----------------------------------------
# Step 8B: Focus on Grid Cells with More Crime
# ----------------------------------------
# Identify grid cells with at least 20 total crimes across all months
top_cells <- cube_long %>%
  group_by(grid_id) %>%
  summarise(total_crimes = sum(n_crimes, na.rm = TRUE)) %>%
  filter(total_crimes >= 20) %>%
  slice_max(total_crimes, n = 6) %>%
  pull(grid_id)

# ----------------------------------------
# Step 8C: Filter Grid Cells with Activity + Variation
# ----------------------------------------
# Find grid cells with at least 20 total crimes and a standard deviation > 1
top_cells_refined <- cube_long %>%
  group_by(grid_id) %>%
  summarise(
    total_crimes = sum(n_crimes, na.rm = TRUE),
    crime_sd = sd(n_crimes, na.rm = TRUE)
  ) %>%
  filter(total_crimes >= 20, crime_sd > 1) %>%
  slice_max(total_crimes, n = 6) %>%
  pull(grid_id)

# ----------------------------------------
# Step 9: Detect Temporal Trends in Grid Cells (Include Zero-Crime Cells)
# ----------------------------------------
# Reshape cube and fill missing months
crime_long <- cube_matrix %>%
  pivot_longer(-grid_id, names_to = "date", values_to = "crime_count") %>%
  mutate(date = as.Date(date)) %>%
  complete(grid_id, date, fill = list(crime_count = 0))

# Add small value to avoid all-zero time series errors
crime_long <- crime_long %>%
  mutate(crime_count = crime_count + 0.001)

# Run trend regression for each grid cell
trend_summary <- crime_long %>%
  group_by(grid_id) %>%
  summarise(
    slope = coef(lm(crime_count ~ as.numeric(date)))[2],
    p_val = summary(lm(crime_count ~ as.numeric(date)))$coefficients[2, "Pr(>|t|)"],
    .groups = "drop"
  ) %>%
  mutate(
    trend_type = case_when(
      p_val > 0.1 ~ "No Trend",
      slope > 0   ~ "Increasing",
      slope < 0   ~ "Decreasing",
      TRUE ~ "No Trend"
    )
  )

# Ensure every grid cell has a trend label
trend_summary_full <- la_grid %>%
  st_drop_geometry() %>%
  select(grid_id) %>%
  left_join(trend_summary, by = "grid_id") %>%
  mutate(trend_type = replace_na(trend_type, "No Trend"))

# ----------------------------------------
# Step 10: Map Spatiotemporal Trends
# ----------------------------------------
la_grid_trend <- la_grid %>%
  left_join(trend_summary_full, by = "grid_id")

tmap_mode("plot")

tm_shape(la_grid_trend) +
  tm_fill("trend_type", palette = "Dark2", title = "Crime Trend") +
  tm_borders(lwd = 0.2, col = "gray") +
  tm_layout(
    title = "Crime Trends (2015–2025)",
    legend.outside = TRUE
  )

# ----------------------------------------
# Step 11: Recalculate Fire Exposure for All Grid Cells
# ----------------------------------------
fires <- st_read("") %>%
  st_transform(5070)

la_grid_trend <- la_grid_trend %>%
  mutate(
    fire_exposed = ifelse(lengths(st_intersects(., st_union(fires))) > 0, "Yes", "No")
  )

# ----------------------------------------
# Step 12: Prepare Data for Trend × Fire Exposure Analysis
# ----------------------------------------
grid_with_trend <- la_grid_trend

trend_analysis_df <- grid_with_trend %>%
  filter(!is.na(trend_type), !is.na(fire_exposed))

# Confirm structure
str(trend_analysis_df$trend_type)
str(trend_analysis_df$fire_exposed)

# ----------------------------------------
# Step 13: Create Contingency Table and Run Chi-Square Test
# ----------------------------------------
trend_exposure_table <- table(
  Fire_Exposed = trend_analysis_df$fire_exposed,
  Crime_Trend = trend_analysis_df$trend_type
)

print(trend_exposure_table)

chi_result <- chisq.test(trend_exposure_table)

print(chi_result)
print(chi_result$expected)

##################################################################
## Secondary Analysis ####
##################################################################

# ----------------------------------------
# Step 1: Load Packages
# ----------------------------------------
library(sf)
library(dplyr)
library(tmap)
library(lubridate)
library(tidyr)
library(readr)
library(ggplot2)
library(purrr)
library(Kendall)
library(gt)
library(zoo)

# ----------------------------------------
# Step 2: Load and Reproject LA County Boundary
# ----------------------------------------
la_county <- st_read("") %>%
  st_transform(5070)  # Equal-area projection for spatial accuracy

# Break multi-polygon into individual parts
la_parts <- st_cast(la_county, "POLYGON")

# Retain only the largest polygon (assumed mainland)
mainland_la <- la_parts %>%
  mutate(area = st_area(.)) %>%
  slice_max(order_by = area, n = 1) %>%
  select(-area)

# ----------------------------------------
# Step 3: Create 1km x 1km grid
# ----------------------------------------
la_grid <- st_make_grid(mainland_la, cellsize = 1000, square = TRUE) %>%
  st_sf(grid_id = 1:length(.)) %>%
  st_intersection(mainland_la)  # Clip grid to boundary

# ----------------------------------------
# Step 4: Load and Combine Crime Data
# ----------------------------------------
crime_city <- st_read("") %>%
  st_transform(5070) %>%
  st_zm(drop = TRUE)

crime_county <- st_read("") %>%
  st_transform(5070) %>%
  st_zm(drop = TRUE)

all_crime <- bind_rows(crime_city, crime_county)

# ----------------------------------------
# Step 5: Parse dates 
# ----------------------------------------
all_crime <- all_crime %>%
  mutate(
    date_parsed = parse_date_time(
      date_cc,
      orders = c("mdy HMS", "mdy HM", "mdy", "Ymd HMS", "Ymd HM", "Ymd")
    )
  ) %>%
  filter(!is.na(date_parsed))

# ----------------------------------------
# Step 6: MONTHLY cube (pre-fill to avoid values_fill errors)
# ----------------------------------------
crime_joined_m <- st_join(all_crime, la_grid, left = FALSE)

crime_counts_m <- crime_joined_m %>%
  st_drop_geometry() %>%
  mutate(year_month = floor_date(date_parsed, unit = "month")) %>%
  group_by(grid_id, year_month) %>%
  summarise(n_crimes = n(), .groups = "drop")

all_months <- seq(min(crime_counts_m$year_month), max(crime_counts_m$year_month), by = "1 month")

monthly_complete <- tidyr::complete(
  crime_counts_m,
  grid_id,
  year_month = all_months,
  fill = list(n_crimes = 0)
)

write_csv(
  monthly_complete,
  ""
)

cube_matrix_m <- monthly_complete %>%
  mutate(ym_chr = format(year_month, "%Y-%m-01")) %>%
  select(grid_id, ym_chr, n_crimes) %>%
  tidyr::pivot_wider(
    id_cols    = grid_id,
    names_from = ym_chr,
    values_from = n_crimes,
    names_prefix = "m_",
    names_sort   = TRUE
  )

write_csv(
  cube_matrix_m,
  ""
)

# ----------------------------------------
# Step 7: MONTHLY Mann–Kendall trends
# ----------------------------------------
crime_long_m <- cube_matrix_m %>%
  pivot_longer(-grid_id, names_to = "ym_chr", values_to = "crime_count") %>%
  mutate(date = as.Date(sub("^m_", "", ym_chr))) %>%
  select(grid_id, date, crime_count) %>%
  complete(grid_id, date, fill = list(crime_count = 0))

mk_results_m <- crime_long_m %>%
  group_by(grid_id) %>%
  summarise(mk = list(Kendall::MannKendall(crime_count)), .groups = "drop") %>%
  mutate(
    tau   = map_dbl(mk, ~ .x$tau),
    p_val = map_dbl(mk, ~ .x$sl)
  ) %>%
  select(-mk) %>%
  mutate(
    trend_type = case_when(
      p_val <= 0.05 & tau > 0  ~ "Increasing",
      p_val <= 0.05 & tau < 0  ~ "Decreasing",
      TRUE                     ~ "No Trend"
    )
  )

trend_summary_full_m <- la_grid %>%
  st_drop_geometry() %>%
  select(grid_id) %>%
  left_join(mk_results_m, by = "grid_id") %>%
  mutate(trend_type = tidyr::replace_na(trend_type, "No Trend"))

# ----------------------------------------
# Step 8: Fire exposure flag (efficient intersects)
# ----------------------------------------
fires <- st_read("") %>%
  st_transform(5070)

la_grid_trend_m <- la_grid %>%
  left_join(trend_summary_full_m, by = "grid_id") %>%
  mutate(fire_exposed = ifelse(lengths(st_intersects(., fires)) > 0, "Yes", "No"))

# ----------------------------------------
# Step 9: QUARTERLY sensitivity analysis — de-dup + pre-fill (CORRECTED)
# ----------------------------------------
crime_joined_q <- st_join(all_crime, la_grid, left = FALSE)

crime_counts_q <- crime_joined_q %>%
  st_drop_geometry() %>%
  mutate(quarter = floor_date(date_parsed, unit = "quarter")) %>%
  group_by(grid_id, quarter) %>%
  summarise(n_crimes = n(), .groups = "drop")

# Ensure every grid_id has every quarter; fill zeros
all_quarters <- seq(min(crime_counts_q$quarter), max(crime_counts_q$quarter), by = "3 months")

quarterly_complete <- tidyr::complete(
  crime_counts_q,
  grid_id,
  quarter = all_quarters,
  fill = list(n_crimes = 0)
)

# Clean quarter label and collapse any accidental duplicates
quarterly_clean <- quarterly_complete %>%
  mutate(q_chr = format(quarter, "%Y-Q%q")) %>%
  group_by(grid_id, q_chr) %>%
  summarise(n_crimes = sum(n_crimes), .groups = "drop")

# Wide matrix — unique keys; no list-cols
cube_matrix_q <- quarterly_clean %>%
  tidyr::pivot_wider(
    id_cols    = grid_id,
    names_from = q_chr,
    values_from = n_crimes,
    names_prefix = "q_",
    names_sort   = TRUE
  )

# QUARTERLY MK trends
crime_long_q <- cube_matrix_q %>%
  pivot_longer(-grid_id, names_to = "q_chr", values_to = "crime_count") %>%
  mutate(
    date = as.Date(zoo::as.yearqtr(sub("^q_", "", q_chr), format = "%Y-Q%q")),
    crime_count = as.numeric(crime_count)
  ) %>%
  arrange(grid_id, date) %>%
  tidyr::complete(grid_id, date, fill = list(crime_count = 0))

mk_results_q <- crime_long_q %>%
  group_by(grid_id) %>%
  summarise(mk = list(Kendall::MannKendall(crime_count)), .groups = "drop") %>%
  mutate(
    tau   = purrr::map_dbl(mk, ~ .x$tau),
    p_val = purrr::map_dbl(mk, ~ .x$sl),
    trend_type_q = dplyr::case_when(
      p_val <= 0.05 & tau > 0  ~ "Increasing",
      p_val <= 0.05 & tau < 0  ~ "Decreasing",
      TRUE                     ~ "No Trend"
    )
  ) %>%
  select(-mk)

# Join exposure; KEEP all grid cells and fill any NAs
la_grid_trend_q <- la_grid %>%
  st_drop_geometry() %>%
  select(grid_id) %>%
  left_join(mk_results_q, by = "grid_id") %>%
  left_join(st_drop_geometry(la_grid_trend_m)[, c("grid_id", "fire_exposed")], by = "grid_id") %>%
  mutate(
    trend_type_q = tidyr::replace_na(trend_type_q, "No Trend"),
    fire_exposed = tidyr::replace_na(fire_exposed, "No")
  )

# Contingency table + tests
trend_exposure_table_q <- table(
  Fire_Exposed = la_grid_trend_q$fire_exposed,
  Crime_Trend  = la_grid_trend_q$trend_type_q
)
chi_result_q <- suppressWarnings(chisq.test(trend_exposure_table_q))

# If chi-square has small expected counts, also compute Fisher's exact (optional)
fisher_result_q <- tryCatch(fisher.test(trend_exposure_table_q), error = function(e) NULL)

print(trend_exposure_table_q)
print(chi_result_q)
if (!is.null(fisher_result_q)) print(fisher_result_q)

# ----------------------------------------
# Step 10: Cramér’s V for association strength (2x3 + 2x2)
# ----------------------------------------
cramers_v <- function(tab) {
  k <- min(nrow(tab), ncol(tab))
  chisq <- suppressWarnings(chisq.test(tab))$statistic
  as.numeric(sqrt(chisq / (sum(tab) * (k - 1))))
}

cv_2x3_m <- cramers_v(trend_exposure_table_m)
cv_2x3_q <- cramers_v(trend_exposure_table_q)
cv_2x2_m <- cramers_v(tab2_m)
cv_2x2_q <- cramers_v(tab2_q)

cat("\n--- Cramér's V ---\n")
cat(sprintf("Monthly 2x3: %.4f | Quarterly 2x3: %.4f\n", cv_2x3_m, cv_2x3_q))
cat(sprintf("Monthly 2x2: %.4f | Quarterly 2x2: %.4f\n", cv_2x2_m, cv_2x2_q))

# ----------------------------------------
# Step 11: Proportions with 95% CIs (Trend vs No Trend)
# ----------------------------------------
# Monthly exposed trend proportion
n_exp_m        <- sum(tab2_m["Yes", ])
prop_exp_m     <- tab2_m["Yes", "Trend"] / n_exp_m
prop_exp_m_ci  <- prop.test(tab2_m["Yes", "Trend"], n_exp_m)

# Quarterly exposed trend proportion
n_exp_q        <- sum(tab2_q["Yes", ])
prop_exp_q     <- tab2_q["Yes", "Trend"] / n_exp_q
prop_exp_q_ci  <- prop.test(tab2_q["Yes", "Trend"], n_exp_q)

cat("\n--- Proportion of Trend in Fire-Exposed (w/ 95% CI) ---\n")
cat(sprintf("Monthly:  %.4f  (n=%d)\n", prop_exp_m, n_exp_m))
print(prop_exp_m_ci)
cat(sprintf("Quarterly: %.4f  (n=%d)\n", prop_exp_q, n_exp_q))
print(prop_exp_q_ci)

# Optional: difference in proportions (exposed vs unexposed) — Monthly & Quarterly
prop_diff_test <- function(tab2) {
  x <- c(tab2["Yes", "Trend"], tab2["No", "Trend"])
  n <- c(sum(tab2["Yes", ]),    sum(tab2["No", ]))
  prop.test(x, n, correct = FALSE)
}
cat("\n--- Difference in Trend Proportions (Exposed vs Unexposed) ---\n")
print(prop_diff_test(tab2_m))
print(prop_diff_test(tab2_q))

