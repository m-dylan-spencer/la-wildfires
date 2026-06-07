# ============================================================
# SUPPLEMENT (GWR ANALYSIS-ONLY): Wildfire x Crime (LA County), 2015–2025
# ============================================================

# ----------------------------------------
# Step 1: Load Packages
# ----------------------------------------
library(sf)
library(dplyr)
library(GWmodel)
library(ggplot2)
library(viridis)

# ----------------------------------------
# Step 2: Load Fire Perimeters
# ----------------------------------------
fires <- st_read("") %>%
  st_transform(5070)  # Use equal-area projection for spatial analysis

# ----------------------------------------
# Step 3: Create 2km Buffer Around Fires (study area extent)
# ----------------------------------------
fires_buffer <- st_buffer(fires, dist = 2000)
study_area <- st_union(fires_buffer)

# ----------------------------------------
# Step 4: Generate Grid Over Study Area
# ----------------------------------------
grid <- st_make_grid(study_area, cellsize = 1000, square = TRUE) %>%
  st_sf(grid_id = 1:length(.))

# ----------------------------------------
# Step 5: Load and Filter Crime Data
# ----------------------------------------

crime_city  <- st_read("") %>% st_transform(5070)
crime_county <- st_read("") %>% st_transform(5070)
all_crime <- bind_rows(crime_city, crime_county)

# ----------------------------------------
# Step 6: Count Crimes per Grid Cell
# ----------------------------------------
grid_with_crime <- st_join(grid, all_crime, left = FALSE) %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(crime_count = n()) %>%
  ungroup()

grid <- left_join(grid, grid_with_crime, by = "grid_id") %>%
  mutate(crime_count = tidyr::replace_na(crime_count, 0))

# ----------------------------------------
# Step 7: Tag Fire Exposure (Grid intersects fire perimeter)
# ----------------------------------------
grid <- grid %>%
  mutate(fire_exposed = lengths(st_intersects(., st_union(fires))) > 0)

# ----------------------------------------
# Step 8: Add Population Density (via tidycensus)
# ----------------------------------------
library(tidycensus)
library(tidyr)

# Download 2020 ACS 5-year data for total population by block group
la_pop <- get_acs(
  geography = "block group",
  variables = "B01003_001",  # Total population
  year = 2020,
  state = "CA",
  county = "Los Angeles",
  geometry = TRUE
)

la_pop <- st_transform(la_pop, 5070)

grid <- st_join(grid, la_pop %>% select(GEOID, estimate), left = TRUE)

# Compute population density (people per square meter)
grid <- grid %>%
  rename(population = estimate) %>%
  mutate(
    pop_density = population / as.numeric(sf::st_area(grid)),
    pop_density = tidyr::replace_na(pop_density, 0)
  )
grid <- grid %>%
  mutate(pop_density_km2 = pop_density * 1e6)

# ----------------------------------------
# Step 9: GWR on crime RATE with areal-weighted population
# ----------------------------------------
library(units)

sf::sf_use_s2(FALSE)

grid   <- sf::st_make_valid(grid)
la_pop <- sf::st_make_valid(la_pop)

la_pop_aw <- la_pop %>%
  dplyr::mutate(bg_area = sf::st_area(sf::st_geometry(.))) %>%
  dplyr::select(GEOID, estimate, bg_area)   # 'estimate' = total pop

grid_for_pop <- grid %>% dplyr::select(grid_id)

grid_pop_parts <- sf::st_intersection(grid_for_pop, la_pop_aw) %>%
  dplyr::mutate(
    int_area   = sf::st_area(sf::st_geometry(.)),
    pop_share  = as.numeric(int_area / bg_area),
    pop_alloc  = as.numeric(estimate) * pmin(pmax(pop_share, 0), 1)
  ) %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(grid_id) %>%
  dplyr::summarise(pop_aw = sum(pop_alloc, na.rm = TRUE), .groups = "drop")  # <-- use 'pop_aw'

grid <- grid %>%
  dplyr::select(-dplyr::any_of(c("population", "estimate"))) %>%  # <-- drop old cols if present
  dplyr::left_join(grid_pop_parts, by = "grid_id") %>%
  dplyr::mutate(population = tidyr::replace_na(as.numeric(pop_aw), 0)) %>%  # <-- create clean numeric
  dplyr::select(-pop_aw)

sf::sf_use_s2(TRUE)

grid <- grid %>%
  mutate(
    fire_exposed     = as.numeric(fire_exposed),                 # 0/1 numeric
    crime_rate_1000  = ifelse(population > 0, (crime_count / population) * 1000, 0),
    pop_density_km2  = as.numeric(pop_density_km2)
  )

is_bad <- !is.finite(grid$crime_rate_1000) |
  !is.finite(grid$fire_exposed)    |
  !is.finite(grid$pop_density_km2)
grid_clean <- grid[!is_bad, ]
stopifnot(nrow(grid_clean) > 0)

grid_sp <- as(grid_clean, "Spatial")
dmat    <- GWmodel::gw.dist(dp.locat = coordinates(grid_sp), longlat = FALSE)

set.seed(42)
bw <- GWmodel::bw.gwr(
  formula  = crime_rate_1000 ~ fire_exposed + pop_density_km2,
  data     = grid_sp,
  approach = "AICc",
  kernel   = "bisquare",
  adaptive = TRUE,
  dMat     = dmat
)

gwr_model <- GWmodel::gwr.basic(
  formula  = crime_rate_1000 ~ fire_exposed + pop_density_km2,
  data     = grid_sp,
  bw       = bw,
  kernel   = "bisquare",
  adaptive = TRUE,
  dMat     = dmat
)

gwr_output <- sf::st_as_sf(gwr_model$SDF)

coef_col <- grep("^fire_exposed", names(gwr_output), value = TRUE)[1]
r2_col   <- if ("Local_R2" %in% names(gwr_output)) "Local_R2" else "localR2"

gwr_output$coef_fire <- as.numeric(gwr_output[[coef_col]])                
gwr_output$local_r2  <- pmax(0, pmin(1, as.numeric(gwr_output[[r2_col]])))

####################################################################
#### GWR and Damage Analysis ###
####################################################################

# ----------------------------------------
# Step 0: Prereqs (assumes you already ran the cube + damage code)
# ----------------------------------------
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(GWmodel)
library(tigris)
library(tidycensus)
library(units)

stopifnot(st_crs(la_grid)$epsg == 5070)

# ----------------------------------------
# Step 1: Define study area from fires & subset grid
# ----------------------------------------
fires <- st_transform(large_fires, 5070)
study_area <- st_union(st_buffer(fires, 2000))  # 2km buffer

grid_gwr <- la_grid %>%
  st_intersection(study_area)  # keep IDs consistent

# ----------------------------------------
# Step 2: Build response (crime counts) on the SAME grid
# ----------------------------------------
# Safer point->grid count: join crimes to grid, then tally, then left_join back
crime_city   <- st_read("", quiet = TRUE) %>% st_transform(5070)
crime_county <- st_read("", quiet = TRUE) %>% st_transform(5070)
all_crime <- dplyr::bind_rows(crime_city, crime_county)

crime_to_grid <- st_join(all_crime, grid_gwr, left = FALSE) %>%
  st_drop_geometry() %>%
  count(grid_id, name = "crime_count")

grid_gwr <- grid_gwr %>%
  left_join(crime_to_grid, by = "grid_id") %>%
  mutate(crime_count = tidyr::replace_na(crime_count, 0L))

# ----------------------------------------
# Step 3: Add exposure metrics (binary + continuous + distance)
# ----------------------------------------
# Binary fire_exposed (intersects original fire perimeters)
grid_gwr <- grid_gwr %>%
  mutate(fire_exposed = lengths(st_intersects(., st_union(fires))) > 0)

# Continuous: pct_damaged already merged earlier into la_grid; keep if present
if (!"pct_damaged" %in% names(grid_gwr)) {
  grid_gwr$pct_damaged <- 0
}

# Distance to nearest fire edge (meters)
grid_centroids <- st_centroid(grid_gwr)
# both in EPSG:5070
dist_to_fire <- st_distance(grid_centroids, st_union(fires))  # n x 1
grid_gwr$dist_fire_m <- as.numeric(units::set_units(dist_to_fire, "m"))

# ----------------------------------------
# Step 4: Areal-weighted population → grid (ACS 2020)
# ----------------------------------------
options(tigris_use_cache = TRUE)
census_api_key(Sys.getenv("CENSUS_API_KEY"), install = FALSE, overwrite = FALSE)

la_bg <- get_acs(
  geography = "block group",
  variables = "B01003_001",
  year = 2020,
  state = "CA",
  county = "Los Angeles",
  geometry = TRUE
) %>% st_transform(5070) %>% st_make_valid()

# Areal weighting (population proportion by overlap)
aw <- st_intersection(
  grid_gwr %>% select(grid_id),
  la_bg %>% select(GEOID, pop = estimate)
) %>%
  mutate(
    area_overlap = st_area(.),
    area_bg = st_area(la_bg[match(GEOID, la_bg$GEOID), ])
  ) %>%
  st_drop_geometry() %>%
  mutate(weight = as.numeric(area_overlap / area_bg)) %>%
  group_by(grid_id) %>%
  summarise(pop_aw = sum(pop * pmin(pmax(weight, 0, na.rm = TRUE), 1), na.rm = TRUE), .groups = "drop")

grid_gwr <- grid_gwr %>%
  left_join(aw, by = "grid_id") %>%
  mutate(pop_aw = tidyr::replace_na(pop_aw, 0),
         pop_density_km2 = ifelse(as.numeric(st_area(.)) > 0,
                                  pop_aw / (as.numeric(st_area(.)) / 1e6), 0))

# ----------------------------------------
# Step 5: GWR (prefer Poisson; fallback to Gaussian log-counts)
# ----------------------------------------
# Clean inputs, create scaled vars, and distance matrix
grid_sp <- as(grid_gwr, "Spatial")

# exposure as 0/1 numeric
grid_sp@data$fire_exposed <- as.integer(grid_sp@data$fire_exposed)

# scale: pct_damaged per 10 percentage points (e.g., 0.30 -> 3)
grid_sp@data$pct_damaged_10 <- grid_sp@data$pct_damaged * 10

# scale continuous covariates for numerical stability
grid_sp@data$pop_density_km2_s <- as.numeric(scale(grid_sp@data$pop_density_km2))
grid_sp@data$dist_fire_m_s     <- as.numeric(scale(grid_sp@data$dist_fire_m))

# response for fallback
grid_sp@data$log_crime <- log1p(grid_sp@data$crime_count)

# remove rows with any NA in variables used
keep_cols <- c("crime_count","log_crime","fire_exposed","pct_damaged_10",
               "pop_density_km2_s","dist_fire_m_s")
grid_sp <- grid_sp[stats::complete.cases(grid_sp@data[, keep_cols]), ]

# distance matrix (centroids; EPSG:5070 => planar)
coords <- sp::coordinates(grid_sp)
dMat   <- GWmodel::gw.dist(dp.locat = coords)

# Prefer Poisson GWR if available (ggwr & bw.ggwr); else fallback to Gaussian
use_poisson <- isTRUE("ggwr" %in% getNamespaceExports("GWmodel")) &&
  isTRUE("bw.ggwr" %in% getNamespaceExports("GWmodel"))

if (use_poisson) {
  # ----- Poisson GWR -----
  bw <- GWmodel::bw.ggwr(
    formula = crime_count ~ fire_exposed + pct_damaged_10 + pop_density_km2_s + dist_fire_m_s,
    data    = grid_sp,
    family  = "poisson",
    kernel  = "bisquare",
    adaptive = TRUE,
    dMat    = dMat
  )
  
  gwr_fit <- GWmodel::ggwr(
    formula = crime_count ~ fire_exposed + pct_damaged_10 + pop_density_km2_s + dist_fire_m_s,
    data    = grid_sp,
    bw      = bw,
    family  = "poisson",
    kernel  = "bisquare",
    adaptive = TRUE,
    dMat    = dMat
  )
  
  gwr_family <- "poisson"
} else {
  # ----- Gaussian GWR on log(counts) (fallback) -----
  bw <- GWmodel::bw.gwr(
    formula = log_crime ~ fire_exposed + pct_damaged_10 + pop_density_km2_s + dist_fire_m_s,
    data    = grid_sp,
    kernel  = "bisquare",
    adaptive = TRUE,
    dMat    = dMat
  )
  
  gwr_fit <- GWmodel::gwr.basic(
    formula = log_crime ~ fire_exposed + pct_damaged_10 + pop_density_km2_s + dist_fire_m_s,
    data    = grid_sp,
    bw      = bw,
    kernel  = "bisquare",
    adaptive = TRUE,
    dMat    = dMat
  )
  
  gwr_family <- "gaussian_log"
}

gwr_sf <- sf::st_as_sf(gwr_fit$SDF)]

# ----------------------------------------
# Step 5a: SAFE scaling + predictor diagnostics (before bandwidth search)
# ----------------------------------------

# 1) Safe scaler: returns zeros if sd == 0 (prevents NaN from scale())
safe_scale <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x[is.finite(x)], na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    rep(0, length(x))
  } else {
    as.numeric(scale(x))
  }
}

grid_sp <- as(grid_gwr, "Spatial")

grid_sp@data$fire_exposed <- as.integer(grid_sp@data$fire_exposed)

grid_sp@data$pct_damaged_10 <- grid_sp@data$pct_damaged * 10
grid_sp@data$pop_density_km2_s <- safe_scale(grid_sp@data$pop_density_km2)
grid_sp@data$dist_fire_m_s     <- safe_scale(grid_sp@data$dist_fire_m)

grid_sp@data$log_crime <- log1p(grid_sp@data$crime_count)

# 2) Drop zero-variance or all-NA predictors automatically
preds <- c("fire_exposed","pct_damaged_10","pop_density_km2_s","dist_fire_m_s")
keep <- vapply(preds, function(v) {
  x <- grid_sp@data[[v]]
  fin <- is.finite(x)
  if (!any(fin)) return(FALSE)
  # variance among finite values
  var_x <- stats::var(x[fin], na.rm = TRUE)
  is.finite(var_x) && var_x > 0
}, logical(1))

if (!all(keep)) {
  message("Dropping near-constant/all-NA predictors: ",
          paste(preds[!keep], collapse = ", "))
  preds <- preds[keep]
}

# 3) Complete cases only on the columns
keep_cols <- unique(c("crime_count","log_crime", preds))
grid_sp <- grid_sp[stats::complete.cases(grid_sp@data[, keep_cols, drop = FALSE]), ]

# 4) ensure all finite
stopifnot(all(sapply(keep_cols, function(v) all(is.finite(grid_sp@data[[v]])))))

# 5) Distance matrix
coords <- sp::coordinates(grid_sp)
dMat   <- GWmodel::gw.dist(dp.locat = coords)

# ----------------------------------------
# Step 5b: Robust bandwidth search with guarded Poisson; fallback if needed
# ----------------------------------------

use_poisson <- isTRUE("ggwr" %in% getNamespaceExports("GWmodel")) &&
  isTRUE("bw.ggwr" %in% getNamespaceExports("GWmodel"))

form_poisson <- as.formula(
  paste0("crime_count ~ ", paste(preds, collapse = " + "))
)
form_gauss <- as.formula(
  paste0("log_crime ~ ", paste(preds, collapse = " + "))
)

gwr_family <- NA_character_
gwr_fit <- NULL

if (use_poisson) {
  bw_try <- try({
    GWmodel::bw.ggwr(
      formula = form_poisson,
      data    = grid_sp,
      family  = "poisson",
      kernel  = "bisquare",
      adaptive = TRUE,
      dMat    = dMat
      # (optional) you can also add: approach = "AICc"  # if your GWmodel supports it
    )
  }, silent = TRUE)
  
  if (!inherits(bw_try, "try-error") && is.finite(bw_try) && bw_try > 0) {
    gwr_fit <- GWmodel::ggwr(
      formula = form_poisson,
      data    = grid_sp,
      bw      = bw_try,
      family  = "poisson",
      kernel  = "bisquare",
      adaptive = TRUE,
      dMat    = dMat
    )
    gwr_family <- "poisson"
  } else {
    message("Poisson GWR bandwidth search failed (CV=Inf/NA). Falling back to Gaussian log-counts.")
  }
}

# Fallback to Gaussian log-counts if Poisson failed
if (is.na(gwr_family)) {
  bw_try2 <- GWmodel::bw.gwr(
    formula = form_gauss,
    data    = grid_sp,
    kernel  = "bisquare",
    adaptive = TRUE,
    dMat    = dMat,
    approach = "AICc"
  )
  gwr_fit <- GWmodel::gwr.basic(
    formula = form_gauss,
    data    = grid_sp,
    bw      = bw_try2,
    kernel  = "bisquare",
    adaptive = TRUE,
    dMat    = dMat
  )
  gwr_family <- "gaussian_log"
}

gwr_sf <- sf::st_as_sf(gwr_fit$SDF)




