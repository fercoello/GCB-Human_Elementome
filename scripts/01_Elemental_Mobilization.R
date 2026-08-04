# ==============================================================================
# Project: GCB — Anthropogenic Effects on the Elementome
#
# Script:  01_Elemental_Mobilization.R
#
# Purpose:
#   Reconstruct and quantify the annual anthropogenic mobilization of elements
#   associated with mining, construction materials, coal consumption, and oil
#   consumption. The script also propagates uncertainty through Monte Carlo
#   simulations and generates the main summary tables and figures.
#
# Authors:
#   Fernando Coello Sanz
#
# Affiliations:
#   CREAF
#
# Contact:
#   f.coello@creaf.uab.cat
#
# Repository:
#   To complete
#
# Created:
#   2026-01-30
#
# Last updated:
#   2026-07-30
#
# R version:
#   R 4.5.1
#
# Main dependencies:
#   tidyverse, readxl, readr, here, zoo, purrr, patchwork
#
# Inputs:
#   data/data_FM1.xlsx
#   data/data_REO.xlsx
#   data/data_construction.csv
#   data/data_fossil_fuel.xlsx
#   data/elemental_concentration/data_elemental_concentration.xlsx
#
# Outputs:
#   results/el_mob_construction_avg.csv
#   results/el_mob_construction_draws.csv
#   results/mob_FM1.csv
#   results/anthropogenic_mobilization_summary.csv
#   figures/Figure_2.pdf
#
# Reproducibility:
#   Random-number seeds are set before each Monte Carlo analysis. Run the script
#   from the root directory of the RStudio project. File paths are constructed
#   relative to the project root using here::here().
#
# Notes:
#   Concentrations are standardized as mass fractions in kg kg-1 before being
#   multiplied by material production or consumption. Mobilized elemental masses
#   are expressed in Tg unless stated otherwise.
#
# License:
#   CC BY 4.0
#
# Citation:
#   [To update]
# ==============================================================================

# 00. Setup ----------------------------------------------------------------
# Purpose:
#   1. Load the packages required by the workflow.
#   2. define project paths;
#   3. create output directories;
#   4. declare all important constants and assumptions in one place.

# 1. Loading libraries
library(tidyverse)
library(readxl)
library(zoo)
library(patchwork)
library(here)

# Reproducibility
SEED <- 27L

# Number of Monte Carlo draws used for final flux calculations.
N_SIM <- 500L

# A larger number is used when constructing the average clay composition.
N_COMPOSITION_SIM <- 500L

# Normal-distribution quantiles used repeatedly.
Z_01 <- qnorm(0.01)
Z_25 <- qnorm(0.25)
Z_75 <- qnorm(0.75)
Z_99 <- qnorm(0.99)

# Project paths -------------------------------------------------------------
DATA_DIR <- here("data")
ELEMENT_CONCENTRATION_DIR <- here("data", "elemental_concentration")
RESULTS_DIR <- here("results")
TABULAR_RESULTS_DIR <- here("results", "tabular_results")
DRAW_RESULTS_DIR <- here("results", "draws")
FIGURE_DIR <- here("figures")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABULAR_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DRAW_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# Mass-conversion assumptions ----------------------------------------------
# IEA conversion factors from TWh to Tg
COAL_MASS_FACTOR <- 0.1228
OIL_MASS_FACTOR  <- 0.086

# Construction-material assumptions ----------------------------------------

# Proportion of carbonates extraction (Eggleston, 2006)
LIMESTONE_SHARE <- 0.85
DOLOMITE_SHARE  <- 0.15

# Global Proportion of different clays (Ito and Wagai (2017)
CLAY_WEIGHTS <- tibble(
  clay = c("kaolinite", "smectite", "illite"),
  weight_raw = c(27.26, 14.40, 35.50)
) %>%
  mutate(weight = weight_raw / sum(weight_raw)) %>%
  select(clay, weight)

# Minimum number of observations required to prefer the material-specific
# limestone or dolomite estimate over the UCC fallback.
MIN_CARBONATE_N <- 50L

# Iron-ore grade assumptions ------------------------------------------------
FE_GRADE_1900 <- 0.60
FE_GRADE_2018 <- 0.4481
FE_GRADE_START_YEAR <- 1900L
FE_GRADE_END_YEAR   <- 2018L

# REE assumptions -----------------------------------------------------------
REE_ELEMENTS <- c(
  "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb",
  "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Y"
)

# Row assignations for China minings operations
CHINA_BAYAN_OBO_ROWS <- 1L
CHINA_SICHUAN_ROWS   <- 2:3
CHINA_OTHER_ROWS     <- 4:9

CHINA_DEPOSIT_SHARES <- c(
  Bayan_Obo = 0.55,
  Sichuan   = 0.27,
  Others    = 0.18
)

# Values written as "<x" in the REE table are treated as x, reproducing the
# original workflow. Use 0.5 instead of 1 if half the reporting limit is
# preferred.
REE_LESS_THAN_MULTIPLIER <- 1


# 01. Helper functions ------------------------------------------------------
#
# This file contains small functions reused across several parts of the
# workflow. Keeping them here prevents slightly different versions of the
# same calculation from appearing in multiple sections.

# Restrict values to a specified interval.
clip <- function(x, lower = -Inf, upper = Inf) {
  pmin(pmax(x, lower), upper)
}

# Return NA when every value is missing; otherwise sum the observed values.
#
# This prevents sum(c(NA, NA), na.rm = TRUE) from being interpreted as a
# genuine zero.
sum_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
}

# Collapse unique non-missing text values into one traceable label.
paste_unique <- function(x, separator = "; ") {
  values <- unique(x[!is.na(x) & x != ""])
  if (length(values) == 0) NA_character_ else paste(values, collapse = separator)
}

# Convert numbers stored as text into numeric values.
#
# parse_number() also handles strings containing units or annotations.
to_numeric <- function(x) {
  readr::parse_number(as.character(x))
}

# Replace strings beginning with "<" by NA.
#
# Used when below-detection-limit measurements should not be used.
below_detection_to_na <- function(x) {
  x <- as.character(x)
  x[str_detect(x, "^\\s*<")] <- NA_character_
  x
}

# Parse REE composition values.
#
# "NA" becomes 0, while "<x" becomes x multiplied by the configurable
# reporting-limit multiplier.
parse_ree_value <- function(x, less_than_multiplier = REE_LESS_THAN_MULTIPLIER) {
  x_chr <- str_trim(as.character(x))
  
  case_when(
    is.na(x_chr) | x_chr %in% c("", "NA", "N/A") ~ 0,
    str_detect(x_chr, "^<") ~
      to_numeric(str_remove(x_chr, "^<")) * less_than_multiplier,
    TRUE ~ to_numeric(x_chr)
  )
}

# Convert composition units to a mass fraction in kg/kg.
to_mass_fraction <- function(value, unit) {
  unit_clean <- unit %>%
    as.character() %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all("μ", "µ")
  
  case_when(
    unit_clean %in% c(
      "%", "% w/w", "percent", "perc",
      "wt%", "wt.%", "weight %", "wt%_element"
    ) ~ value / 100,
    
    unit_clean %in% c("ppm", "mg/kg", "ug/g", "µg/g") ~ value * 1e-6,
    unit_clean == "ng/g" ~ value * 1e-9,
    unit_clean == "g/kg" ~ value * 1e-3,
    unit_clean %in% c("kg/kg", "fraction") ~ value,
    TRUE ~ NA_real_
  )
}

# Convert the 25th and 75th percentiles of a lognormal distribution into
# random draws.
draw_lognormal_from_iqr <- function(n, q25, q75, epsilon = 1e-15) {
  q25 <- pmax(q25, epsilon)
  q75 <- pmax(q75, epsilon)
  
  sigma <- (log(q75) - log(q25)) / (Z_75 - Z_25)
  mu <- (log(q25) + log(q75)) / 2
  
  rlnorm(n, meanlog = mu, sdlog = sigma)
}

# Summary statistics used for all Monte Carlo mobilization tables.
summarise_mc <- function(data, value_column = "mob_sim") {
  value_column <- rlang::ensym(value_column)
  
  data %>%
    summarise(
      mob_mean   = mean(!!value_column, na.rm = TRUE),
      mob_median = median(!!value_column, na.rm = TRUE),
      mob_sd     = sd(!!value_column, na.rm = TRUE),
      mob_q25    = quantile(!!value_column, 0.25, na.rm = TRUE),
      mob_q75    = quantile(!!value_column, 0.75, na.rm = TRUE),
      mob_IQR    = mob_q75 - mob_q25,
      .groups = "drop"
    )
}

# Check that an input table has all columns expected by the workflow.
check_required_columns <- function(data, required, object_name) {
  missing <- setdiff(required, names(data))
  
  if (length(missing) > 0) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
  
  invisible(data)
}

# Standardise coal-type names used in different source sheets.
standardise_coal_type <- function(x) {
  case_when(
    x == "Anthracite" ~ "Anthracite",
    x == "Metallurgical coal" ~ "Metallurgical",
    x == "Metallurgical" ~ "Metallurgical",
    x == "Bituminous" ~ "Bituminous",
    x == "Subbituminous" ~ "Subbituminous",
    x == "Lignite" ~ "Lignite",
    x == "perc_Anthracite" ~ "Anthracite",
    x == "perc_Metallurgical" ~ "Metallurgical",
    x == "perc_Bituminous" ~ "Bituminous",
    x == "perc_Subbituminous" ~ "Subbituminous",
    x == "perc_Lignite" ~ "Lignite",
    TRUE ~ x
  )
}

# Calculate the summary statistics of a uniform distribution.
uniform_stats <- function(min_value, max_value) {
  tibble(
    median = (min_value + max_value) / 2,
    sd = (max_value - min_value) / sqrt(12),
    q25 = min_value + 0.25 * (max_value - min_value),
    q75 = min_value + 0.75 * (max_value - min_value)
  )
}

# 02. Load and minimally clean input data ----------------------------------
#
# This section reads data and performs transformations that are intrinsic
# to the source table, such as interpolating the fossil-fuel consumption
# series. Material-composition calculations are carried out in later scripts.

## Mining and ore-production data -------------------------------------------

data_mining <- read_excel(
  "data/data_FM1.xlsx",
  sheet = "FM1_long_data",
  col_types = "text"
) %>%
  mutate(
    Year = parse_double(Year),
    
    Mining_Production = parse_number(
      Mining_Production,
      na = c("", "NA", "N/A")
    ),
    
    Note = str_squish(Note),
    Note = na_if(Note, ""),
    Note = na_if(Note, "NA"),
    
    Finish = str_squish(Finish)
  )

## Rare-earth production and composition ------------------------------------
data_REO_production <- read_excel(
  file.path(DATA_DIR, "data_REO.xlsx"),
  sheet = "REO_country_production"
)

data_REO_content <- read_excel(
  file.path(DATA_DIR, "data_REO.xlsx"),
  sheet = "REO_content"
)

## Construction-material mass series ----------------------------------------
data_construction <- read_csv(
  file.path(DATA_DIR, "data_construction.csv"),
  show_col_types = FALSE
)

## Fossil-fuel consumption --------------------------------------------------
fossil_consumption_raw <- read_excel(
  file.path(DATA_DIR, "data_fossil_fuel.xlsx"),
  sheet = "Fossil_Fuel_Consumption"
)

check_required_columns(
  fossil_consumption_raw,
  c("Product", "Year", "Value"),
  "Fossil_Fuel_Consumption"
)

data_coal <- fossil_consumption_raw %>%
  filter(Product == "coal") %>%
  arrange(Year) %>%
  mutate(
    Value_interp = zoo::na.approx(
      Value,
      x = Year,
      na.rm = FALSE
    ),
    # The source mass is converted to Tg using the factor declared in setup.
    Coal_Tg = Value_interp * COAL_MASS_FACTOR
  )

data_oil <- fossil_consumption_raw %>%
  filter(Product == "oil") %>%
  arrange(Year) %>%
  mutate(
    Value_interp = zoo::na.approx(
      Value,
      x = Year,
      na.rm = FALSE
    ),
    Oil_Tg = Value_interp * OIL_MASS_FACTOR
  )

## Construction-material elemental compositions -----------------------------
element_composition_file <- file.path(
  ELEMENT_CONCENTRATION_DIR,
  "data_elemental_concentration.xlsx"
)

el_UCC_raw <- read_excel(
  element_composition_file,
  sheet = "UCC"
)

el_sand_gravel_raw <- read_excel(
  element_composition_file,
  sheet = "river_sand_gravel"
)

el_limestone_raw <- read_excel(
  element_composition_file,
  sheet = "limestone_best"
) %>%
  rename(Element = 1)

el_dolomite_raw <- read_excel(
  element_composition_file,
  sheet = "dolomite_best"
) %>%
  rename(Element = 1)

el_clays_raw <- read_excel(
  element_composition_file,
  sheet = "clay_best"
)

## Coal composition inputs --------------------------------------------------
el_major_coal_raw <- read_excel(
  file.path(DATA_DIR, "data_fossil_fuel.xlsx"),
  sheet = "MajorEl_concentration_coal"
)

el_trace_coal_raw <- read_excel(
  file.path(DATA_DIR, "data_fossil_fuel.xlsx"),
  sheet = "TraceEl_concentration_coal"
)

coal_type_consumption_raw <- read_excel(
  file.path(DATA_DIR, "data_fossil_fuel.xlsx"),
  sheet = "Coal_Type_Consumption"
)

## Oil composition input ----------------------------------------------------
#
# The spelling of this sheet name reproduces the current Excel workbook.
el_oil_raw <- read_excel(
  file.path(DATA_DIR, "data_fossil_fuel.xlsx"),
  sheet = "Element_conentration_crude_oil"
)

## 03. Construction-material elemental composition -------------------------
#
# This section prepares elemental-composition tables for:
#   1. upper continental crust (UCC);
#   2. sand and gravel;
#   3. limestone and dolomite;
#   4. a weighted clay mixture.
#
# All final composition values are converted to kg element per kg material.

## 03.1 Upper continental crust ---------------------------------------------
#
if (nrow(el_UCC_raw) < 86) {
  stop("The UCC sheet contains fewer than 86 rows; check the input workbook.")
}

el_UCC <- el_UCC_raw %>%
  slice(11:86) %>%
  mutate(
    Upper_crust_concentration = as.numeric(Upper_crust_concentration),
    Upper_crust_concentration_sd = as.numeric(Upper_crust_concentration_sd),
    CV = Upper_crust_concentration_sd / Upper_crust_concentration
  )

# Use the median observed coefficient of variation to fill missing UCC SDs.
ucc_cv_reference <- median(el_UCC$CV, na.rm = TRUE)

el_UCC <- el_UCC %>%
  mutate(
    Upper_crust_concentration_sd = coalesce(
      Upper_crust_concentration_sd,
      Upper_crust_concentration * ucc_cv_reference
    )
  )

# UCC table in mass-fraction units, reused as a fallback in several sections.
el_UCC_fraction <- el_UCC %>%
  transmute(
    Element,
    med_frac = to_mass_fraction(Upper_crust_concentration, Unit),
    sd_frac = to_mass_fraction(Upper_crust_concentration_sd, Unit),
    q25_frac = clip(med_frac + Z_25 * sd_frac, lower = 0),
    q75_frac = clip(med_frac + Z_75 * sd_frac, lower = 0),
    unit_best = "kg/kg",
    source = "UCC_Rudnick"
  )

## 03.2 Oxide-to-element conversion factors ---------------------------------
#
# Factors are the mass of the target element divided by the molar mass of the
# oxide. They are applied to both the central estimate and its uncertainty.

oxide_to_element <- tribble(
  ~Oxide,   ~Element_out, ~factor,
  "SiO2",   "Si",          28.085 / 60.084,
  "Al2O3",  "Al",          (2 * 26.982) / 101.961,
  "Fe2O3T", "Fe",          (2 * 55.845) / 159.687,
  "Fe2O3",  "Fe",          (2 * 55.845) / 159.687,
  "FeO",    "Fe",          55.845 / 71.844,
  "MnO",    "Mn",          54.938 / 70.937,
  "CaO",    "Ca",          40.078 / 56.077,
  "MgO",    "Mg",          24.305 / 40.304,
  "K2O",    "K",           (2 * 39.098) / 94.196,
  "Na2O",   "Na",          (2 * 22.990) / 61.979,
  "TiO2",   "Ti",          47.867 / 79.866,
  "P2O5T",  "P",           (2 * 30.974) / 141.944,
  "P2O5",   "P",           (2 * 30.974) / 141.944
)

## 03.3 Sand and gravel -----------------------------------------------------
#
# Prefer river-sediment values where available. Missing elements are filled
# with UCC concentrations.

sand_gravel_source <- el_sand_gravel_raw %>%
  transmute(
    Element,
    mean_value = to_numeric(mean),
    sd_value = to_numeric(sd),
    unit = Units
  ) %>%
  left_join(
    oxide_to_element,
    by = c("Element" = "Oxide")
  ) %>%
  transmute(
    Element = coalesce(Element_out, Element),
    mean_value = if_else(
      !is.na(factor),
      mean_value * factor,
      mean_value
    ),
    sd_value = if_else(
      !is.na(factor),
      sd_value * factor,
      sd_value
    ),
    unit
  )

el_sand_gravel <- el_UCC %>%
  full_join(sand_gravel_source, by = "Element") %>%
  mutate(
    across(where(is.character), ~ na_if(.x, "NA")),
    best_value = coalesce(mean_value, Upper_crust_concentration),
    best_value_sd = coalesce(sd_value, Upper_crust_concentration_sd),
    unit_best = coalesce(unit, Unit),
    source = if_else(
      !is.na(mean_value),
      "Sediments_Muller",
      "UCC_Rudnick"
    ),
    CV = best_value_sd / best_value
  )

## Fill remaining missing SDs with the mean CV among available sediment/UCC
# values. This reproduces the original general approach.
sand_gravel_cv_reference <- mean(el_sand_gravel$CV, na.rm = TRUE)

el_sand_gravel <- el_sand_gravel %>%
  mutate(
    best_value_sd = coalesce(
      best_value_sd,
      best_value * sand_gravel_cv_reference
    ),
    
    # For a normal distribution, q25 and q75 are mean ± 0.67449 SD.
    best_value_q25 = clip(
      best_value + Z_25 * best_value_sd,
      lower = 0
    ),
    best_value_q75 = clip(
      best_value + Z_75 * best_value_sd,
      lower = 0
    ),
    
    med_frac = to_mass_fraction(best_value, unit_best),
    sd_frac = to_mass_fraction(best_value_sd, unit_best),
    q25_frac = to_mass_fraction(best_value_q25, unit_best),
    q75_frac = to_mass_fraction(best_value_q75, unit_best)
  ) %>%
  select(
    Element,
    best_value,
    best_value_sd,
    best_value_q25,
    best_value_q75,
    unit_best,
    source,
    med_frac,
    sd_frac,
    q25_frac,
    q75_frac
  ) %>%
  filter(Element != "LOI")

## 03.4 Limestone -----------------------------------------------------------
#
# Material-specific values are used only when more than MIN_CARBONATE_N
# observations are available. Otherwise UCC is used as a transparent fallback.

el_limestone <- el_UCC %>%
  full_join(el_limestone_raw, by = "Element") %>%
  mutate(
    across(where(is.character), ~ na_if(.x, "NA")),
    use_limestone = !is.na(n) & n > MIN_CARBONATE_N,
    
    best_value = if_else(
      use_limestone,
      as.numeric(median),
      Upper_crust_concentration
    ),
    unit_best = if_else(use_limestone, unit, Unit),
    
    best_value_q25 = if_else(
      use_limestone,
      as.numeric(q25),
      clip(
        Upper_crust_concentration +
          Z_25 * Upper_crust_concentration_sd,
        lower = 0
      )
    ),
    best_value_q75 = if_else(
      use_limestone,
      as.numeric(q75),
      clip(
        Upper_crust_concentration +
          Z_75 * Upper_crust_concentration_sd,
        lower = 0
      )
    ),
    
    source = if_else(
      use_limestone,
      "Limestone_Gard",
      "UCC_Rudnick"
    ),
    n_best = if_else(use_limestone, as.numeric(n), NA_real_)
  ) %>%
  select(
    Element,
    best_value,
    best_value_q25,
    best_value_q75,
    unit_best,
    source,
    n_best
  )

## 03.5 Dolomite ------------------------------------------------------------
#
# Priority order:
#   dolomite-specific value -> limestone value -> UCC value.

el_dolomite <- el_limestone %>%
  full_join(el_dolomite_raw, by = "Element") %>%
  mutate(
    across(where(is.character), ~ na_if(.x, "NA")),
    use_dolomite = !is.na(n) & n > MIN_CARBONATE_N,
    
    best_value = if_else(
      use_dolomite,
      as.numeric(median),
      best_value
    ),
    unit_best = if_else(use_dolomite, unit, unit_best),
    best_value_q25 = if_else(
      use_dolomite,
      as.numeric(q25),
      best_value_q25
    ),
    best_value_q75 = if_else(
      use_dolomite,
      as.numeric(q75),
      best_value_q75
    ),
    source = if_else(
      use_dolomite,
      "Dolomite_Gard",
      source
    ),
    n_best = if_else(
      use_dolomite,
      as.numeric(n),
      n_best
    )
  ) %>%
  select(
    Element,
    best_value,
    best_value_q25,
    best_value_q75,
    unit_best,
    source,
    n_best
  )

## 03.6 Clay-mineral source tables ------------------------------------------
#
# Post-Archean Australian Shale (PAAS) is used as the first fallback for
# mineral-specific clay compositions. UCC is the second fallback.

prepare_oxide_composition <- function(
    data,
    mean_column,
    sd_column,
    mean_name,
    sd_name,
    unit_name
) {
  mean_column <- rlang::ensym(mean_column)
  sd_column <- rlang::ensym(sd_column)
  
  data %>%
    transmute(
      Oxide = Element,
      mean_raw = as.numeric(!!mean_column),
      sd_raw = as.numeric(!!sd_column),
      unit_raw = Units
    ) %>%
    left_join(oxide_to_element, by = "Oxide") %>%
    transmute(
      Element = coalesce(Element_out, Oxide),
      mean_value = if_else(
        !is.na(factor),
        mean_raw * factor,
        mean_raw
      ),
      sd_value = if_else(
        !is.na(factor),
        sd_raw * factor,
        sd_raw
      ),
      unit_value = unit_raw
    ) %>%
    rename(
      !!mean_name := mean_value,
      !!sd_name := sd_value,
      !!unit_name := unit_value
    ) %>%
    filter(!is.na(.data[[mean_name]]))
}

el_PAAS <- prepare_oxide_composition(
  el_clays_raw,
  Shales,
  Shales_sd,
  "mean_shale",
  "sd_shale",
  "unit_shale"
)

el_kaolinite_source <- prepare_oxide_composition(
  el_clays_raw,
  Kaolinite1,
  Sd_Kaolinite1,
  "mean_mineral",
  "sd_mineral",
  "unit_mineral"
)

el_smectite_source <- prepare_oxide_composition(
  el_clays_raw,
  Na_montmorillonite_SW,
  Na_montmorillonite_sd,
  "mean_mineral",
  "sd_mineral",
  "unit_mineral"
)

el_illite_source <- prepare_oxide_composition(
  el_clays_raw,
  Illite,
  Illite_sd,
  "mean_mineral",
  "sd_mineral",
  "unit_mineral"
)

# Complete one clay mineral using the hierarchy:
# mineral-specific composition -> PAAS -> UCC.
complete_clay_mineral <- function(mineral_data, mineral_name) {
  output <- el_UCC %>%
    full_join(el_PAAS, by = "Element") %>%
    full_join(mineral_data, by = "Element") %>%
    mutate(
      source = case_when(
        !is.na(mean_mineral) ~ paste0("CMS_", mineral_name),
        !is.na(mean_shale) ~ "Shale_PAAS",
        !is.na(Upper_crust_concentration) ~ "UCC_Rudnick",
        TRUE ~ NA_character_
      ),
      best_value = coalesce(
        mean_mineral,
        mean_shale,
        Upper_crust_concentration
      ),
      best_value_sd = case_when(
        !is.na(mean_mineral) ~ sd_mineral,
        !is.na(mean_shale) ~ sd_shale,
        TRUE ~ Upper_crust_concentration_sd
      ),
      unit_best = case_when(
        !is.na(mean_mineral) ~ unit_mineral,
        !is.na(mean_shale) ~ unit_shale,
        TRUE ~ Unit
      )
    ) %>%
    select(
      Element,
      best_value,
      best_value_sd,
      unit_best,
      source
    ) %>%
    filter(!is.na(best_value))
  
  # Derive a scalar CV from elements with both a mean and SD.
  cv_reference <- output %>%
    filter(best_value > 0, !is.na(best_value_sd)) %>%
    summarise(
      CV = mean(best_value_sd / best_value, na.rm = TRUE)
    ) %>%
    pull(CV)
  
  output %>%
    mutate(
      best_value_sd = coalesce(
        best_value_sd,
        best_value * cv_reference
      ),
      med_frac = to_mass_fraction(best_value, unit_best),
      sd_frac = to_mass_fraction(best_value_sd, unit_best),
      q25_frac = clip(med_frac + Z_25 * sd_frac, lower = 0),
      q75_frac = clip(med_frac + Z_75 * sd_frac, lower = 0)
    ) %>%
    filter(!is.na(med_frac))
}

el_kaolinite <- complete_clay_mineral(
  el_kaolinite_source,
  "kaolinite"
)

el_smectite <- complete_clay_mineral(
  el_smectite_source,
  "smectite"
)

el_illite <- complete_clay_mineral(
  el_illite_source,
  "illite"
)

## 03.7 Weighted average clay composition -----------------------------------
#
# The three clay-mineral distributions are sampled independently. Each draw is
# multiplied by the corresponding mineral weight and then summed.

clay_composition_long <- bind_rows(
  el_kaolinite %>%
    transmute(
      Element,
      clay = "kaolinite",
      med_frac,
      sd_frac
    ),
  el_smectite %>%
    transmute(
      Element,
      clay = "smectite",
      med_frac,
      sd_frac
    ),
  el_illite %>%
    transmute(
      Element,
      clay = "illite",
      med_frac,
      sd_frac
    )
) %>%
  left_join(CLAY_WEIGHTS, by = "clay")

set.seed(SEED)

el_clay <- clay_composition_long %>%
  group_by(Element) %>%
  group_modify(
    ~ {
      group_data <- .x %>%
        filter(!is.na(med_frac), !is.na(weight)) %>%
        mutate(
          # Renormalisation matters only if one component is unexpectedly
          # missing for an element.
          weight = weight / sum(weight)
        )
      
      component_draws <- map2(
        group_data$med_frac,
        group_data$sd_frac,
        ~ {
          sd_use <- if_else(
            is.na(.y) | .y < 0,
            0,
            .y
          )
          
          rnorm(
            N_COMPOSITION_SIM,
            mean = .x,
            sd = sd_use
          ) %>%
            clip(lower = 0, upper = 1)
        }
      )
      
      mixture_draws <- map2(
        component_draws,
        group_data$weight,
        ~ .x * .y
      ) %>%
        reduce(`+`) %>%
        clip(lower = 0, upper = 1)
      
      tibble(
        med_frac = median(mixture_draws, na.rm = TRUE),
        sd_frac = sd(mixture_draws, na.rm = TRUE),
        q25_frac = quantile(mixture_draws, 0.25, na.rm = TRUE),
        q75_frac = quantile(mixture_draws, 0.75, na.rm = TRUE)
      )
    }
  ) %>%
  ungroup() %>%
  filter(Element != "LOI")

# 04. Fossil-fuel elemental composition -----------------------------------
#
# This script prepares:
#   1. a coal composition table;
#   2. annual Monte Carlo estimates of coal C and S composition;
#   3. an oil composition table.
#
# Final compositions are stored as kg element per kg fuel.

## 04.1 Coal-type shares through time ---------------------------------------

coal_type_shares <- coal_type_consumption_raw %>%
  rename(
    Anthracite = 3,
    Metallurgical = 4,
    Bituminous = 5,
    Subbituminous = 6,
    Lignite = 7
  ) %>%
  mutate(
    perc_Anthracite = Anthracite / Tota_Coal_Mst,
    perc_Metallurgical = Metallurgical / Tota_Coal_Mst,
    perc_Bituminous = Bituminous / Tota_Coal_Mst,
    perc_Subbituminous = Subbituminous / Tota_Coal_Mst,
    perc_Lignite = Lignite / Tota_Coal_Mst,
    share_sum = perc_Anthracite +
      perc_Metallurgical +
      perc_Bituminous +
      perc_Subbituminous +
      perc_Lignite
  ) %>%
  select(
    Year,
    perc_Anthracite,
    perc_Metallurgical,
    perc_Bituminous,
    perc_Subbituminous,
    perc_Lignite,
    share_sum
  )

# Keep the share-sum diagnostic. Values should normally be close to one.
coal_share_diagnostic <- coal_type_shares %>%
  summarise(
    min_share_sum = min(share_sum, na.rm = TRUE),
    max_share_sum = max(share_sum, na.rm = TRUE)
  )

## 04.2 Coal major-element table --------------------------------------------

el_major_coal <- el_major_coal_raw %>%
  rename(Coal_type = 1) %>%
  mutate(
    Coal_type = standardise_coal_type(Coal_type),
    Mean = coalesce(
      as.numeric(Unique),
      (as.numeric(Min) + as.numeric(Max)) / 2
    )
  ) %>%
  select(Coal_type, Element, Min, Mean, Max)

## 04.3 Coal trace-element composition --------------------------------------
#
# Source priority reproduces the original workflow:
#   - Si: Finkelman;
#   - all other elements: IEA where available, then Finkelman;
#   - missing values: UCC fallback.
#
# IEA minima and maxima are interpreted as the 1st and 99th percentiles of an
# approximately normal distribution.

el_coal_trace <- el_trace_coal_raw %>%
  mutate(
    Element = str_trim(Element),
    
    Crude_coal_IEA_clean =
      below_detection_to_na(Crude_coal_IEA),
    Crude_coal_min_IEA_clean =
      below_detection_to_na(Crude_coal_min_IEA),
    Crude_coal_max_IEA_clean =
      below_detection_to_na(Crude_coal_max_IEA),
    
    Coal_Finkelman_mean_num =
      to_numeric(Coal_Finkelman_mean),
    Coal_Finkelman_sd_num =
      to_numeric(Coal_Finkelman_sd),
    
    # In the source table, IEA zeros are interpreted as missing.
    Crude_coal_IEA_num =
      na_if(to_numeric(Crude_coal_IEA_clean), 0),
    Crude_coal_min_IEA_num =
      na_if(to_numeric(Crude_coal_min_IEA_clean), 0),
    Crude_coal_max_IEA_num =
      na_if(to_numeric(Crude_coal_max_IEA_clean), 0)
  ) %>%
  mutate(
    IEA_mean_num = case_when(
      !is.na(Crude_coal_min_IEA_num) &
        !is.na(Crude_coal_max_IEA_num) ~
        (Crude_coal_min_IEA_num +
           Crude_coal_max_IEA_num) / 2,
      
      TRUE ~ Crude_coal_IEA_num
    ),
    
    IEA_sd_num = case_when(
      !is.na(Crude_coal_min_IEA_num) &
        !is.na(Crude_coal_max_IEA_num) ~
        (Crude_coal_max_IEA_num -
           Crude_coal_min_IEA_num) /
        (Z_99 - Z_01),
      
      TRUE ~ NA_real_
    ),
    
    IEA_q25_num = IEA_mean_num + Z_25 * IEA_sd_num,
    IEA_q75_num = IEA_mean_num + Z_75 * IEA_sd_num
  ) %>%
  mutate(
    best_value = case_when(
      Element == "Si" &
        !is.na(Coal_Finkelman_mean_num) ~
        Coal_Finkelman_mean_num,
      
      !is.na(IEA_mean_num) ~ IEA_mean_num,
      !is.na(Coal_Finkelman_mean_num) ~ Coal_Finkelman_mean_num,
      TRUE ~ NA_real_
    ),
    
    best_value_q25 = case_when(
      Element == "Si" &
        !is.na(Coal_Finkelman_mean_num) &
        !is.na(Coal_Finkelman_sd_num) ~
        Coal_Finkelman_mean_num +
        Z_25 * Coal_Finkelman_sd_num,
      
      !is.na(IEA_q25_num) ~ IEA_q25_num,
      
      !is.na(Coal_Finkelman_mean_num) &
        !is.na(Coal_Finkelman_sd_num) ~
        Coal_Finkelman_mean_num +
        Z_25 * Coal_Finkelman_sd_num,
      
      TRUE ~ NA_real_
    ),
    
    best_value_q75 = case_when(
      Element == "Si" &
        !is.na(Coal_Finkelman_mean_num) &
        !is.na(Coal_Finkelman_sd_num) ~
        Coal_Finkelman_mean_num +
        Z_75 * Coal_Finkelman_sd_num,
      
      !is.na(IEA_q75_num) ~ IEA_q75_num,
      
      !is.na(Coal_Finkelman_mean_num) &
        !is.na(Coal_Finkelman_sd_num) ~
        Coal_Finkelman_mean_num +
        Z_75 * Coal_Finkelman_sd_num,
      
      TRUE ~ NA_real_
    ),
    
    source = case_when(
      Element == "Si" &
        !is.na(Coal_Finkelman_mean_num) ~
        "Coal_Finkelman",
      
      !is.na(IEA_mean_num) ~ "Coal_IEA",
      !is.na(Coal_Finkelman_mean_num) ~ "Coal_Finkelman",
      TRUE ~ NA_character_
    ),
    
    unit_best = Units
  ) %>%
  mutate(
    med_frac = to_mass_fraction(best_value, unit_best),
    q25_frac = clip(
      to_mass_fraction(best_value_q25, unit_best),
      lower = 0
    ),
    q75_frac = clip(
      to_mass_fraction(best_value_q75, unit_best),
      lower = 0
    ),
    sd_frac = clip(
      (q75_frac - q25_frac) / (Z_75 - Z_25),
      lower = 0
    ),
    unit_best = "kg/kg"
  ) %>%
  select(
    Element,
    med_frac,
    sd_frac,
    q25_frac,
    q75_frac,
    unit_best,
    source
  )

## Add UCC values for trace elements missing from the coal sources.
el_coal_with_ucc <- el_UCC_fraction %>%
  rename(
    ucc_med_frac = med_frac,
    ucc_sd_frac = sd_frac,
    ucc_q25_frac = q25_frac,
    ucc_q75_frac = q75_frac,
    ucc_source = source
  ) %>%
  full_join(el_coal_trace, by = "Element") %>%
  transmute(
    Element,
    med_frac = coalesce(med_frac, ucc_med_frac),
    sd_frac = coalesce(sd_frac, ucc_sd_frac),
    q25_frac = coalesce(q25_frac, ucc_q25_frac),
    q75_frac = coalesce(q75_frac, ucc_q75_frac),
    unit_best = "kg/kg",
    source = coalesce(source, ucc_source)
  ) %>%
  arrange(Element)

## 04.4 Annual C and S composition of coal ----------------------------------
#
# For each coal type, the reported minimum and maximum are treated as the
# limits of a uniform distribution. Coal-type draws are weighted by the annual
# shares of anthracite, metallurgical, bituminous, subbituminous, and lignite.

simulate_coal_major_element <- function(element_symbol) {
  parameters <- el_major_coal %>%
    filter(Element == element_symbol) %>%
    transmute(
      Coal_type = standardise_coal_type(Coal_type),
      Min = as.numeric(Min),
      Max = as.numeric(Max)
    )
  
  set.seed(SEED)
  
  composition_draws <- parameters %>%
    mutate(
      draws = map2(
        Min,
        Max,
        ~ runif(
          N_COMPOSITION_SIM,
          min = .x,
          max = .y
        )
      )
    ) %>%
    select(Coal_type, draws) %>%
    unnest_longer(draws, indices_to = "sim") %>%
    rename(content = draws)
  
  shares_long <- coal_type_shares %>%
    select(-share_sum) %>%
    pivot_longer(
      cols = starts_with("perc_"),
      names_to = "Coal_type",
      values_to = "share"
    ) %>%
    mutate(
      Coal_type = standardise_coal_type(Coal_type)
    )
  
  annual_draws <- shares_long %>%
    inner_join(composition_draws, by = "Coal_type") %>%
    mutate(weighted_content = share * content) %>%
    group_by(Year, sim) %>%
    summarise(
      annual_content = sum_or_na(weighted_content),
      .groups = "drop"
    )
  
  annual_summary <- annual_draws %>%
    group_by(Year) %>%
    summarise(
      mean = mean(annual_content, na.rm = TRUE),
      median = median(annual_content, na.rm = TRUE),
      sd = sd(annual_content, na.rm = TRUE),
      q2.5 = quantile(annual_content, 0.025, na.rm = TRUE),
      q25 = quantile(annual_content, 0.25, na.rm = TRUE),
      q75 = quantile(annual_content, 0.75, na.rm = TRUE),
      q97.5 = quantile(annual_content, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Element = element_symbol)
  
  list(
    draws = annual_draws,
    summary = annual_summary
  )
}

coal_C_mc <- simulate_coal_major_element("C")
coal_S_mc <- simulate_coal_major_element("S")

coal_major_yearly_summary <- bind_rows(
  coal_C_mc$summary,
  coal_S_mc$summary
)

# The original workflow averages the annual composition summaries to obtain
# one constant C fraction and one constant S fraction for the complete coal
# series. Keep this choice explicit because using annual values directly would
# produce a different mobilization series.
coal_major_average <- coal_major_yearly_summary %>%
  group_by(Element) %>%
  summarise(
    med_frac = mean(mean, na.rm = TRUE) / 100,
    q25_frac = mean(q25, na.rm = TRUE) / 100,
    q75_frac = mean(q75, na.rm = TRUE) / 100,
    sd_frac = (q75_frac - q25_frac) / (Z_75 - Z_25),
    source = "Coal_MC_average",
    unit_best = "kg/kg",
    .groups = "drop"
  )

el_coal_final <- el_coal_with_ucc %>%
  filter(!Element %in% c("C", "S")) %>%
  bind_rows(coal_major_average) %>%
  arrange(Element)

# 04.5 Oil composition -----------------------------------------------------
#
# Priority logic:
#   1. typical minimum and maximum -> uniform distribution;
#   2. multiple reported values -> uniform range spanning all values;
#   3. one reported value -> normal distribution using a reference CV;
#   4. no oil value -> UCC fallback.

oil_composition_raw <- el_oil_raw %>%
  rowwise() %>%
  mutate(
    summary = list({
      has_typical <-
        !is.na(Crude_oil_min_typical) &&
        !is.na(Crude_oil_max_typical)
      
      has_yang <- !is.na(Crude_oil_Yang)
      has_sen <- !is.na(Crude_oil_Sen)
      has_samples <-
        !is.na(Crude_oil_min_samples) ||
        !is.na(Crude_oil_max_samples)
      
      all_values <- c(
        Crude_oil_Sen,
        Crude_oil_Yang,
        Crude_oil_min_samples,
        Crude_oil_max_samples
      )
      
      observed_values <- all_values[!is.na(all_values)]
      
      only_yang <-
        has_yang &&
        length(observed_values) == 1 &&
        isTRUE(
          all.equal(
            as.numeric(observed_values),
            as.numeric(Crude_oil_Yang)
          )
        )
      
      only_sen <-
        has_sen &&
        length(observed_values) == 1 &&
        isTRUE(
          all.equal(
            as.numeric(observed_values),
            as.numeric(Crude_oil_Sen)
          )
        )
      
      if (has_typical) {
        stats <- uniform_stats(
          Crude_oil_min_typical,
          Crude_oil_max_typical
        )
        
        stats %>%
          mutate(
            source = "Oil_typical",
            distribution = "uniform",
            range_min = Crude_oil_min_typical,
            range_max = Crude_oil_max_typical
          )
        
      } else if (only_yang) {
        tibble(
          median = Crude_oil_Yang,
          sd = NA_real_,
          q25 = NA_real_,
          q75 = NA_real_,
          source = "Oil_Yang",
          distribution = "normal_cv",
          range_min = NA_real_,
          range_max = NA_real_
        )
        
      } else if (only_sen) {
        tibble(
          median = Crude_oil_Sen,
          sd = NA_real_,
          q25 = NA_real_,
          q75 = NA_real_,
          source = "Oil_Sen",
          distribution = "normal_cv",
          range_min = NA_real_,
          range_max = NA_real_
        )
        
      } else if (length(observed_values) >= 2) {
        range_min <- min(observed_values)
        range_max <- max(observed_values)
        stats <- uniform_stats(range_min, range_max)
        
        source_label <- case_when(
          has_sen & has_yang & has_samples ~
            "Oil_Sen+Yang+samples",
          has_sen & has_yang ~
            "Oil_Sen+Yang",
          has_sen & has_samples ~
            "Oil_Sen+samples",
          has_yang & has_samples ~
            "Oil_Yang+samples",
          has_sen ~
            "Oil_Sen_range",
          has_yang ~
            "Oil_Yang_range",
          TRUE ~
            "Oil_samples"
        )
        
        stats %>%
          mutate(
            source = source_label,
            distribution = "uniform",
            range_min = range_min,
            range_max = range_max
          )
        
      } else if (length(observed_values) == 1) {
        tibble(
          median = observed_values,
          sd = NA_real_,
          q25 = NA_real_,
          q75 = NA_real_,
          source = case_when(
            has_sen ~ "Oil_Sen",
            has_yang ~ "Oil_Yang",
            TRUE ~ "Oil_single"
          ),
          distribution = "normal_cv",
          range_min = NA_real_,
          range_max = NA_real_
        )
        
      } else {
        tibble(
          median = NA_real_,
          sd = NA_real_,
          q25 = NA_real_,
          q75 = NA_real_,
          source = NA_character_,
          distribution = NA_character_,
          range_min = NA_real_,
          range_max = NA_real_
        )
      }
    })
  ) %>%
  unnest(summary) %>%
  ungroup() %>%
  transmute(
    Element,
    med_raw = median,
    sd_raw = sd,
    q25_raw = q25,
    q75_raw = q75,
    range_min_raw = range_min,
    range_max_raw = range_max,
    unit_raw = "mg/kg",
    source,
    distribution
  )

# Reference CV from trace elements with uniform ranges.
oil_major_elements <- c("C", "H", "O", "N", "S", "Cl")

oil_trace_cv_reference <- oil_composition_raw %>%
  filter(
    !Element %in% oil_major_elements,
    distribution == "uniform",
    !is.na(med_raw),
    !is.na(sd_raw),
    med_raw > 0
  ) %>%
  summarise(
    CV = median(sd_raw / med_raw, na.rm = TRUE)
  ) %>%
  pull(CV)

# Add uncertainty to single-value oil observations.
oil_composition_completed <- oil_composition_raw %>%
  mutate(
    sd_raw = if_else(
      distribution == "normal_cv" & !is.na(med_raw),
      med_raw * oil_trace_cv_reference,
      sd_raw
    ),
    q25_raw = if_else(
      distribution == "normal_cv" & !is.na(med_raw),
      clip(
        med_raw + Z_25 * med_raw * oil_trace_cv_reference,
        lower = 0
      ),
      q25_raw
    ),
    q75_raw = if_else(
      distribution == "normal_cv" & !is.na(med_raw),
      clip(
        med_raw + Z_75 * med_raw * oil_trace_cv_reference,
        lower = 0
      ),
      q75_raw
    ),
    source = case_when(
      source == "Oil_Yang" ~ "Oil_Yang_CVtrace",
      source == "Oil_Sen" ~ "Oil_Sen_CVtrace",
      source == "Oil_single" ~ "Oil_single_CVtrace",
      TRUE ~ source
    )
  )

# Convert oil values to kg/kg and fill missing elements with UCC.
el_oil_final <- oil_composition_completed %>%
  mutate(
    from_oil = !is.na(med_raw),
    med_frac_oil = to_mass_fraction(med_raw, unit_raw),
    sd_frac_oil = to_mass_fraction(sd_raw, unit_raw),
    q25_frac_oil = clip(
      to_mass_fraction(q25_raw, unit_raw),
      lower = 0
    ),
    q75_frac_oil = clip(
      to_mass_fraction(q75_raw, unit_raw),
      lower = 0
    ),
    range_min_frac_oil = if_else(
      !is.na(range_min_raw),
      clip(
        to_mass_fraction(range_min_raw, unit_raw),
        lower = 0
      ),
      NA_real_
    ),
    range_max_frac_oil = if_else(
      !is.na(range_max_raw),
      clip(
        to_mass_fraction(range_max_raw, unit_raw),
        lower = 0
      ),
      NA_real_
    )
  ) %>%
  full_join(
    el_UCC_fraction %>%
      rename(
        med_frac_ucc = med_frac,
        sd_frac_ucc = sd_frac,
        q25_frac_ucc = q25_frac,
        q75_frac_ucc = q75_frac,
        source_ucc = source
      ),
    by = "Element"
  ) %>%
  transmute(
    Element,
    med_frac = coalesce(med_frac_oil, med_frac_ucc),
    sd_frac = coalesce(sd_frac_oil, sd_frac_ucc),
    q25_frac = coalesce(q25_frac_oil, q25_frac_ucc),
    q75_frac = coalesce(q75_frac_oil, q75_frac_ucc),
    range_min_frac = range_min_frac_oil,
    range_max_frac = range_max_frac_oil,
    unit_best = "kg/kg",
    source = case_when(
      coalesce(from_oil, FALSE) ~ source,
      !is.na(source_ucc) ~ source_ucc,
      TRUE ~ source
    ),
    distribution = case_when(
      coalesce(from_oil, FALSE) ~ distribution,
      !is.na(med_frac_ucc) ~ "normal_ucc",
      TRUE ~ distribution
    )
  ) %>%
  arrange(Element)

# 05. Construction-related elemental mobilization -------------------------
#
# General calculation:
#
#   elemental mobilization =
#     construction-material mass × elemental mass fraction
#
# Both material mass and composition are sampled when uncertainty estimates
# are available. All mobilization values are assumed to remain in Tg because
# material mass is supplied in Tg and composition is expressed as kg/kg.

## 05.1 Sand and gravel -----------------------------------------------------

aggregate_mass <- data_construction %>%
  transmute(
    Year,
    mass_mean = aggregates_best,
    mass_sd = aggregates_sd
  )

set.seed(SEED)

aggregate_mass_draws <- aggregate_mass %>%
  crossing(sim = seq_len(N_SIM)) %>%
  mutate(
    mass_sim = rnorm(
      n(),
      mean = mass_mean,
      sd = mass_sd
    ),
    mass_sim = clip(mass_sim, lower = 0)
  ) %>%
  select(Year, sim, mass_sim)

aggregate_composition_draws <- el_sand_gravel %>%
  select(Element, med_frac, sd_frac) %>%
  crossing(sim = seq_len(N_SIM)) %>%
  group_by(Element) %>%
  mutate(
    sd_use = if_else(
      is.na(sd_frac) | sd_frac < 0,
      0,
      sd_frac
    ),
    comp_sim = rnorm(
      n(),
      mean = med_frac[1],
      sd = sd_use[1]
    ),
    comp_sim = clip(comp_sim, lower = 0, upper = 1)
  ) %>%
  ungroup() %>%
  select(Element, sim, comp_sim)

mob_aggregate_sims <- aggregate_mass_draws %>%
  inner_join(aggregate_composition_draws, by = "sim") %>%
  transmute(
    Year,
    Element,
    sim,
    mob_sim = mass_sim * comp_sim,
    material = "sand and gravel"
  )

elem_mob_aggregate <- mob_aggregate_sims %>%
  group_by(Year, Element) %>%
  summarise_mc() %>%
  mutate(
    unit = "Tg",
    material = "sand and gravel"
  )

## 05.2 Limestone and dolomite mixture -------------------------------------
#
# The construction carbonate stream is represented by a weighted mixture of
# limestone and dolomite. Medians and quartiles are mixed using fixed material
# shares declared in 00_setup.R.

el_carbonate_mix <- el_limestone %>%
  select(
    Element,
    unit_limestone = unit_best,
    med_limestone = best_value,
    q25_limestone = best_value_q25,
    q75_limestone = best_value_q75
  ) %>%
  full_join(
    el_dolomite %>%
      select(
        Element,
        unit_dolomite = unit_best,
        med_dolomite = best_value,
        q25_dolomite = best_value_q25,
        q75_dolomite = best_value_q75
      ),
    by = "Element"
  ) %>%
  filter(
    !is.na(med_limestone),
    !is.na(q25_limestone),
    !is.na(q75_limestone),
    !is.na(med_dolomite),
    !is.na(q25_dolomite),
    !is.na(q75_dolomite)
  ) %>%
  mutate(
    # The units are expected to agree. This column remains available for
    # diagnostic checks.
    units_match =
      unit_limestone == unit_dolomite |
      (is.na(unit_limestone) & is.na(unit_dolomite)),
    
    unit = coalesce(unit_limestone, unit_dolomite),
    
    med = LIMESTONE_SHARE * med_limestone +
      DOLOMITE_SHARE * med_dolomite,
    
    q25 = LIMESTONE_SHARE * q25_limestone +
      DOLOMITE_SHARE * q25_dolomite,
    
    q75 = LIMESTONE_SHARE * q75_limestone +
      DOLOMITE_SHARE * q75_dolomite,
    
    med_frac = to_mass_fraction(med, unit),
    q25_frac = to_mass_fraction(q25, unit),
    q75_frac = to_mass_fraction(q75, unit)
  ) %>%
  select(
    Element,
    unit,
    units_match,
    med,
    q25,
    q75,
    med_frac,
    q25_frac,
    q75_frac
  )

carbonate_mass <- data_construction %>%
  transmute(
    Year,
    mass_mean = limestone_best,
    mass_sd = limestone_sd
  )

set.seed(SEED)

carbonate_mass_draws <- carbonate_mass %>%
  crossing(sim = seq_len(N_SIM)) %>%
  mutate(
    mass_sim = rnorm(
      n(),
      mean = mass_mean,
      sd = mass_sd
    ),
    mass_sim = clip(mass_sim, lower = 0)
  ) %>%
  select(Year, sim, mass_sim)

# A lognormal distribution is used because carbonate concentrations may be
# right-skewed. The distribution is parameterised by q25 and q75.
carbonate_composition_draws <- el_carbonate_mix %>%
  select(Element, q25_frac, q75_frac) %>%
  group_by(Element) %>%
  group_modify(
    ~ tibble(
      sim = seq_len(N_SIM),
      comp_sim = draw_lognormal_from_iqr(
        N_SIM,
        q25 = .x$q25_frac[1],
        q75 = .x$q75_frac[1]
      ) %>%
        clip(lower = 0, upper = 1)
    )
  ) %>%
  ungroup()

mob_carbonate_sims <- carbonate_mass_draws %>%
  inner_join(carbonate_composition_draws, by = "sim") %>%
  transmute(
    Year,
    Element,
    sim,
    mob_sim = mass_sim * comp_sim,
    material = "limestone and dolomites"
  )

elem_mob_carbonate <- mob_carbonate_sims %>%
  group_by(Year, Element) %>%
  summarise_mc() %>%
  mutate(
    unit = "Tg",
    material = "limestone and dolomites"
  )

## 05.3 Clay ----------------------------------------------------------------

clay_mass <- data_construction %>%
  transmute(
    Year,
    mass_mean = clay_best,
    mass_sd = clay_sd
  )

set.seed(SEED)

clay_mass_draws <- clay_mass %>%
  crossing(sim = seq_len(N_SIM)) %>%
  mutate(
    mass_sim = rnorm(
      n(),
      mean = mass_mean,
      sd = mass_sd
    ),
    mass_sim = clip(mass_sim, lower = 0)
  ) %>%
  select(Year, sim, mass_sim)

clay_composition_draws <- el_clay %>%
  select(Element, med_frac, sd_frac) %>%
  crossing(sim = seq_len(N_SIM)) %>%
  group_by(Element) %>%
  mutate(
    sd_use = if_else(
      is.na(sd_frac) | sd_frac < 0,
      0,
      sd_frac
    ),
    comp_sim = rnorm(
      n(),
      mean = med_frac[1],
      sd = sd_use[1]
    ),
    comp_sim = clip(comp_sim, lower = 0, upper = 1)
  ) %>%
  ungroup() %>%
  select(Element, sim, comp_sim)

mob_clay_sims <- clay_mass_draws %>%
  inner_join(clay_composition_draws, by = "sim") %>%
  transmute(
    Year,
    Element,
    sim,
    mob_sim = mass_sim * comp_sim,
    material = "clay"
  )

elem_mob_clay <- mob_clay_sims %>%
  group_by(Year, Element) %>%
  summarise_mc() %>%
  mutate(
    unit = "Tg",
    material = "clay"
  )

## 05.4 Total construction mobilization ------------------------------------
#
# Draws with the same simulation number are added so uncertainty propagates
# through the sum of the three construction-material streams.

mob_construction_component_sims <- bind_rows(
  mob_aggregate_sims,
  mob_carbonate_sims,
  mob_clay_sims
)

mob_construction_sims <- mob_construction_component_sims %>%
  group_by(Year, Element, sim) %>%
  summarise(
    mob_sim = sum_or_na(mob_sim),
    .groups = "drop"
  )

mob_construction_total <- mob_construction_sims %>%
  group_by(Year, Element) %>%
  summarise_mc()

# Save both the compact summary and the full Monte Carlo draws.
write_csv(
  mob_construction_total,
  file.path(TABULAR_RESULTS_DIR, "element_mobilization_construction_summary.csv")
)

write_csv(
  mob_construction_sims,
  file.path(DRAW_RESULTS_DIR, "element_mobilization_construction_draws.csv")
)

# 06. Mining-derived elemental mobilization (FM1) --------------------------
#
# The mining database contains a mixture of:
#   - commodities already expressed as elemental mass;
#   - compounds that require stoichiometric conversion;
#   - combined commodity groups;
#   - time-dependent ore-grade conversions;
#   - REO and PGM totals that require allocation among elements.
#
# This section converts all supported commodities into elemental FM1 values.

## 06.1 General cleaning ----------------------------------------------------

mining_clean <- data_mining %>%
  mutate(
    across(
      where(is.character),
      ~ na_if(str_trim(.x), "NA")
    ),
    Year = as.numeric(Year),
    Mining_Production = to_numeric(Mining_Production),
    Note = na_if(str_trim(Note), "")
  ) %>%
  rename(Commodity = Element)

check_required_columns(
  mining_clean,
  c(
    "Commodity",
    "Year",
    "Mining_Production",
    "Note",
    "Finish"
  ),
  "FM1_long_data"
)

## 06.2 Direct and fixed-factor conversions --------------------------------
#
# factor = kg target element per kg reported commodity.
#
# A factor of one means that the source already reports elemental production.

mining_conversion_table <- tribble(
  ~Commodity,     ~Element, ~factor,   ~conversion_note,
  "Ag",           "Ag",     1,         "Already elemental",
  "Al",           "Al",     1,         "Already elemental",
  "As",           "As",     1,         "Already elemental",
  "Au",           "Au",     1,         "Already elemental",
  "B2O3",         "B",      0.31058,   "B2O3 to B",
  "BaSO4",        "Ba",     0.58835,   "BaSO4 to Ba",
  "Be",           "Be",     1,         "Already elemental",
  "Bi",           "Bi",     1,         "Already elemental",
  "Br",           "Br",     1,         "Already elemental",
  "Cd",           "Cd",     1,         "Already elemental",
  "Co",           "Co",     1,         "Already elemental",
  "Cr",           "Cr",     1,         "Already elemental",
  "Cs",           "Cs",     1,         "Already elemental",
  "Cu",           "Cu",     1,         "Already elemental",
  "F",            "F",      0.48667,   "CaF2 to F; confirm source commodity",
  "Ga",           "Ga",     1,         "Already elemental",
  "Ge",           "Ge",     1,         "Already elemental",
  "He",           "He",     1,         "Already elemental",
  "Hf",           "Hf",     0.01531,   "Zirconium concentrate to Hf",
  "In",           "In",     1,         "Already elemental",
  "K",            "K",      0.83,      "K2O to K",
  "Li",           "Li",     1,         "Already elemental",
  "Mn",           "Mn",     1,         "Already elemental",
  "Mo",           "Mo",     1,         "Already elemental",
  "N",            "N",      1,         "Already elemental",
  "Nb",           "Nb",     1,         "Already elemental",
  "Ni",           "Ni",     1,         "Already elemental",
  "P",            "P",      1,         "Already elemental",
  "Pb",           "Pb",     1,         "Already elemental",
  "Re",           "Re",     1,         "Already elemental",
  "S",            "S",      1,         "Already elemental",
  "Sb",           "Sb",     1,         "Already elemental",
  "Se",           "Se",     1,         "Already elemental",
  "Si",           "Si",     1,         "Already elemental",
  "Sn",           "Sn",     1,         "Already elemental",
  "Sr",           "Sr",     1,         "Already elemental",
  "Ta",           "Ta",     1,         "Already elemental",
  "Te",           "Te",     1,         "Already elemental",
  "Th",           "Th",     1,         "Already elemental",
  "U",            "U",      1,         "Already elemental",
  "V",            "V",      1,         "Already elemental",
  "W",            "W",      1,         "Already elemental",
  "Zn",           "Zn",     1,         "Already elemental",
  "Zr",           "Zr",     0.5512,    "Zirconium concentrate to Zr"
)

mining_direct <- mining_clean %>%
  inner_join(mining_conversion_table, by = "Commodity") %>%
  mutate(
    FM1 = Mining_Production * factor
  ) %>%
  select(
    Element,
    Year,
    Mining_Production,
    Note,
    Finish,
    FM1,
    Commodity,
    conversion_note
  )

# 06.3 Iron: time-dependent ore grade --------------------------------------
#
# The iron series reports iron-ore mass. Ore grade is held at 60% through
# 1900, declines linearly to 44.81% in 2018, and is held constant thereafter.

mining_Fe <- mining_clean %>%
  filter(Commodity == "Fe") %>%
  mutate(
    Fe_grade = case_when(
      Year <= FE_GRADE_START_YEAR ~ FE_GRADE_1900,
      
      Year > FE_GRADE_START_YEAR &
        Year <= FE_GRADE_END_YEAR ~
        FE_GRADE_1900 -
        (Year - FE_GRADE_START_YEAR) *
        (FE_GRADE_1900 - FE_GRADE_2018) /
        (FE_GRADE_END_YEAR - FE_GRADE_START_YEAR),
      
      Year > FE_GRADE_END_YEAR ~ FE_GRADE_2018,
      TRUE ~ NA_real_
    ),
    FM1 = Mining_Production * Fe_grade,
    Element = "Fe",
    conversion_note = paste0(
      "Iron ore multiplied by time-dependent Fe grade: ",
      FE_GRADE_1900,
      " in ",
      FE_GRADE_START_YEAR,
      " to ",
      FE_GRADE_2018,
      " in ",
      FE_GRADE_END_YEAR
    )
  ) %>%
  select(
    Element,
    Year,
    Mining_Production,
    Note,
    Finish,
    FM1,
    Commodity,
    conversion_note,
    Fe_grade
  )

## 06.4 Combined commodities: Mg and Na -------------------------------------
#
# sum_or_na() preserves NA when every contributing commodity is missing.

combine_mining_commodities <- function(
    data,
    factor_table,
    output_element,
    output_note
) {
  data %>%
    inner_join(factor_table, by = "Commodity") %>%
    mutate(
      component_FM1 = Mining_Production * factor
    ) %>%
    group_by(Year) %>%
    summarise(
      FM1 = sum_or_na(component_FM1),
      Mining_Production = sum_or_na(Mining_Production),
      Note = paste_unique(Note),
      Finish = if_else(
        all(is.na(Finish)),
        NA_character_,
        "Yes"
      ),
      Commodity = paste_unique(Commodity),
      .groups = "drop"
    ) %>%
    mutate(
      Element = output_element,
      conversion_note = output_note
    ) %>%
    select(
      Element,
      Year,
      Mining_Production,
      Note,
      Finish,
      FM1,
      Commodity,
      conversion_note
    )
}

magnesium_factors <- tribble(
  ~Commodity,       ~factor,
  "Mg_compounds",   0.60303,
  "Mg_metal",       1
)

mining_Mg <- combine_mining_commodities(
  mining_clean,
  magnesium_factors,
  output_element = "Mg",
  output_note = "Mg compounds converted to Mg and added to Mg metal"
)

sodium_factors <- tribble(
  ~Commodity, ~factor,
  "Na2CO3",   0.43400,
  "NaCl",     0.39339
)

mining_Na <- combine_mining_commodities(
  mining_clean,
  sodium_factors,
  output_element = "Na",
  output_note = "Na2CO3 and NaCl converted to Na and summed"
)

## 06.5 Rare-earth-element allocation ---------------------------------------
#
# REO production is allocated among individual REEs using:
#   - country-specific compositions from 1990 onward;
#   - the mean composition of Australia and the United States before 1990.
#
# The Chinese composition is a weighted mixture of Bayan Obo, Sichuan, and
# other deposits, reproducing the original spreadsheet assumptions.

check_required_columns(
  data_REO_content,
  c("Country", REE_ELEMENTS),
  "REO_content"
)

# Helper: calculate the mean composition of selected rows.
mean_ree_rows <- function(data, rows, region_name, share) {
  data[rows, , drop = FALSE] %>%
    select(all_of(REE_ELEMENTS)) %>%
    mutate(
      across(
        everything(),
        ~ parse_ree_value(.x)
      )
    ) %>%
    summarise(
      across(
        everything(),
        ~ mean(.x, na.rm = TRUE)
      )
    ) %>%
    mutate(
      Region = region_name,
      Share = share
    )
}

china_bayan_obo <- mean_ree_rows(
  data_REO_content,
  CHINA_BAYAN_OBO_ROWS,
  "Bayan Obo",
  CHINA_DEPOSIT_SHARES[["Bayan_Obo"]]
)

china_sichuan <- mean_ree_rows(
  data_REO_content,
  CHINA_SICHUAN_ROWS,
  "Sichuan",
  CHINA_DEPOSIT_SHARES[["Sichuan"]]
)

china_other <- mean_ree_rows(
  data_REO_content,
  CHINA_OTHER_ROWS,
  "Others",
  CHINA_DEPOSIT_SHARES[["Others"]]
)

china_ree_composition <- bind_rows(
  china_bayan_obo,
  china_sichuan,
  china_other
) %>%
  summarise(
    across(
      all_of(REE_ELEMENTS),
      ~ weighted.mean(.x, w = Share, na.rm = TRUE)
    )
  ) %>%
  mutate(Country = "China")

other_country_ree_composition <- data_REO_content %>%
  filter(
    Country %in% c(
      "Australia",
      "United States",
      "Russia",
      "India"
    )
  ) %>%
  select(Country, all_of(REE_ELEMENTS)) %>%
  mutate(
    across(
      all_of(REE_ELEMENTS),
      parse_ree_value
    )
  )

china_ree_composition_clean <- china_ree_composition %>%
  mutate(
    across(
      all_of(REE_ELEMENTS),
      as.numeric
    )
  )

# Combine the two cleaned tables
reo_country_composition <- bind_rows(
  other_country_ree_composition,
  china_ree_composition_clean
)

# Complete annual REO production for the main countries.
#
# rule = 2 extends the first/last available value to the edges. This is an
# explicit extrapolation assumption.
reo_production_1990_2024 <- data_REO_production %>%
  mutate(
    Country = if_else(
      Country == "former USSR",
      "Russia",
      Country
    ),
    Year = as.numeric(Year),
    REO_production = as.numeric(REO_production)
  ) %>%
  group_by(Country) %>%
  complete(Year = 1990:2024) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(
    REO_production = zoo::na.approx(
      REO_production,
      x = Year,
      na.rm = FALSE,
      rule = 2
    )
  ) %>%
  ungroup()

# Convert each country's total REO production into individual elements.
reo_country_element_production <- reo_production_1990_2024 %>%
  left_join(reo_country_composition, by = "Country") %>%
  mutate(
    across(
      all_of(REE_ELEMENTS),
      ~ REO_production * (.x / 100)
    )
  ) %>%
  pivot_longer(
    cols = all_of(REE_ELEMENTS),
    names_to = "Element",
    values_to = "Element_production"
  )

# Derive the global REE composition in each year.
reo_global_composition_1990_2024 <- reo_country_element_production %>%
  group_by(Year, Element) %>%
  summarise(
    Element_production =
      sum_or_na(Element_production),
    .groups = "drop"
  ) %>%
  group_by(Year) %>%
  mutate(
    Total_REE_production =
      sum_or_na(Element_production),
    Element_content =
      Element_production / Total_REE_production
  ) %>%
  ungroup() %>%
  select(Year, Element, Element_content)

# Before 1990, use the average composition of Australia and the United States.
reo_composition_pre_1990 <- data_REO_content %>%
  filter(Country %in% c("Australia", "United States")) %>%
  select(Country, all_of(REE_ELEMENTS)) %>%
  mutate(
    across(
      all_of(REE_ELEMENTS),
      ~ parse_ree_value(.x)
    )
  ) %>%
  summarise(
    across(
      all_of(REE_ELEMENTS),
      ~ mean(.x, na.rm = TRUE) / 100
    )
  ) %>%
  crossing(Year = 1850:1989) %>%
  pivot_longer(
    cols = all_of(REE_ELEMENTS),
    names_to = "Element",
    values_to = "Element_content"
  )

reo_global_composition <- bind_rows(
  reo_composition_pre_1990,
  reo_global_composition_1990_2024
)

mining_REE <- mining_clean %>%
  filter(Commodity %in% REE_ELEMENTS) %>%
  transmute(
    Element = Commodity,
    Year,
    Mining_Production,
    Note,
    Finish,
    Commodity
  ) %>%
  left_join(
    reo_global_composition,
    by = c("Year", "Element")
  ) %>%
  mutate(
    FM1 = Mining_Production * Element_content,
    conversion_note =
      "Total REO production allocated using annual elemental composition"
  ) %>%
  select(
    Element,
    Year,
    Mining_Production,
    Note,
    Finish,
    FM1,
    Commodity,
    conversion_note,
    Element_content
  )

## 06.6 Platinum-group-metal allocation -------------------------------------
#
# When Note is missing, production is treated as element-specific.
# When Note is "PGM" or "other PGM", the reported group total is multiplied by
# the corresponding allocation factor.

pgm_allocation_factors <- tribble(
  ~Element, ~Note,        ~factor,
  "Pd",     "PGM",        0.45209,
  "Pt",     "PGM",        0.42477,
  "Ir",     "PGM",        0.015874,
  "Rh",     "PGM",        0.04918,
  "Ru",     "PGM",        0.06462,
  "Ir",     "other PGM",  0.12242,
  "Rh",     "other PGM",  0.39912,
  "Ru",     "other PGM",  0.47188
)

PGM_ELEMENTS <- c("Pd", "Pt", "Ir", "Rh", "Ru")

pgm_raw_diagnostic <- mining_clean %>%
  filter(Commodity %in% PGM_ELEMENTS) %>%
  count(
    Commodity,
    Note,
    sort = TRUE,
    name = "n_rows"
  )

print(pgm_raw_diagnostic, n = Inf)

mining_PGM <- mining_clean %>%
  filter(Commodity %in% PGM_ELEMENTS) %>%
  transmute(
    Element = Commodity,
    Year,
    Mining_Production,
    Note,
    Finish,
    Commodity
  ) %>%
  left_join(
    pgm_allocation_factors,
    by = c("Element", "Note")
  ) %>%
  mutate(
    factor = case_when(
      is.na(Note) ~ 1,
      !is.na(factor) ~ factor,
      TRUE ~ NA_real_
    ),
    FM1 = Mining_Production * factor,
    conversion_note = case_when(
      is.na(Note) ~ "Already element-specific",
      Note == "PGM" ~ "Allocated from total PGM production",
      Note == "other PGM" ~ "Allocated from other-PGM production",
      TRUE ~ "Unknown PGM allocation"
    )
  ) %>%
  select(
    Element,
    Year,
    Mining_Production,
    Note,
    Finish,
    FM1,
    Commodity,
    conversion_note,
    factor
  )

# PGM-factor diagnostics. Values may differ slightly from one because of
# rounding or because an element such as Os is not allocated.
pgm_factor_diagnostic <- pgm_allocation_factors %>%
  group_by(Note) %>%
  summarise(
    factor_sum = sum(factor),
    .groups = "drop"
  )

## 06.7 Final FM1 table ------------------------------------------------------
#
# The use of bind_rows() avoids the repeated rbind() blocks in the original
# code and makes omissions easier to detect.

df_mining <- bind_rows(
  mining_direct,
  mining_Fe,
  mining_Mg,
  mining_Na,
  mining_REE,
  mining_PGM
) %>%
  arrange(Element, Year)

mob_FM1 <- df_mining %>%
  select(Element, Year, FM1)

# Identify source commodities that were not converted.
converted_commodities <- unique(
  c(
    mining_conversion_table$Commodity,
    "Fe",
    magnesium_factors$Commodity,
    sodium_factors$Commodity,
    REE_ELEMENTS,
    PGM_ELEMENTS
  )
)

unconverted_mining_commodities <- mining_clean %>%
  filter(!Commodity %in% converted_commodities) %>%
  distinct(Commodity) %>%
  arrange(Commodity)

write_csv(
  df_mining,
  file.path(TABULAR_RESULTS_DIR, "mining_FM1_detailed.csv")
)

write_csv(
  mob_FM1,
  file.path(TABULAR_RESULTS_DIR, "mining_FM1.csv")
)

write_csv(
  unconverted_mining_commodities,
  file.path(TABULAR_RESULTS_DIR, "mining_unconverted_commodities.csv")
)

# 07. Fossil-fuel elemental mobilization ----------------------------------
#
# Coal and oil masses are treated as deterministic after interpolation.
# Elemental-composition uncertainty is propagated using Monte Carlo draws.

## 07.1 Coal ----------------------------------------------------------------

coal_mass <- data_coal %>%
  select(Year, Coal_Tg)

set.seed(SEED)

coal_composition_draws <- el_coal_final %>%
  select(Element, med_frac, sd_frac) %>%
  crossing(sim = seq_len(N_SIM)) %>%
  group_by(Element) %>%
  mutate(
    comp_sim = case_when(
      is.na(med_frac) ~ NA_real_,
      is.na(sd_frac) | sd_frac <= 0 ~ med_frac[1],
      TRUE ~ rnorm(
        n(),
        mean = med_frac[1],
        sd = sd_frac[1]
      )
    ),
    comp_sim = clip(comp_sim, lower = 0, upper = 1)
  ) %>%
  ungroup() %>%
  select(Element, sim, comp_sim)

coal_sims <- coal_mass %>%
  crossing(sim = seq_len(N_SIM)) %>%
  left_join(coal_composition_draws, by = "sim") %>%
  transmute(
    Year,
    Element,
    sim,
    mob_sim = Coal_Tg * comp_sim,
    flux = "FF_coal"
  )

coal_mobilization_summary <- coal_sims %>%
  group_by(Year, Element) %>%
  summarise_mc()

## 07.2 Oil -----------------------------------------------------------------

oil_mass <- data_oil %>%
  select(Year, Oil_Tg)

set.seed(SEED)

oil_composition_draws <- el_oil_final %>%
  crossing(sim = seq_len(N_SIM)) %>%
  group_by(Element) %>%
  mutate(
    comp_sim = case_when(
      is.na(med_frac) ~ NA_real_,
      
      distribution == "uniform" &
        !is.na(range_min_frac) &
        !is.na(range_max_frac) ~
        runif(
          n(),
          min = range_min_frac[1],
          max = range_max_frac[1]
        ),
      
      distribution %in% c("normal_cv", "normal_ucc") &
        !is.na(sd_frac) ~
        rnorm(
          n(),
          mean = med_frac[1],
          sd = sd_frac[1]
        ),
      
      !is.na(med_frac) ~ med_frac[1],
      TRUE ~ NA_real_
    ),
    comp_sim = clip(comp_sim, lower = 0, upper = 1)
  ) %>%
  ungroup() %>%
  select(
    Element,
    sim,
    comp_sim,
    source,
    distribution
  )

oil_sims <- oil_mass %>%
  crossing(sim = seq_len(N_SIM)) %>%
  left_join(oil_composition_draws, by = "sim") %>%
  transmute(
    Year,
    Element,
    sim,
    mob_sim = Oil_Tg * comp_sim,
    flux = "FF_oil"
  )

oil_mobilization_summary <- oil_sims %>%
  group_by(Year, Element) %>%
  summarise_mc()

write_csv(
  coal_mobilization_summary,
  file.path(TABULAR_RESULTS_DIR, "coal_mobilization_summary.csv")
)

write_csv(
  oil_mobilization_summary,
  file.path(TABULAR_RESULTS_DIR, "oil_mobilization_summary.csv")
)

# 08. Total anthropogenic elemental mobilization ---------------------------
#
# The total combines:
#   FM1      = mining or synthetic production;
#   FF_coal  = coal-related mobilization;
#   FF_oil   = oil-related mobilization;
#   FC       = construction-material mobilization.
#
# FM1 is currently deterministic and is repeated across all simulations.

set.seed(SEED)

## 08.1 Prepare each flux in a common structure -----------------------------

fm1_sims <- mob_FM1 %>%
  rename(mob_sim = FM1) %>%
  crossing(sim = seq_len(N_SIM)) %>%
  transmute(
    Element,
    Year,
    sim,
    mob_sim,
    flux = "FM1"
  )

coal_flux_sims <- coal_sims %>%
  select(Element, Year, sim, mob_sim) %>%
  mutate(flux = "FF_coal")

oil_flux_sims <- oil_sims %>%
  select(Element, Year, sim, mob_sim) %>%
  mutate(flux = "FF_oil")

construction_flux_sims <- mob_construction_sims %>%
  select(Element, Year, sim, mob_sim) %>%
  mutate(flux = "FC")

flux_sims <- bind_rows(
  fm1_sims,
  coal_flux_sims,
  oil_flux_sims,
  construction_flux_sims
)

## 08.2 Aggregate fluxes within each simulation -----------------------------
#
# sum_or_na() means:
#   - some known fluxes + some missing fluxes -> sum of known fluxes;
#   - all fluxes missing -> NA, not zero.
#
# This retains the original partial-sum logic while avoiding false zeros when
# no information is available at all.

mob_total_sims <- flux_sims %>%
  group_by(Element, Year, sim) %>%
  summarise(
    total_sim = sum_or_na(mob_sim),
    n_flux_observed = sum(!is.na(mob_sim)),
    .groups = "drop"
  )

## 08.3 Annual summaries ----------------------------------------------------

mob_total <- mob_total_sims %>%
  group_by(Element, Year) %>%
  summarise(
    total_mean = mean(total_sim, na.rm = TRUE),
    total_median = median(total_sim, na.rm = TRUE),
    total_sd = sd(total_sim, na.rm = TRUE),
    total_q25 = quantile(total_sim, 0.25, na.rm = TRUE),
    total_q75 = quantile(total_sim, 0.75, na.rm = TRUE),
    min_fluxes_observed = min(n_flux_observed, na.rm = TRUE),
    max_fluxes_observed = max(n_flux_observed, na.rm = TRUE),
    .groups = "drop"
  )

## 08.4 Log-ratio comparisons -----------------------------------------------
#
# log(M_t2 / M_t1):
#   > 0 means mobilization increased;
#   = 0 means no change;
#   < 0 means mobilization decreased.
#
# Ratios are calculated only when both endpoint values are positive.

comparison_years <- c(1850, 1950, 1960, 1980, 2000, 2020)

mob_ratio_sims <- mob_total_sims %>%
  filter(Year %in% comparison_years) %>%
  select(Element, Year, sim, total_sim) %>%
  pivot_wider(
    names_from = Year,
    values_from = total_sim
  ) %>%
  transmute(
    Element,
    sim,
    `1850–2020` = if_else(
      `1850` > 0 & `2020` > 0,
      log(`2020` / `1850`),
      NA_real_
    ),
    `1950–2020` = if_else(
      `1950` > 0 & `2020` > 0,
      log(`2020` / `1950`),
      NA_real_
    ),
    `1960–1980` = if_else(
      `1960` > 0 & `1980` > 0,
      log(`1980` / `1960`),
      NA_real_
    ),
    `1980–2000` = if_else(
      `1980` > 0 & `2000` > 0,
      log(`2000` / `1980`),
      NA_real_
    ),
    `2000–2020` = if_else(
      `2000` > 0 & `2020` > 0,
      log(`2020` / `2000`),
      NA_real_
    )
  ) %>%
  pivot_longer(
    cols = matches("–"),
    names_to = "comparison",
    values_to = "log_ratio"
  ) %>%
  filter(
    !is.na(log_ratio),
    !Element %in% c("Re", "He")
  )

PERIOD_LEVELS <- c(
  "1850–2020",
  "1950–2020",
  "1960–1980",
  "1980–2000",
  "2000–2020"
)

mob_ratio_sims <- mob_ratio_sims %>%
  mutate(
    comparison = factor(
      comparison,
      levels = PERIOD_LEVELS
    )
  )

mob_ratio_summary <- mob_ratio_sims %>%
  group_by(Element, comparison) %>%
  summarise(
    q25 = quantile(log_ratio, 0.25, na.rm = TRUE),
    median = median(log_ratio, na.rm = TRUE),
    q75 = quantile(log_ratio, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

## 08.5 Save results ---------------------------------------------------------

write_csv(
  mob_total,
  file.path(TABULAR_RESULTS_DIR, "total_mobilization_summary.csv")
)

write_csv(
  mob_total_sims,
  file.path(DRAW_RESULTS_DIR, "total_mobilization_draws.csv")
)

write_csv(
  mob_ratio_summary,
  file.path(TABULAR_RESULTS_DIR, "mobilization_log_ratio_summary.csv")
)


# 09. Figures and diagnostic plots -----------------------------------------
#
# This section creates the example construction plot and the multi-period
# mobilization comparison figure from the original workflow.

# 09.1 Example construction-mobilization trajectories ---------------------

construction_example_plot <- mob_construction_total %>%
  filter(Element %in% c("Fe", "Ca")) %>%
  ggplot(
    aes(
      x = Year,
      y = mob_mean,
      colour = Element,
      fill = Element
    )
  ) +
  geom_ribbon(
    aes(
      ymin = mob_q25,
      ymax = mob_q75
    ),
    colour = NA,
    alpha = 0.30
  ) +
  geom_line(linewidth = 1.2) +
  labs(
    x = NULL,
    y = "Element mobilization from construction (Tg)",
    colour = "Element",
    fill = "Element"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 0.8
    )
  )


## 09.2 Construction-only endpoint comparisons -----------------------------

construction_year_pairs <- tribble(
  ~start_year, ~end_year, ~period,
  1850,        2020,      "1850–2020",
  1950,        2020,      "1950–2020"
)

construction_log_ratio_sims <- pmap_dfr(
  construction_year_pairs,
  function(start_year, end_year, period) {
    mob_construction_sims %>%
      filter(Year %in% c(start_year, end_year)) %>%
      select(Year, Element, sim, mob_sim) %>%
      pivot_wider(
        names_from = Year,
        values_from = mob_sim
      ) %>%
      filter(
        .data[[as.character(start_year)]] > 0,
        .data[[as.character(end_year)]] > 0
      ) %>%
      mutate(
        log_ratio = log(
          .data[[as.character(end_year)]] /
            .data[[as.character(start_year)]]
        ),
        period = period
      )
  }
)

construction_log_ratio_summary <- construction_log_ratio_sims %>%
  group_by(period, Element) %>%
  summarise(
    q25 = quantile(log_ratio, 0.25, na.rm = TRUE),
    median = median(log_ratio, na.rm = TRUE),
    q75 = quantile(log_ratio, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

construction_element_order <- construction_log_ratio_summary %>%
  filter(period == "1850–2020") %>%
  arrange(median) %>%
  pull(Element)

construction_ratio_plot <- construction_log_ratio_summary %>%
  mutate(
    Element = factor(
      Element,
      levels = construction_element_order
    )
  ) %>%
  ggplot(
    aes(
      x = Element,
      y = median
    )
  ) +
  geom_errorbar(
    aes(
      ymin = q25,
      ymax = q75
    ),
    width = 0.2
  ) +
  geom_point(size = 2) +
  facet_wrap(
    ~ period,
    scales = "free_y"
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "log(mobilization in end year / mobilization in start year)"
  ) +
  theme_bw(base_size = 14)

## 09.3 Total mobilization comparisons -------------------------------------

PERIOD_COLOURS <- c(
  "1850–2020" = "#4C9F70",
  "1950–2020" = "#3B8EA5",
  "1960–1980" = "#E0C34A",
  "1980–2000" = "#F08A5D",
  "2000–2020" = "#7A4E5D"
)

mob_ratio_medians <- mob_ratio_sims %>%
  group_by(Element, comparison) %>%
  summarise(
    median_ratio = median(log_ratio, na.rm = TRUE),
    .groups = "drop"
  )

make_period_plot <- function(focal_period) {
  focal_draws <- mob_ratio_sims %>%
    filter(comparison == focal_period)
  
  # Each panel is ordered by its own median values.
  element_order <- mob_ratio_medians %>%
    filter(comparison == focal_period) %>%
    arrange(median_ratio) %>%
    pull(Element)
  
  focal_draws <- focal_draws %>%
    mutate(
      Element = factor(
        Element,
        levels = element_order
      )
    )
  
  # Long-period panels show no overlays. Short-period panels show the medians
  # of the other short periods as small vertical ticks.
  overlay_medians <- if (
    focal_period %in% c("1850–2020", "1950–2020")
  ) {
    mob_ratio_medians %>%
      slice(0)
  } else {
    mob_ratio_medians %>%
      filter(
        comparison != focal_period,
        !comparison %in% c(
          "1850–2020",
          "1950–2020"
        )
      )
  }
  
  overlay_medians <- overlay_medians %>%
    mutate(
      Element = factor(
        Element,
        levels = element_order
      ),
      Element_number = as.numeric(Element)
    ) %>%
    filter(!is.na(Element_number))
  
  ggplot(
    focal_draws,
    aes(
      x = Element,
      y = log_ratio
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.45,
      colour = "grey45"
    ) +
    geom_boxplot(
      fill = PERIOD_COLOURS[[focal_period]],
      width = 0.68,
      outlier.shape = NA,
      linewidth = 0.3,
      alpha = 0.6
    ) +
    geom_segment(
      data = overlay_medians,
      aes(
        x = Element_number,
        xend = Element_number,
        y = median_ratio - 0.01,
        yend = median_ratio + 0.01,
        colour = comparison
      ),
      inherit.aes = FALSE,
      linewidth = 1.4,
      lineend = "round"
    ) +
    coord_flip() +
    scale_colour_manual(
      values = PERIOD_COLOURS,
      breaks = PERIOD_LEVELS,
      name = "Median of other periods"
    ) +
    labs(
      x = NULL,
      y = expression(log(M[t2] / M[t1])),
      title = focal_period
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11,
        hjust = 0.5
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(
        colour = "grey88",
        linewidth = 0.35
      ),
      axis.text.y = element_text(
        size = 8.3,
        colour = "grey15"
      ),
      axis.text.x = element_text(
        size = 10,
        colour = "grey15"
      )
    )
}

period_plots <- map(
  PERIOD_LEVELS,
  make_period_plot
)

total_ratio_figure <- wrap_plots(
  period_plots,
  ncol = 5,
  guides = "collect"
) &
  theme(
    legend.position = "top",
    legend.justification = "center"
  )

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "Figure_2_Draft.pdf"
  ),
  plot = total_ratio_figure,
  width = 18,
  height = 12,
  device = cairo_pdf,
  dpi = 300
)














# 10. Results tables: Selected-year mobilization tables -----------------------------------------
#
# Purpose:
#   Create separate summary tables for:
#     1. Mining mobilization (FM1)
#     2. Fossil-fuel mobilization (coal + oil)
#     3. Construction mobilization
#
#   For each element and selected year, the tables report:
#     - Mean estimated mobilization
#     - Standard deviation
#
# Units:
#   Tg of element mobilized per year
# ==============================================================================


# 1. Define the years included in the tables -------------------------------

selected_years <- c(
  1850,
  1950,
  1960,
  1980,
  2000,
  2020
)


# 2. Helper function for safe summation ------------------------------------

# This function returns NA when all values are missing.
#
# This avoids interpreting an entirely missing mobilization value as zero.

sum_or_na <- function(x) {
  
  if (all(is.na(x))) {
    return(NA_real_)
  }
  
  sum(x, na.rm = TRUE)
}


# 3. Helper function for summarising mobilization draws --------------------

summarise_selected_years <- function(draws, selected_years) {
  
  # Keep a record of every element included in the input dataset.
  all_elements <- sort(unique(draws$Element))
  
  draws %>%
    
    # Retain only the years required for the summary tables.
    filter(Year %in% selected_years) %>%
    
    # Calculate the mean estimate and uncertainty for each element and year.
    group_by(Element, Year) %>%
    summarise(
      estimate_Tg = mean(mob_sim, na.rm = TRUE),
      sd_Tg       = sd(mob_sim, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    
    # Replace NaN values produced when all simulations are missing.
    mutate(
      estimate_Tg = if_else(is.nan(estimate_Tg), NA_real_, estimate_Tg),
      sd_Tg       = if_else(is.nan(sd_Tg), NA_real_, sd_Tg)
    ) %>%
    
    # Ensure that every element has one row for every selected year.
    complete(
      Element = all_elements,
      Year = selected_years
    ) %>%
    
    arrange(Element, Year)
}

# ==============================================================================
# Mining mobilization
# ==============================================================================

mining_selected_draws <- fm1_sims %>%
  
  # Combine multiple mining entries if an element has more than one source.
  group_by(Element, Year, sim) %>%
  summarise(
    mob_sim = sum_or_na(mob_sim),
    .groups = "drop"
  )

table_mining_long <- summarise_selected_years(
  draws = mining_selected_draws,
  selected_years = selected_years
) %>%
  mutate(
    mobilization = "Mining",
    unit = "Tg"
  ) %>%
  select(
    mobilization,
    Element,
    Year,
    estimate_Tg,
    sd_Tg,
    unit
  )

table_mining_long

# Fossil-fuel mobilization

fossil_fuel_selected_draws <- bind_rows(
  
  coal_sims %>%
    select(
      Element,
      Year,
      sim,
      mob_sim
    ),
  
  oil_sims %>%
    select(
      Element,
      Year,
      sim,
      mob_sim
    )
  
) %>%
  
  # Add coal and oil mobilization within each simulation.
  group_by(Element, Year, sim) %>%
  summarise(
    mob_sim = sum_or_na(mob_sim),
    .groups = "drop"
  )

table_fossil_fuels_long <- summarise_selected_years(
  draws = fossil_fuel_selected_draws,
  selected_years = selected_years
) %>%
  mutate(
    mobilization = "Fossil fuels",
    unit = "Tg"
  ) %>%
  select(
    mobilization,
    Element,
    Year,
    estimate_Tg,
    sd_Tg,
    unit
  )

table_fossil_fuels_long

construction_selected_draws <- construction_flux_sims %>%
  
  # This aggregation is included as a safeguard against duplicate rows.
  group_by(Element, Year, sim) %>%
  summarise(
    mob_sim = sum_or_na(mob_sim),
    .groups = "drop"
  )

table_construction_long <- summarise_selected_years(
  draws = construction_selected_draws,
  selected_years = selected_years
) %>%
  mutate(
    mobilization = "Construction",
    unit = "Tg"
  ) %>%
  select(
    mobilization,
    Element,
    Year,
    estimate_Tg,
    sd_Tg,
    unit
  )

table_construction_long

# Convert long summary tables to wide format
make_wide_mobilization_table <- function(summary_table) {
  
  summary_table %>%
    select(
      Element,
      Year,
      estimate_Tg,
      sd_Tg
    ) %>%
    pivot_wider(
      names_from = Year,
      values_from = c(
        estimate_Tg,
        sd_Tg
      ),
      names_glue = "{Year}_{.value}",
      names_vary = "slowest"
    ) %>%
    arrange(Element)
} 
  
table_mining_wide <- make_wide_mobilization_table(
    table_mining_long
  )
  
table_fossil_fuels_wide <- make_wide_mobilization_table(
    table_fossil_fuels_long
  )
  
table_construction_wide <- make_wide_mobilization_table(
    table_construction_long
  )
}
