# ==============================================================================
# Project: GCB — Anthropogenic Effects on the Elementome
#
# Script: 01_Elemental_Mobilization.R
#
# Purpose:
#   Reconstruct and quantify annual anthropogenic elemental mobilization
#   associated with mining, construction materials, coal consumption, and oil
#   consumption, while propagating uncertainty through Monte Carlo simulation.
#
# Analytical framework:
#   Total anthropogenic elemental mobilization is reconstructed as the sum of:
#     - mining-derived elemental production (FM1);
#     - coal-related mobilization (FF_coal);
#     - oil-related mobilization (FF_oil); and
#     - construction-related mobilization (FC).
#
#   Mining mobilization is treated as deterministic. Uncertainty in fossil-fuel
#   elemental composition and construction-material mass/composition is
#   propagated using Monte Carlo simulation.
#
# Units:
#   Annual elemental mobilization is expressed in Tg yr-1. Elemental
#   compositions are internally standardized as mass fractions (kg/kg).
#
# Temporal coverage:
#   Annual mobilization series span 1800-2024 where source data are available;
#   principal historical comparisons focus on 1850-2020.
#
# Author:
#   Fernando Coello Sanz
#
# Affiliation:
#   CREAF
#
# Contact:
#   f.coello@creaf.uab.cat
#
# Repository:
#   https://github.com/fercoello/GCB-Human_Elementome
#
# Created:
#   2026-01-30
#
# Last updated:
#   2026-10-07
#
# R version:
#   R 4.5.1
#
# Main dependencies:
#   tidyverse, readxl, here, zoo
#
# Inputs:
#   data/data_FM1_mining_production.xlsx
#   data/data_REO.xlsx
#   data/data_construction.csv
#   data/data_fossil_fuel.xlsx
#   data/elemental_concentration/data_elemental_concentration.xlsx
#
# Outputs:
#   results/tabular_results/element_mobilization_construction_summary.csv
#   results/tabular_results/mining_FM1.csv
#   results/tabular_results/coal_mobilization_summary.csv
#   results/tabular_results/oil_mobilization_summary.csv
#   results/tabular_results/total_mobilization_summary.csv
#   results/tabular_results/mobilization_log_ratio_summary.csv
#   results/tabular_results/selected_year_mobilization_table.csv
#   results/draws/element_mobilization_construction_draws.csv
#   results/draws/total_mobilization_draws.csv
#
# Reproducibility:
#   A fixed random-number seed is used for Monte Carlo simulations. Run the
#   script from the root directory of the RStudio project. File paths are
#   constructed relative to the project root using here::here().
#
# License:
#   Code: MIT License
#   Data and documentation: CC BY 4.0
#
# Citation:
#   See repository citation information.
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

# Simulations for coal composition
N_COAL_COMPOSITION_SIM <- 5000L

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
# IEA conversion factors from TWh to Tg. See reference supplementary material
COAL_MASS_FACTOR <- 0.1228
OIL_MASS_FACTOR  <- 0.086

# Construction-material assumptions ----------------------------------------

# Proportion of carbonates extraction (Eggleston, 2006)
LIMESTONE_SHARE <- 0.85
DOLOMITE_SHARE  <- 0.15

# Global Proportion of different clays (Ito and Wagai (2017))
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
  "data/data_FM1_mining_production.xlsx",
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
# Rare-earth oxide (REO) production and elemental-composition data were
# obtained from the United States Geological Survey (USGS).
#
# `REO_country_production` contains country-level REO production data used to
# estimate country- and year-specific REE composition. Data were obtained from
# Table T8 'RARE EARTHS: WORLD MINE PRODUCTION, BY COUNTRY OR LOCALITY' of the USGS Mineral Yearbooks XLSX releases:
# https://www.usgs.gov/centers/national-minerals-information-center/
# rare-earths-statistics-and-information
# Country-level production data are available from 1990 onward in several Excel files.
#
# `REO_content` contains the elemental composition of REO deposits/mines used
# to estimate country- and year-specific REE composition. Data were obtained
# from Table T2: 'RARE EARTH CONTENTS OF SELECTED SOURCE MINERALS' of the USGS Mineral Yearbooks XLSX releases available from the
# same source.
#
# Values were transcribed/compiled from the original USGS sources without
# modification. Complete references, temporal coverage, and the procedure used
# to combine production and composition data are described in the
# Supplementary Methods.
#
# Original third-party data are not redistributed in this repository.

data_REO_production <- read_excel(
  file.path(DATA_DIR, "data_REO.xlsx"),
  sheet = "REO_country_production"
)

data_REO_content <- read_excel(
  file.path(DATA_DIR, "data_REO.xlsx"),
  sheet = "REO_content"
)

## Construction-material mass series ----------------------------------------

# Annual global construction-material mass fluxes were reconstructed for
# aggregates, clay, and limestone using a dedicated preprocessing workflow.
#
# The reconstruction combines:
#   - historical concrete, asphalt, and brick production from Plank et al. (2022);
#   - USGS cement-production statistics;
#   - United Nations bitumen non-energy-use statistics;
#   - global population data from Our World in Data;
#   - virgin aggregate fractions from Krausmann et al. (2017);
#   - historical and World Bank railway-length statistics;
#   - country-specific railway-gauge information; and
#   - lime-production data from the USGS.
#
# The preprocessing script reconstructs annual material demand associated with
# concrete, asphalt, road and building bases, railway ballast, bricks, cement,
# and lime, and propagates the corresponding uncertainty.
#
# `data_construction.csv` is the resulting harmonized analysis-ready dataset.
# It contains annual best estimates and SDs for:
#   - aggregates;
#   - clay; and
#   - limestone.
#
# Masses are expressed in Tg yr-1.
#
# The full reconstruction procedure, assumptions, source references, and
# validation checks are documented in the Supplementary Methods and in
# `00_data_construction_creation.R`.

data_construction <- read_csv(
  file.path(DATA_DIR, "data_construction.csv"),
  show_col_types = FALSE
)

## Fossil-fuel consumption --------------------------------------------------

# Historical global fossil-fuel consumption was obtained from the Our World in
# Data compilation of coal, oil, and natural-gas consumption (Ritchie and
# Rosado, 2017), based on data from Smil (2017) and the Energy Institute
# (2025). The series extends back to 1800.
#
# Ritchie, H., Rosado, P., 2017. Fossil fuels. https://ourworldindata.org/fossil-fuels (accessed January 2026)
#
#
# Missing values in the historical consumption series are linearly interpolated
# between available observations.
#
# Global annual consumption is assumed to approximate global annual production.
# This assumption is appropriate for the present mass-flow reconstruction
# because changes in global stocks and trade imbalances are small relative to
# total annual flows and tend to balance at the global scale over longer time
# periods.
#
# Source values are expressed as energy consumption (TWh). Coal and oil are
# converted to mass using IEA equivalence factors:
#   - coal: 0.1228 Tg TWh-1
#   - oil:  0.0860 Tg TWh-1
# IEA, 2026. Unit Converter. https://www.iea.org/data-and-statistics/data-tools/unit-converter (accessed January 2026)
#
# These conversion factors are declared in the Setup section as
# COAL_MASS_FACTOR and OIL_MASS_FACTOR.
#
# Complete source references and methodological details are provided in the
# Supplementary Methods.

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

# Elemental-composition data used to estimate elemental mobilization associated
# with construction materials.
#
# The workbook contains both preprocessed material-specific compositions and
# source tables that are further processed within this script.
#
# Limestone, dolomite, and aggregate compositions are already prepared before
# being loaded here. In contrast, the UCC and clay datasets are further
# processed below to derive the final material-specific concentrations used in
# the mobilization calculations.
#
# Complete references, assumptions, selection criteria, and processing steps
# are described in the Supplementary Methods and in the README sheet of the
# workbook.

element_composition_file <- file.path(
  ELEMENT_CONCENTRATION_DIR,
  "data_elemental_concentration.xlsx"
)


### Upper Continental Crust ----------------------------------------------------

# Upper Continental Crust (UCC) elemental-composition data.
#
# These values are used later in this script as fallback concentrations when
# material-specific elemental data are unavailable or do not meet the minimum
# number of observations required for a material-specific estimate.

el_UCC_raw <- read_excel(
  element_composition_file,
  sheet = "UCC"
)


### Aggregates: river sand and gravel -----------------------------------------

# Analysis-ready elemental composition used to represent construction
# aggregates.
#
# Aggregate composition is approximated using the coarse fraction (>125 µm) of
# global river sediments from Müller et al. (2021). The aggregation and
# weighting procedure used to derive these values was performed prior to this
# script.

el_sand_gravel_raw <- read_excel(
  element_composition_file,
  sheet = "river_sand_gravel"
)


### Limestone -----------------------------------------------------------------

# Analysis-ready elemental composition of limestone derived from the
# whole-rock geochemical database compiled by Gard et al. (2019).
#
# The preprocessing, filtering, summary statistics, and oxide-to-element
# conversions used to generate these values are implemented in
# `00_data_elemental_composition.R`.

el_limestone_raw <- read_excel(
  element_composition_file,
  sheet = "limestone"
) %>%
  rename(
    Element = 1
  )


### Dolomite ------------------------------------------------------------------

# Analysis-ready elemental composition of dolomite derived from the
# whole-rock geochemical database compiled by Gard et al. (2019).
#
# As for limestone, preprocessing is performed in
# `00_data_elemental_composition.R`. Limestone and dolomite are subsequently
# combined within this script according to the assumed 85:15 contribution to
# global carbonate-rock extraction.

el_dolomite_raw <- read_excel(
  element_composition_file,
  sheet = "dolomite"
) %>%
  rename(
    Element = 1
  )


### Clay ----------------------------------------------------------------------

# Clay mineral-composition data used to derive the representative clay
# composition within this script.
#
# The input includes compositions for kaolinite, smectite, illite/mica, and
# post-Archean shale fallback values. These source compositions are combined
# below according to the assumed global relative abundance of the principal
# clay-mineral groups.
#
# Missing mineral-specific elemental concentrations and uncertainty estimates
# are also resolved later in this script according to the procedures described
# in the Supplementary Methods.

el_clays_raw <- read_excel(
  element_composition_file,
  sheet = "clay"
)

## Fossil-fuel elemental compositions ---------------------------------------

# Elemental-composition data used to estimate elemental mobilization associated
# with global coal and crude-oil consumption.
#
# The workbook contains separate coal datasets for:
#   1. major-element concentrations (C and S) by coal type;
#   2. trace-element concentrations; and
#   3. the historical contribution of different coal types to global
#      production.
#
# Crude-oil elemental concentrations are stored separately and are further
# processed below according to data availability and uncertainty.
#
# Representative concentrations, uncertainty distributions, fallback values,
# and annual weighted coal compositions are calculated later in this script.
# Complete references and methodological details are provided in the
# Supplementary Methods and in the workbook README.


### Coal: major elements -------------------------------------------------------

# Carbon and sulfur concentrations for the principal coal types considered in
# the analysis: anthracite, lignite, bituminous, sub-bituminous, and
# metallurgical coal.
#
# Concentration ranges and, where available, representative central estimates
# were compiled from Speight (2012), Breeze (2015), Kagiliery et al. (2019),
# and Nedelin et al. (2024).
#
# These values are combined later with annual coal-type production shares.
# For uncertainty propagation, triangular distributions are used when a
# representative central estimate is available together with minimum and
# maximum values; otherwise, uniform distributions bounded by the reported
# range are used.

el_coal_major_raw <- read_excel(
  element_composition_file,
  sheet = "coal_major"
)


### Coal: trace elements -------------------------------------------------------

# Trace-element concentration data for coal.
#
# Representative ranges were compiled primarily from Clarke and Sloss (1992),
# with additional values from Finkelman (1999) and other sources described in
# the Supplementary Methods.
#
# Where coal-specific concentrations are unavailable or unsuitable, fallback
# estimates are assigned later in the script, including UCC concentrations for
# a small number of elements.

el_coal_trace_raw <- read_excel(
  element_composition_file,
  sheet = "coal_trace"
)


### Historical coal-type production -------------------------------------------

# Annual global production by coal type, compiled from the U.S. Energy
# Information Administration (EIA).
#
# These data are converted below to annual proportions of anthracite, lignite,
# bituminous, sub-bituminous, and metallurgical coal. The resulting proportions
# are used to weight coal-type-specific C and S concentrations and reconstruct
# the annual elemental composition of the global coal mix.

coal_type_production_raw <- read_excel(
  element_composition_file,
  sheet = "coal_type_production"
)


### Crude oil -----------------------------------------------------------------

# Major- and trace-element concentration data for crude oil.
#
# Representative concentration ranges were compiled primarily from Coker
# (2018) and Hsu and Robinson (2019). For elements not covered by those
# compilations, concentration bounds were derived from published measurements
# of crude-oil samples, including Akinlua and Torto (2006), Dittert et al.
# (2009), Pereira et al. (2009), Sen and Peucker-Ehrenbrink (2012),
# Luz et al. (2013), Akinlua et al. (2015), Yang et al. (2017), and
# Seeger et al. (2019).
#
# Final representative values and uncertainty are derived later in this script
# according to a hierarchy based on data availability:
#
#   1. For elements with reported concentration ranges, a uniform distribution
#      bounded by the minimum and maximum is assumed.
#
#   2. For elements represented by at least two oil-specific measurements, the
#      minimum and maximum of those observations define a uniform distribution.
#
#   3. For elements represented by a single oil-specific value, uncertainty is
#      estimated using a reference coefficient of variation derived from trace
#      elements with bounded oil-specific ranges.
#
#   4. Remaining gaps are filled using UCC concentrations and associated
#      uncertainty as a last-resort proxy.

el_oil_raw <- read_excel(
  element_composition_file,
  sheet = "oil"
)


# 03. Construction-material elemental composition --------------------------
#
# This section derives the elemental compositions used to convert annual
# construction-material production into mobilized elemental masses.
#
# Final compositions are prepared for:
#   1. Upper Continental Crust (UCC), used as a fallback reference;
#   2. sand and gravel, representing construction aggregates;
#   3. limestone and dolomite, representing carbonate extraction; and
#   4. a globally weighted clay-mineral mixture.
#
# Material-specific concentrations are preferred whenever sufficiently
# supported by the available data. Where material-specific information is
# unavailable or insufficient, fallback compositions are applied according
# to the hierarchy described below.
#
# All final concentration estimates are converted to mass fractions
# (kg element per kg material) before being used in the mobilization analysis.

## 03.1 Upper Continental Crust ---------------------------------------------

# Upper Continental Crust (UCC) concentrations provide the final fallback
# composition for elements lacking sufficiently supported material-specific
# estimates.
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

# For elements lacking a reported UCC SD, uncertainty is estimated using the
# median coefficient of variation among elements for which both concentration
# and SD are available.ucc_cv_reference <- median(el_UCC$CV, na.rm = TRUE)

el_UCC <- el_UCC %>%
  mutate(
    Upper_crust_concentration_sd = coalesce(
      Upper_crust_concentration_sd,
      Upper_crust_concentration * ucc_cv_reference
    )
  )

# Convert UCC concentrations and uncertainty to mass fractions.
#
# The 25th and 75th percentiles are approximated assuming a normal
# distribution and are constrained to non-negative concentrations.

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

# Construction aggregates are represented using the elemental composition of
# coarse (>125 µm) global river sediments from Müller et al. (2021).
#
# The source dataset reports some major elements as oxides; these are converted
# to elemental concentrations using the stoichiometric factors defined above.
#
# For elements without a coarse-sediment estimate, UCC concentrations are used
# as fallback values.

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

# For values lacking an SD, uncertainty is estimated using the mean
# coefficient of variation across the available sediment/UCC estimates.

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

# Limestone concentrations are taken from the Gard et al. (2019) whole-rock
# database summaries generated in `00_data_elemental_composition.R`.
#
# Material-specific median and IQR estimates are used when at least
# MIN_CARBONATE_N observations are available for the element.
# Otherwise, the corresponding UCC concentration and uncertainty are used.

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

# Dolomite concentrations are also derived from the Gard et al. (2019)
# whole-rock database.
#
# The fallback hierarchy is:
#   1. dolomite-specific composition when n >= MIN_CARBONATE_N;
#   2. limestone composition when sufficient dolomite data are unavailable;
#   3. UCC composition when neither carbonate-specific estimate is available.
#
# This hierarchy ensures that carbonate-specific information is retained
# whenever possible before reverting to the general crustal composition.

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

## 03.6 Clay-mineral source compositions ------------------------------------

# Representative clay composition is reconstructed from three mineral groups:
# kaolinite, smectite, and illite/mica.
#
# Mineral-specific compositions are obtained from the clay source table loaded
# above. Major elements reported as oxides are converted to elemental
# concentrations using the same stoichiometric factors applied elsewhere.
#
# Missing mineral-specific concentrations are filled according to the
# following hierarchy:
#
#   mineral-specific value -> Post-Archean Australian Shale (PAAS) -> UCC
#
# PAAS therefore provides a geochemically closer fallback for clay minerals,
# while UCC is retained as the final fallback when neither mineral-specific
# nor shale values are available.

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

## 03.7 Weighted global clay composition ------------------------------------

# The final representative clay composition is obtained by combining
# kaolinite, smectite, and illite/mica according to their estimated global
# relative abundances defined in CLAY_WEIGHTS.
#
# Uncertainty is propagated by Monte Carlo simulation. For each element and
# clay-mineral group, concentrations are sampled independently from normal
# distributions defined by the estimated mean and SD. Negative concentrations
# are truncated to zero and mass fractions are constrained to a maximum of 1.
#
# Each simulated mineral concentration is multiplied by its global clay-mineral
# weight, and the weighted components are summed to obtain a distribution of
# global clay composition. The median, SD, and interquartile range of this
# distribution are retained for the mobilization analysis.

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

# 04. Fossil-fuel elemental composition ------------------------------------
#
# This section derives the elemental compositions used to estimate elemental
# mobilization associated with global coal and crude-oil consumption.
#
# The workflow:
#   1. reconstructs annual global production shares of the major coal types;
#   2. prepares coal C and S concentration ranges by coal type;
#   3. derives representative trace-element concentrations for coal;
#   4. propagates uncertainty in annual coal C and S composition using
#      Monte Carlo simulations; and
#   5. derives representative crude-oil elemental compositions and associated
#      uncertainty according to data availability.
#
# Upper Continental Crust (UCC) concentrations are used as a final fallback
# where fossil-fuel-specific elemental data are unavailable.
#
# Final elemental compositions are expressed as kg element per kg fuel.

## 04.1 Coal-type shares through time ---------------------------------------

# Annual global production by coal type is converted to relative proportions.
# These proportions are subsequently used to weight coal-type-specific C and S
# concentrations and reconstruct the annual composition of the global coal mix.

coal_type_shares <- coal_type_production_raw %>%
  rename(
    Anthracite = 3,
    Metallurgical = 4,
    Bituminous = 5,
    Subbituminous = 6,
    Lignite = 7
  ) %>%
  mutate(
    perc_Anthracite =
      Anthracite / Tota_Coal_Mst,
    
    perc_Metallurgical =
      Metallurgical / Tota_Coal_Mst,
    
    perc_Bituminous =
      Bituminous / Tota_Coal_Mst,
    
    perc_Subbituminous =
      Subbituminous / Tota_Coal_Mst,
    
    perc_Lignite =
      Lignite / Tota_Coal_Mst,
    
    # Diagnostic: the five proportions should sum approximately to one.
    share_sum =
      perc_Anthracite +
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


coal_share_diagnostic <- coal_type_shares %>%
  summarise(
    min_share_sum = min(
      share_sum,
      na.rm = TRUE
    ),
    
    max_share_sum = max(
      share_sum,
      na.rm = TRUE
    )
  )

## 04.2 Coal major-element composition --------------------------------------

# Carbon and sulfur concentrations are compiled separately for the major coal
# types considered in the analysis.
#
# When a representative value is explicitly reported in the source table it is
# retained as the central estimate. Otherwise, the midpoint of the reported
# minimum and maximum is used as a descriptive central value.

el_major_coal <- el_coal_major_raw %>%
  rename(
    Coal_type = 1
  ) %>%
  mutate(
    Coal_type =
      standardise_coal_type(Coal_type),
    
    has_reported_best =
      !is.na(as.numeric(Unique)),
    
    Mean = coalesce(
      as.numeric(Unique),
      (
        as.numeric(Min) +
          as.numeric(Max)
      ) / 2
    )
  ) %>%
  select(
    Coal_type,
    Element,
    Min,
    Mean,
    Max,
    has_reported_best
  )
## 04.3 Coal trace-element composition --------------------------------------
#
# Source priority reproduces the original workflow:
#   - Si: Finkelman;
#   - all other elements: IEA where available, then Finkelman;
#   - missing values: UCC fallback.
#
# IEA minima and maxima are interpreted as the 1st and 99th percentiles of an
# approximately normal distribution.

el_coal_trace <- el_coal_trace_raw %>%
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

# Annual global coal C and S compositions are reconstructed by combining
# coal-type-specific concentration distributions with annual coal-type
# production shares.
#
# When a literature-derived representative concentration is available together
# with minimum and maximum values, a triangular distribution is used with the
# representative value as its mode.
#
# When only minimum and maximum values are available, a uniform distribution
# bounded by those values is used.

### Draw random values from a triangular distribution -------------------------
rtriangular <- function(
    n,
    min,
    mode,
    max
) {
  
  if (
    is.na(min) |
    is.na(mode) |
    is.na(max)
  ) {
    return(
      rep(
        NA_real_,
        n
      )
    )
  }
  
  if (
    min > mode |
    mode > max
  ) {
    stop(
      "Triangular distribution requires min <= mode <= max."
    )
  }
  
  if (
    min == max
  ) {
    return(
      rep(
        min,
        n
      )
    )
  }
  
  u <- runif(
    n
  )
  
  threshold <-
    (mode - min) /
    (max - min)
  
  ifelse(
    u < threshold,
    
    min +
      sqrt(
        u *
          (max - min) *
          (mode - min)
      ),
    
    max -
      sqrt(
        (1 - u) *
          (max - min) *
          (max - mode)
      )
  )
}

simulate_coal_major_element <- function(
    element_symbol
) {
  
  parameters <- el_major_coal %>%
    filter(
      Element == element_symbol
    ) %>%
    transmute(
      Coal_type =
        standardise_coal_type(Coal_type),
      
      Min =
        as.numeric(Min),
      
      Mean =
        as.numeric(Mean),
      
      Max =
        as.numeric(Max),
      
      has_reported_best
    )
  
  
  set.seed(
    SEED
  )
  
  
  composition_draws <- parameters %>%
    mutate(
      draws = pmap(
        list(
          Min,
          Mean,
          Max,
          has_reported_best
        ),
        
        function(
    Min,
    Mean,
    Max,
    has_reported_best
        ) {
          
          if (has_reported_best) {
            
            rtriangular(
              N_COAL_COMPOSITION_SIM,
              min = Min,
              mode = Mean,
              max = Max
            )
            
          } else {
            
            runif(
              N_COAL_COMPOSITION_SIM,
              min = Min,
              max = Max
            )
          }
        }
      )
    ) %>%
    
    select(
      Coal_type,
      draws
    ) %>%
    
    unnest_longer(
      draws,
      indices_to = "sim"
    ) %>%
    
    rename(
      content = draws
    )
  
  
  shares_long <- coal_type_shares %>%
    select(
      -share_sum
    ) %>%
    
    pivot_longer(
      cols = starts_with(
        "perc_"
      ),
      names_to = "Coal_type",
      values_to = "share"
    ) %>%
    
    mutate(
      Coal_type =
        standardise_coal_type(
          Coal_type
        )
    )
  
  
  annual_draws <- shares_long %>%
    inner_join(
      composition_draws,
      by = "Coal_type"
    ) %>%
    
    mutate(
      weighted_content =
        share * content
    ) %>%
    
    group_by(
      Year,
      sim
    ) %>%
    
    summarise(
      annual_content =
        sum_or_na(
          weighted_content
        ),
      
      .groups = "drop"
    )
  
  
  annual_summary <- annual_draws %>%
    group_by(
      Year
    ) %>%
    
    summarise(
      mean =
        mean(
          annual_content,
          na.rm = TRUE
        ),
      
      median =
        median(
          annual_content,
          na.rm = TRUE
        ),
      
      sd =
        sd(
          annual_content,
          na.rm = TRUE
        ),
      
      q2.5 =
        quantile(
          annual_content,
          0.025,
          na.rm = TRUE
        ),
      
      q25 =
        quantile(
          annual_content,
          0.25,
          na.rm = TRUE
        ),
      
      q75 =
        quantile(
          annual_content,
          0.75,
          na.rm = TRUE
        ),
      
      q97.5 =
        quantile(
          annual_content,
          0.975,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    ) %>%
    
    mutate(
      Element =
        element_symbol
    )
  
  
  list(
    draws = annual_draws,
    summary = annual_summary
  )
}

coal_C_mc <- simulate_coal_major_element(
  "C"
)

coal_S_mc <- simulate_coal_major_element(
  "S"
)

coal_major_yearly_summary <- bind_rows(
  coal_C_mc$summary,
  coal_S_mc$summary
)


# The original workflow reduces the annual C and S composition reconstructions
# to one representative composition for each element by averaging the annual
# summaries across the available period.
#
# This constant representative composition is then applied to the complete
# historical coal-consumption series.
#
# Using the annual coal compositions directly would constitute a different
# methodology and would produce a different elemental-mobilization series.

coal_major_average <- coal_major_yearly_summary %>%
  group_by(
    Element
  ) %>%
  
  summarise(
    med_frac =
      mean(
        mean,
        na.rm = TRUE
      ) / 100,
    
    q25_frac =
      mean(
        q25,
        na.rm = TRUE
      ) / 100,
    
    q75_frac =
      mean(
        q75,
        na.rm = TRUE
      ) / 100,
    
    sd_frac =
      (
        q75_frac -
          q25_frac
      ) /
      (
        Z_75 -
          Z_25
      ),
    
    source =
      "Coal_MC_average",
    
    unit_best =
      "kg/kg",
    
    .groups =
      "drop"
  )

el_coal_final <- el_coal_with_ucc %>%
  filter(!Element %in% c("C", "S")) %>%
  bind_rows(coal_major_average) %>%
  arrange(Element)

## 04.5 Oil composition ------------------------------------------------------
#
# Crude-oil elemental composition is reconstructed using a hierarchical
# approach that depends on the amount of oil-specific information available
# for each element.
#
# Priority order:
#   1. If a representative minimum and maximum concentration are available,
#      assume a uniform distribution bounded by those values.
#
#   2. If no typical range is available but at least two oil-specific values
#      are reported, define a uniform distribution using the minimum and
#      maximum of all available values.
#
#   3. If only one oil-specific value is available, retain that concentration
#      as the central estimate and estimate its uncertainty later using a
#      reference coefficient of variation derived from trace elements.
#
#   4. If no oil-specific concentration is available, use the corresponding
#      Upper Continental Crust (UCC) concentration as a final fallback.
#
# The final composition table is expressed as kg element per kg oil.


## 04.5.1 Derive oil-specific central values and uncertainty -----------------

oil_composition_raw <- el_oil_raw %>%
  rowwise() %>%
  mutate(
    summary = list({
      
      has_typical <-
        !is.na(Crude_oil_min_typical) &&
        !is.na(Crude_oil_max_typical)
      
      has_yang <-
        !is.na(Crude_oil_Yang)
      
      has_sen <-
        !is.na(Crude_oil_Sen)
      
      has_samples <-
        !is.na(Crude_oil_min_samples) ||
        !is.na(Crude_oil_max_samples)
      
      
      all_values <- c(
        Crude_oil_Sen,
        Crude_oil_Yang,
        Crude_oil_min_samples,
        Crude_oil_max_samples
      )
      
      observed_values <-
        all_values[
          !is.na(all_values)
        ]
      
      
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
      
      

      # Case 1: representative typical range available

      
      if (has_typical) {
        
        # Minimum and maximum define a uniform distribution.
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
        
        

        # Case 2a: Yang provides the only oil-specific concentration

        
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
        
        

        # Case 2b: Sen provides the only oil-specific concentration

        
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
        
        

        # Case 3: at least two oil-specific concentrations available

        
      } else if (length(observed_values) >= 2) {
        
       
        range_min <-
          min(observed_values)
        
        range_max <-
          max(observed_values)
        
        stats <- uniform_stats(
          range_min,
          range_max
        )
        
        
        source_label <- case_when(
          has_sen &
            has_yang &
            has_samples ~
            "Oil_Sen+Yang+samples",
          
          has_sen &
            has_yang ~
            "Oil_Sen+Yang",
          
          has_sen &
            has_samples ~
            "Oil_Sen+samples",
          
          has_yang &
            has_samples ~
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
        
        

        # Case 4: one oil-specific concentration available

        
      } else if (length(observed_values) == 1) {
        
       
        tibble(
          median = observed_values,
          sd = NA_real_,
          q25 = NA_real_,
          q75 = NA_real_,
          
          source = case_when(
            has_sen ~
              "Oil_Sen",
            
            has_yang ~
              "Oil_Yang",
            
            TRUE ~
              "Oil_single"
          ),
          
          distribution = "normal_cv",
          range_min = NA_real_,
          range_max = NA_real_
        )
        
        

        # Case 5: no oil-specific information

        
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
  
  unnest(
    summary
  ) %>%
  
  ungroup() %>%
  
  transmute(
    Element,
    med_raw = median,
    sd_raw = sd,
    q25_raw = q25,
    q75_raw = q75,
    range_min_raw = range_min,
    range_max_raw = range_max,
    
    # Oil concentration data are expressed in mg/kg in the source table.
    unit_raw = "mg/kg",
    
    source,
    distribution
  )

## 04.5.2 Reference uncertainty for single-value estimates ------------------

# Elements represented by a single oil-specific concentration do not have an
# empirical uncertainty range. Their uncertainty is therefore estimated using
# a reference coefficient of variation (CV).
#
# The reference CV is calculated as the median relative uncertainty among
# trace elements represented by oil-specific uniform ranges.
#
# Major elements are excluded because their relative variability is not
# considered representative of trace-element uncertainty.

oil_major_elements <- c(
  "C",
  "H",
  "O",
  "N",
  "S",
  "Cl"
)

oil_trace_cv_reference <- oil_composition_raw %>%
  filter(
    !Element %in% oil_major_elements,
    distribution == "uniform",
    !is.na(med_raw),
    !is.na(sd_raw),
    med_raw > 0
  ) %>%
  
  summarise(
    CV = median(
      sd_raw / med_raw,
      na.rm = TRUE
    )
  ) %>%
  
  pull(
    CV
  )

## 04.5.3 Add uncertainty to single-value oil estimates ---------------------

# Apply the reference trace-element CV to elements represented by only one
# oil-specific concentration.
#
# The reported concentration is treated as the central value. SD is calculated
# from the reference CV, and the 25th and 75th percentiles are estimated
# assuming an approximately normal distribution.
#
# The lower quartile is constrained to zero to avoid non-physical negative
# concentrations.

oil_composition_completed <- oil_composition_raw %>%
  mutate(
    sd_raw = if_else(
      distribution == "normal_cv" &
        !is.na(med_raw),
      
      med_raw *
        oil_trace_cv_reference,
      
      sd_raw
    ),
    
    q25_raw = if_else(
      distribution == "normal_cv" &
        !is.na(med_raw),
      
      clip(
        med_raw +
          Z_25 *
          med_raw *
          oil_trace_cv_reference,
        lower = 0
      ),
      
      q25_raw
    ),
    
    q75_raw = if_else(
      distribution == "normal_cv" &
        !is.na(med_raw),
      
      clip(
        med_raw +
          Z_75 *
          med_raw *
          oil_trace_cv_reference,
        lower = 0
      ),
      
      q75_raw
    ),
    
    # Update provenance labels to indicate that uncertainty for these values
    # was derived using the reference trace-element CV.
    source = case_when(
      source == "Oil_Yang" ~
        "Oil_Yang_CVtrace",
      
      source == "Oil_Sen" ~
        "Oil_Sen_CVtrace",
      
      source == "Oil_single" ~
        "Oil_single_CVtrace",
      
      TRUE ~
        source
    )
  )

## 04.5.4 Convert to mass fractions and apply UCC fallback ------------------

# Convert all oil-specific concentrations to kg element per kg oil.
#
# Elements without an oil-specific estimate are then completed using the UCC
# composition and associated uncertainty.

el_oil_final <- oil_composition_completed %>%
  mutate(
    from_oil =
      !is.na(med_raw),
    
    med_frac_oil =
      to_mass_fraction(
        med_raw,
        unit_raw
      ),
    
    sd_frac_oil =
      to_mass_fraction(
        sd_raw,
        unit_raw
      ),
    
    q25_frac_oil = clip(
      to_mass_fraction(
        q25_raw,
        unit_raw
      ),
      lower = 0
    ),
    
    q75_frac_oil = clip(
      to_mass_fraction(
        q75_raw,
        unit_raw
      ),
      lower = 0
    ),
    
    range_min_frac_oil = if_else(
      !is.na(range_min_raw),
      
      clip(
        to_mass_fraction(
          range_min_raw,
          unit_raw
        ),
        lower = 0
      ),
      
      NA_real_
    ),
    
    range_max_frac_oil = if_else(
      !is.na(range_max_raw),
      
      clip(
        to_mass_fraction(
          range_max_raw,
          unit_raw
        ),
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
    
    med_frac =
      coalesce(
        med_frac_oil,
        med_frac_ucc
      ),
    
    sd_frac =
      coalesce(
        sd_frac_oil,
        sd_frac_ucc
      ),
    
    q25_frac =
      coalesce(
        q25_frac_oil,
        q25_frac_ucc
      ),
    
    q75_frac =
      coalesce(
        q75_frac_oil,
        q75_frac_ucc
      ),
    
    # Range bounds are retained only for oil-specific uniform distributions.
    range_min_frac =
      range_min_frac_oil,
    
    range_max_frac =
      range_max_frac_oil,
    
    unit_best =
      "kg/kg",
    
    source = case_when(
      coalesce(
        from_oil,
        FALSE
      ) ~
        source,
      
      !is.na(source_ucc) ~
        source_ucc,
      
      TRUE ~
        source
    ),
    
    distribution = case_when(
      coalesce(
        from_oil,
        FALSE
      ) ~
        distribution,
      
      !is.na(med_frac_ucc) ~
        "normal_ucc",
      
      TRUE ~
        distribution
    )
  ) %>%
  
  arrange(
    Element
  )

# 05. Construction-related elemental mobilization --------------------------
#
# Annual elemental mobilization associated with construction is calculated as:
#
#   elemental mobilization =
#     construction-material mass × elemental mass fraction
#
# Uncertainty in both annual material mass and material elemental composition
# is propagated using Monte Carlo simulations.
#
# Construction-material masses are expressed in Tg yr-1, whereas elemental
# compositions are expressed as mass fractions (kg element per kg material).
# Their product therefore gives elemental mobilization directly in Tg yr-1.
#
# Three construction-material streams are considered:
#   1. sand and gravel;
#   2. carbonate rocks (limestone and dolomite); and
#   3. clay.
#
# The resulting simulated elemental fluxes are subsequently summed to obtain
# total construction-related mobilization.

SEED <- 27
set.seed(SEED)

## 05.1 Sand and gravel -----------------------------------------------------

# Annual aggregate extraction and its uncertainty are obtained from the
# construction-material reconstruction.
#
# Material mass is represented by a normal distribution defined by the annual
# best estimate and SD. Simulated negative masses are constrained to zero.

aggregate_mass <- data_construction %>%
  transmute(
    Year,
    mass_mean = aggregates_best,
    mass_sd = aggregates_sd
  )


aggregate_mass_draws <- aggregate_mass %>%
  crossing(
    sim = seq_len(N_SIM)
  ) %>%
  mutate(
    mass_sim = rnorm(
      n(),
      mean = mass_mean,
      sd = mass_sd
    ),
    
    mass_sim = clip(
      mass_sim,
      lower = 0
    )
  ) %>%
  select(
    Year,
    sim,
    mass_sim
  )


# Concentrations are constrained between 0 and 1 because they are expressed
# as mass fractions.

aggregate_composition_draws <- el_sand_gravel %>%
  select(
    Element,
    med_frac,
    sd_frac
  ) %>%
  crossing(
    sim = seq_len(N_SIM)
  ) %>%
  group_by(
    Element
  ) %>%
  mutate(
    sd_use = if_else(
      is.na(sd_frac) |
        sd_frac < 0,
      0,
      sd_frac
    ),
    
    comp_sim = rnorm(
      n(),
      mean = med_frac[1],
      sd = sd_use[1]
    ),
    
    comp_sim = clip(
      comp_sim,
      lower = 0,
      upper = 1
    )
  ) %>%
  ungroup() %>%
  select(
    Element,
    sim,
    comp_sim
  )


mob_aggregate_sims <- aggregate_mass_draws %>%
  inner_join(
    aggregate_composition_draws,
    by = "sim"
  ) %>%
  transmute(
    Year,
    Element,
    sim,
    mob_sim = mass_sim * comp_sim,
    material = "sand and gravel"
  )


elem_mob_aggregate <- mob_aggregate_sims %>%
  group_by(
    Year,
    Element
  ) %>%
  summarise_mc() %>%
  mutate(
    unit = "Tg yr-1",
    material = "sand and gravel"
  )

## 05.2 Limestone and dolomite mixture -------------------------------------
#
# Construction carbonates are represented as a fixed mixture of limestone
# and dolomite according to the global extraction shares defined in Setup:
#
#   limestone = 85%
#   dolomite  = 15%
#
# For each element, the median and interquartile-range estimates of limestone
# and dolomite composition are combined using these fixed material shares.
#
# The resulting carbonate composition is then combined with annual carbonate
# extraction and its uncertainty to estimate elemental mobilization.


## Combine limestone and dolomite elemental compositions ---------------------

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
    # Limestone and dolomite values should be expressed in compatible units
    # before they are combined.
    units_match =
      unit_limestone == unit_dolomite |
      (
        is.na(unit_limestone) &
          is.na(unit_dolomite)
      ),
    
    unit =
      coalesce(
        unit_limestone,
        unit_dolomite
      ),
    
    # Weighted carbonate composition.
    med =
      LIMESTONE_SHARE * med_limestone +
      DOLOMITE_SHARE * med_dolomite,
    
    q25 =
      LIMESTONE_SHARE * q25_limestone +
      DOLOMITE_SHARE * q25_dolomite,
    
    q75 =
      LIMESTONE_SHARE * q75_limestone +
      DOLOMITE_SHARE * q75_dolomite,
    
    # Convert concentrations to mass fractions.
    med_frac =
      to_mass_fraction(
        med,
        unit
      ),
    
    q25_frac =
      to_mass_fraction(
        q25,
        unit
      ),
    
    q75_frac =
      to_mass_fraction(
        q75,
        unit
      )
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

### Annual carbonate-rock extraction -----------------------------------------

# Annual carbonate extraction is represented by the reconstructed limestone
# series, which includes both limestone and dolomite within the assumed
# carbonate-rock mixture.
#
# Material mass is expressed in Tg yr-1 and sampled from a normal distribution
# defined by the annual best estimate and SD. Negative simulated masses are
# constrained to zero.

carbonate_mass <- data_construction %>%
  transmute(
    Year,
    mass_mean = limestone_best,
    mass_sd = limestone_sd
  )


carbonate_mass_draws <- carbonate_mass %>%
  crossing(
    sim = seq_len(N_SIM)
  ) %>%
  
  mutate(
    mass_sim = rnorm(
      n(),
      mean = mass_mean,
      sd = mass_sd
    ),
    
    mass_sim = clip(
      mass_sim,
      lower = 0
    )
  ) %>%
  
  select(
    Year,
    sim,
    mass_sim
  )

### Carbonate elemental-composition uncertainty ------------------------------

# Carbonate elemental concentrations are represented by lognormal
# distributions parameterised using the 25th and 75th percentiles of the
# weighted limestone-dolomite composition.
#
# A lognormal distribution is used to maintain positive concentrations and to
# allow for right-skewed compositional uncertainty.

carbonate_composition_draws <- el_carbonate_mix %>%
  select(
    Element,
    q25_frac,
    q75_frac
  ) %>%
  
  group_by(
    Element
  ) %>%
  
  group_modify(
    ~ tibble(
      sim = seq_len(N_SIM),
      
      comp_sim = draw_lognormal_from_iqr(
        N_SIM,
        q25 = .x$q25_frac[1],
        q75 = .x$q75_frac[1]
      ) %>%
        clip(
          lower = 0,
          upper = 1
        )
    )
  ) %>%
  
  ungroup()

### Carbonate-related elemental mobilization ---------------------------------

# Combine simulated annual carbonate mass with simulated elemental
# composition.
#
# Tg material × kg/kg = Tg element.

mob_carbonate_sims <- carbonate_mass_draws %>%
  inner_join(
    carbonate_composition_draws,
    by = "sim"
  ) %>%
  transmute(
    Year,
    Element,
    sim,
    mob_sim = mass_sim * comp_sim,
    material = "limestone and dolomite"
  )

elem_mob_carbonate <- mob_carbonate_sims %>%
  group_by(
    Year,
    Element
  ) %>%
  summarise_mc() %>%
  mutate(
    unit = "Tg yr-1",
    material = "limestone and dolomite"
  )
## 05.3 Clay ----------------------------------------------------------------
#
# Annual clay demand and its uncertainty are obtained from the construction-
# material reconstruction.
#
# Material mass is represented by a normal distribution defined by the annual
# best estimate and SD. Negative simulated masses are constrained to zero.
#
# The representative clay elemental composition derived in Section 03 is
# sampled using normal distributions defined by the estimated central
# concentration and SD.
#
# Material mass is expressed in Tg yr-1 and elemental composition as kg/kg,
# so their product gives elemental mobilization directly in Tg yr-1.


### Annual clay mass -----------------------------------------------------------

clay_mass <- data_construction %>%
  transmute(
    Year,
    mass_mean = clay_best,
    mass_sd = clay_sd
  )


clay_mass_draws <- clay_mass %>%
  crossing(
    sim = seq_len(N_SIM)
  ) %>%
  mutate(
    mass_sim = rnorm(
      n(),
      mean = mass_mean,
      sd = mass_sd
    ),
    
    mass_sim = clip(
      mass_sim,
      lower = 0
    )
  ) %>%
  select(
    Year,
    sim,
    mass_sim
  )


### Clay elemental-composition uncertainty ------------------------------------

clay_composition_draws <- el_clay %>%
  select(
    Element,
    med_frac,
    sd_frac
  ) %>%
  crossing(
    sim = seq_len(N_SIM)
  ) %>%
  group_by(
    Element
  ) %>%
  mutate(
    sd_use = if_else(
      is.na(sd_frac) |
        sd_frac < 0,
      0,
      sd_frac
    ),
    
    comp_sim = rnorm(
      n(),
      mean = med_frac[1],
      sd = sd_use[1]
    ),
    
    comp_sim = clip(
      comp_sim,
      lower = 0,
      upper = 1
    )
  ) %>%
  ungroup() %>%
  select(
    Element,
    sim,
    comp_sim
  )


### Clay-related elemental mobilization ---------------------------------------

mob_clay_sims <- clay_mass_draws %>%
  inner_join(
    clay_composition_draws,
    by = "sim"
  ) %>%
  transmute(
    Year,
    Element,
    sim,
    
    # Tg material × kg/kg = Tg element.
    mob_sim = mass_sim * comp_sim,
    
    material = "clay"
  )


elem_mob_clay <- mob_clay_sims %>%
  group_by(
    Year,
    Element
  ) %>%
  summarise_mc() %>%
  mutate(
    unit = "Tg yr-1",
    material = "clay"
  )

## 05.4 Total construction mobilization ------------------------------------
#
# Elemental mobilization from sand and gravel, carbonate rocks, and clay is
# combined to obtain total construction-related mobilization.
#
# Monte Carlo draws are summed within each year, element, and simulation
# iteration. Using the same simulation identifier preserves the propagated
# uncertainty from the three construction-material streams.
#
# All component mobilization values are expressed in Tg yr-1, so their sum is
# also expressed in Tg yr-1.


## Combine simulated mobilization from all construction materials ------------

mob_construction_component_sims <- bind_rows(
  mob_aggregate_sims,
  mob_carbonate_sims,
  mob_clay_sims
)


## Sum material-specific mobilization within each Monte Carlo iteration -------

mob_construction_sims <- mob_construction_component_sims %>%
  group_by(
    Year,
    Element,
    sim
  ) %>%
  summarise(
    mob_sim = sum_or_na(
      mob_sim
    ),
    .groups = "drop"
  )


## Summarize the Monte Carlo distribution for each year and element -----------

mob_construction_total <- mob_construction_sims %>%
  group_by(
    Year,
    Element
  ) %>%
  summarise_mc() %>%
  mutate(
    unit = "Tg yr-1"
  )


write_csv(
  mob_construction_total,
  file.path(
    TABULAR_RESULTS_DIR,
    "element_mobilization_construction_summary.csv"
  )
)

write_csv(
  mob_construction_sims,
  file.path(
    DRAW_RESULTS_DIR,
    "element_mobilization_construction_draws.csv"
  )
)

# 06. Mining-derived elemental mobilization (FM1) --------------------------
#
# Mining production is converted to elemental mass (FM1) according to the
# form in which each commodity is reported in the source database.
#
# The mining dataset contains:
#   1. commodities already reported as elemental production;
#   2. compounds requiring stoichiometric conversion to elemental mass;
#   3. multiple commodities contributing to the same element;
#   4. iron ore requiring a time-dependent ore-grade correction;
#   5. rare-earth oxide (REO) production requiring allocation among REEs; and
#   6. platinum-group-metal (PGM) totals requiring allocation among elements.
#
# The result of this section is one harmonized annual elemental-production
# series for each supported element.


## 06.1 General cleaning ----------------------------------------------------


mining_clean <- data_mining %>%
  mutate(
    across(
      where(is.character),
      ~ na_if(
        str_trim(.x),
        "NA"
      )
    ),
    
    Year =
      as.numeric(Year),
    
    Mining_Production =
      to_numeric(Mining_Production),
    
    Note =
      na_if(
        str_trim(Note),
        ""
      )
  ) %>%
  rename(
    Commodity = Element
  )


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
# Commodities that are already reported as elemental production receive a
# conversion factor of 1.
#
# For compounds or concentrates, the conversion factor represents:
#
#   kg target element / kg reported commodity
#
# Multiplying reported production by this factor therefore converts commodity
# mass to elemental mass while retaining the original production mass unit.

mining_conversion_table <- tribble(
  ~Commodity, ~Element, ~factor, ~conversion_note,
  "Ag",       "Ag",     1,       "Already elemental",
  "Al",       "Al",     1,       "Already elemental",
  "As",       "As",     1,       "Already elemental",
  "Au",       "Au",     1,       "Already elemental",
  "B2O3",     "B",      0.31058, "B2O3 to B",
  "BaSO4",    "Ba",     0.58835, "BaSO4 to Ba",
  "Be",       "Be",     1,       "Already elemental",
  "Bi",       "Bi",     1,       "Already elemental",
  "Br",       "Br",     1,       "Already elemental",
  "Cd",       "Cd",     1,       "Already elemental",
  "Co",       "Co",     1,       "Already elemental",
  "Cr",       "Cr",     1,       "Already elemental",
  "Cs",       "Cs",     1,       "Already elemental",
  "Cu",       "Cu",     1,       "Already elemental",
  "F",        "F",      0.48667, "CaF2 to F; confirm source commodity",
  "Ga",       "Ga",     1,       "Already elemental",
  "Ge",       "Ge",     1,       "Already elemental",
  "He",       "He",     1,       "Already elemental",
  "Hf",       "Hf",     0.01531, "Zirconium concentrate to Hf",
  "In",       "In",     1,       "Already elemental",
  "K",        "K",      0.83,    "K2O to K",
  "Li",       "Li",     1,       "Already elemental",
  "Mn",       "Mn",     1,       "Already elemental",
  "Mo",       "Mo",     1,       "Already elemental",
  "N",        "N",      1,       "Already elemental",
  "Nb",       "Nb",     1,       "Already elemental",
  "Ni",       "Ni",     1,       "Already elemental",
  "P",        "P",      1,       "Already elemental",
  "Pb",       "Pb",     1,       "Already elemental",
  "Re",       "Re",     1,       "Already elemental",
  "S",        "S",      1,       "Already elemental",
  "Sb",       "Sb",     1,       "Already elemental",
  "Se",       "Se",     1,       "Already elemental",
  "Si",       "Si",     1,       "Already elemental",
  "Sn",       "Sn",     1,       "Already elemental",
  "Sr",       "Sr",     1,       "Already elemental",
  "Ta",       "Ta",     1,       "Already elemental",
  "Te",       "Te",     1,       "Already elemental",
  "Th",       "Th",     1,       "Already elemental",
  "U",        "U",      1,       "Already elemental",
  "V",        "V",      1,       "Already elemental",
  "W",        "W",      1,       "Already elemental",
  "Zn",       "Zn",     1,       "Already elemental",
  "Zr",       "Zr",     0.5512,  "Zirconium concentrate to Zr"
)


mining_direct <- mining_clean %>%
  inner_join(
    mining_conversion_table,
    by = "Commodity"
  ) %>%
  mutate(
    FM1 =
      Mining_Production *
      factor
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
## 06.3 Iron: time-dependent ore grade -------------------------------------
#
# The historical Fe mining series represents iron-ore mass rather than pure
# elemental Fe.
#
# Ore grade is assumed to:
#   - remain constant at FE_GRADE_1900 through FE_GRADE_START_YEAR;
#   - decline linearly to FE_GRADE_2018 by FE_GRADE_END_YEAR; and
#   - remain constant at the 2018 value thereafter.
#
# Elemental Fe production is calculated as:
#
#   iron-ore production × Fe grade.

mining_Fe <- mining_clean %>%
  filter(
    Commodity == "Fe"
  ) %>%
  mutate(
    Fe_grade = case_when(
      Year <= FE_GRADE_START_YEAR ~
        FE_GRADE_1900,
      
      Year > FE_GRADE_START_YEAR &
        Year <= FE_GRADE_END_YEAR ~
        FE_GRADE_1900 -
        (
          Year -
            FE_GRADE_START_YEAR
        ) *
        (
          FE_GRADE_1900 -
            FE_GRADE_2018
        ) /
        (
          FE_GRADE_END_YEAR -
            FE_GRADE_START_YEAR
        ),
      
      Year > FE_GRADE_END_YEAR ~
        FE_GRADE_2018,
      
      TRUE ~
        NA_real_
    ),
    
    FM1 =
      Mining_Production *
      Fe_grade,
    
    Element =
      "Fe",
    
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

## 06.4 Combined commodities: Mg and Na ------------------------------------
#
# Some elements are represented by more than one mining commodity.
#
# Each commodity is first converted to elemental mass and the resulting
# elemental contributions are then summed within each year.
#



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
# Total rare-earth oxide (REO) production is allocated among individual REEs
# using representative elemental compositions for the major producing
# countries considered in the analysis.
#
# From 1990 onward:
#   - country-specific REE compositions are used for Australia, China,
#     the United States, Russia, and India;
#   - China's composition is represented by a weighted mixture of Bayan Obo,
#     Sichuan, and other deposits.
#
# Before 1990, the global REE composition is approximated using the mean
# composition of Australia and the United States.
#
# These countries and deposits are retained because they represent the major
# REO producers considered in the original reconstruction.


check_required_columns(
  data_REO_content,
  c(
    "Country",
    REE_ELEMENTS
  ),
  "REO_content"
)


### Helper function ------------------------------------------------------------


mean_ree_rows <- function(
    data,
    rows,
    region_name,
    share
) {
  
  data[
    rows,
    ,
    drop = FALSE
  ] %>%
    
    select(
      all_of(
        REE_ELEMENTS
      )
    ) %>%
    
    mutate(
      across(
        everything(),
        ~ parse_ree_value(.x)
      )
    ) %>%
    
    summarise(
      across(
        everything(),
        ~ mean(
          .x,
          na.rm = TRUE
        )
      )
    ) %>%
    
    mutate(
      Region = region_name,
      Share = share
    )
}

### Chinese REE composition ----------------------------------------------------

# China is represented as a weighted mixture of three deposit groups. The
# deposit weights are defined in Setup and reproduce the assumptions used in
# the original reconstruction.

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
      all_of(
        REE_ELEMENTS
      ),
      ~ weighted.mean(
        .x,
        w = Share,
        na.rm = TRUE
      )
    )
  ) %>%
  
  mutate(
    Country = "China"
  )

### Other major producing countries -------------------------------------------

# Country-specific compositions for the other major REO producers
# included in the reconstruction.

other_country_ree_composition <- data_REO_content %>%
  filter(
    Country %in% c(
      "Australia",
      "United States",
      "Russia",
      "India"
    )
  ) %>%
  
  select(
    Country,
    all_of(
      REE_ELEMENTS
    )
  ) %>%
  
  mutate(
    across(
      all_of(
        REE_ELEMENTS
      ),
      parse_ree_value
    )
  )


china_ree_composition_clean <- china_ree_composition %>%
  mutate(
    across(
      all_of(
        REE_ELEMENTS
      ),
      as.numeric
    )
  )


reo_country_composition <- bind_rows(
  other_country_ree_composition,
  china_ree_composition_clean
)

### Annual REO production, 1990-2024 ------------------------------------------

# Complete annual REO production series for the selected major producing
# countries.
#
# Missing years between observations are linearly interpolated.
#
# `rule = 2` extends the first and last observed values to the temporal
# boundaries of the 1990-2024 reconstruction. This is therefore an explicit
# extrapolation assumption.

reo_production_1990_2024 <- data_REO_production %>%
  mutate(
    Country = if_else(
      Country == "former USSR",
      "Russia",
      Country
    ),
    
    Year =
      as.numeric(
        Year
      ),
    
    REO_production =
      as.numeric(
        REO_production
      )
  ) %>%
  
  group_by(
    Country
  ) %>%
  
  complete(
    Year = 1990:2024
  ) %>%
  
  arrange(
    Year,
    .by_group = TRUE
  ) %>%
  
  mutate(
    REO_production = zoo::na.approx(
      REO_production,
      x = Year,
      na.rm = FALSE,
      rule = 2
    )
  ) %>%
  
  ungroup()

### Country-level elemental REE production ------------------------------------

# Country-level REO production is distributed among individual REEs according
# to the corresponding country-specific elemental composition.

reo_country_element_production <- reo_production_1990_2024 %>%
  left_join(
    reo_country_composition,
    by = "Country"
  ) %>%
  
  mutate(
    across(
      all_of(
        REE_ELEMENTS
      ),
      ~ REO_production *
        (
          .x / 100
        )
    )
  ) %>%
  
  pivot_longer(
    cols =
      all_of(
        REE_ELEMENTS
      ),
    
    names_to =
      "Element",
    
    values_to =
      "Element_production"
  )

### Annual global REE composition ---------------------------------------------

# Element-specific production is summed across the selected major producing
# countries and converted to an annual global relative composition.

reo_global_composition_1990_2024 <- reo_country_element_production %>%
  group_by(
    Year,
    Element
  ) %>%
  
  summarise(
    Element_production =
      sum_or_na(
        Element_production
      ),
    
    .groups =
      "drop"
  ) %>%
  
  group_by(
    Year
  ) %>%
  
  mutate(
    Total_REE_production =
      sum_or_na(
        Element_production
      ),
    
    Element_content =
      Element_production /
      Total_REE_production
  ) %>%
  
  ungroup() %>%
  
  select(
    Year,
    Element,
    Element_content
  )

### REE composition before 1990 -----------------------------------------------

# Before country-level REO production data are available, the global REE
# composition is approximated using the mean elemental composition of
# Australia and the United States.
#
# This composition is held constant from 1850 through 1989.

reo_composition_pre_1990 <- data_REO_content %>%
  filter(
    Country %in% c(
      "Australia",
      "United States"
    )
  ) %>%
  
  select(
    Country,
    all_of(
      REE_ELEMENTS
    )
  ) %>%
  
  mutate(
    across(
      all_of(
        REE_ELEMENTS
      ),
      ~ parse_ree_value(.x)
    )
  ) %>%
  
  summarise(
    across(
      all_of(
        REE_ELEMENTS
      ),
      ~ mean(
        .x,
        na.rm = TRUE
      ) / 100
    )
  ) %>%
  
  crossing(
    Year = 1850:1989
  ) %>%
  
  pivot_longer(
    cols =
      all_of(
        REE_ELEMENTS
      ),
    
    names_to =
      "Element",
    
    values_to =
      "Element_content"
  )

### Combine historical and modern REE compositions ----------------------------

reo_global_composition <- bind_rows(
  reo_composition_pre_1990,
  reo_global_composition_1990_2024
)

### Allocate total REO production among individual REEs -----------------------

# The mining database stores the same total REO production under each REE
# commodity. Each series is therefore multiplied by the corresponding annual
# elemental share to obtain element-specific mining production.

mining_REE <- mining_clean %>%
  filter(
    Commodity %in%
      REE_ELEMENTS
  ) %>%
  
  transmute(
    Element =
      Commodity,
    
    Year,
    Mining_Production,
    Note,
    Finish,
    Commodity
  ) %>%
  
  left_join(
    reo_global_composition,
    by = c(
      "Year",
      "Element"
    )
  ) %>%
  
  mutate(
    FM1 =
      Mining_Production *
      Element_content,
    
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

## 06.6 Platinum-group-metal allocation ------------------------------------
#
# Platinum-group-metal (PGM) production is reported in different forms in the
# mining database:
#
#   1. element-specific production, indicated by a missing `Note`;
#   2. total PGM production, indicated by `Note == "PGM"`; and
#   3. grouped production of the remaining PGMs, indicated by
#      `Note == "other PGM"`.
#
# Element-specific observations are retained unchanged.
#
# Grouped PGM production is allocated among individual elements using fixed
# proportional factors derived from the original reconstruction.


### Allocation factors ---------------------------------------------------------

pgm_allocation_factors <- tribble(
  ~Element, ~Note,       ~factor,
  "Pd",     "PGM",       0.45209,
  "Pt",     "PGM",       0.42477,
  "Ir",     "PGM",       0.015874,
  "Rh",     "PGM",       0.04918,
  "Ru",     "PGM",       0.06462,
  "Ir",     "other PGM", 0.12242,
  "Rh",     "other PGM", 0.39912,
  "Ru",     "other PGM", 0.47188
)

PGM_ELEMENTS <- c(
  "Pd",
  "Pt",
  "Ir",
  "Rh",
  "Ru"
)
#
# This diagnostic helps distinguish element-specific observations from grouped
# "PGM" and "other PGM" records before allocation factors are applied.

pgm_raw_diagnostic <- mining_clean %>%
  filter(
    Commodity %in% PGM_ELEMENTS
  ) %>%
  count(
    Commodity,
    Note,
    sort = TRUE,
    name = "n_rows"
  )

print(
  pgm_raw_diagnostic,
  n = Inf
)

### Allocate grouped PGM production --------------------------------------------

mining_PGM <- mining_clean %>%
  filter(
    Commodity %in% PGM_ELEMENTS
  ) %>%
  
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
    by = c(
      "Element",
      "Note"
    )
  ) %>%
  
  mutate(
    factor = case_when(
      is.na(Note) ~ 1,
      !is.na(factor) ~ factor,
      TRUE ~ NA_real_
    ),
    
    FM1 =
      Mining_Production *
      factor,
    
    conversion_note = case_when(
      is.na(Note) ~
        "Already element-specific",
      
      Note == "PGM" ~
        "Allocated from total PGM production",
      
      Note == "other PGM" ~
        "Allocated from other-PGM production",
      
      TRUE ~
        "Unknown PGM allocation"
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

# Check the total allocation represented by each grouped PGM category.
#
# Factor sums may differ slightly from one because of rounding and because not
# all platinum-group elements are represented in the allocation table.

pgm_factor_diagnostic <- pgm_allocation_factors %>%
  group_by(
    Note
  ) %>%
  summarise(
    factor_sum =
      sum(
        factor
      ),
    
    difference_from_one =
      factor_sum - 1,
    
    .groups =
      "drop"
  )

## 06.7 Final FM1 table ------------------------------------------------------
#
# Combine all mining-conversion pathways into a single harmonized table.
#
# Each `FM1` value represents annual mining-derived elemental production after
# applying the corresponding direct, stoichiometric, ore-grade, or allocation
# conversion.
#
# `bind_rows()` is used instead of repeated `rbind()` calls to simplify the
# assembly and make missing conversion pathways easier to identify.


df_mining <- bind_rows(
  mining_direct,
  mining_Fe,
  mining_Mg,
  mining_Na,
  mining_REE,
  mining_PGM
) %>%
  arrange(
    Element,
    Year
  )


mob_FM1 <- df_mining %>%
  select(
    Element,
    Year,
    FM1
  )

### Export mining-derived elemental production --------------------------------
#
# `mining_FM1_detailed.csv` retains conversion provenance and intermediate
# information.
#
# `mining_FM1.csv` is the compact element-year table used in downstream
# analyses.
#
# `mining_unconverted_commodities.csv` records any source commodities not
# assigned to a conversion pathway.

write_csv(
  df_mining,
  file.path(
    TABULAR_RESULTS_DIR,
    "mining_FM1_detailed.csv"
  )
)

write_csv(
  mob_FM1,
  file.path(
    TABULAR_RESULTS_DIR,
    "mining_FM1.csv"
  )
)

# 07. Fossil-fuel elemental mobilization -----------------------------------
#
# Annual elemental mobilization associated with fossil-fuel consumption is
# calculated as:
#
#   elemental mobilization =
#     fuel consumption × elemental mass fraction
#
# Coal and oil consumption are treated as deterministic after interpolation,
# whereas uncertainty in elemental composition is propagated using Monte
# Carlo simulations.
#
# Fuel consumption is expressed in Tg yr-1 and elemental composition as
# kg element per kg fuel. Their product therefore gives elemental
# mobilization directly in Tg yr-1.

set.seed(SEED)

## 07.1 Coal ----------------------------------------------------------------
#
# Annual coal consumption is combined with the representative coal elemental
# composition derived in Section 04.
#
# Composition uncertainty is represented using normal distributions defined
# by the final central concentration and SD for each element.
#
# For C and S, these final uncertainty estimates already incorporate the
# coal-type composition and production-share Monte Carlo calculations
# performed in Section 04.

coal_mass <- data_coal %>%
  select(
    Year,
    Coal_Tg
  )


### Coal elemental-composition draws ------------------------------------------

coal_composition_draws <- el_coal_final %>%
  select(
    Element,
    med_frac,
    sd_frac
  ) %>%
  
  crossing(
    sim = seq_len(N_SIM)
  ) %>%
  
  group_by(
    Element
  ) %>%
  
  mutate(
    comp_sim = case_when(
      # Missing central estimates remain missing.
      is.na(med_frac) ~
        NA_real_,
      
      # If no usable uncertainty estimate is available, retain the central
      # concentration deterministically.
      is.na(sd_frac) |
        sd_frac <= 0 ~
        med_frac[1],
      
      # Otherwise sample composition from the final estimated distribution.
      TRUE ~
        rnorm(
          n(),
          mean = med_frac[1],
          sd = sd_frac[1]
        )
    ),
    
    # Elemental mass fractions cannot be negative or exceed one.
    comp_sim = clip(
      comp_sim,
      lower = 0,
      upper = 1
    )
  ) %>%
  
  ungroup() %>%
  
  select(
    Element,
    sim,
    comp_sim
  )

### Coal-related elemental mobilization ---------------------------------------

# Tg coal yr-1 × kg element/kg coal = Tg element yr-1.

coal_sims <- coal_mass %>%
  crossing(
    sim = seq_len(N_SIM)
  ) %>%
  
  left_join(
    coal_composition_draws,
    by = "sim"
  ) %>%
  
  transmute(
    Year,
    Element,
    sim,
    
    mob_sim =
      Coal_Tg *
      comp_sim,
    
    flux =
      "FF_coal"
  )


coal_mobilization_summary <- coal_sims %>%
  group_by(
    Year,
    Element
  ) %>%
  
  summarise_mc() %>%
  
  mutate(
    unit = "Tg yr-1",
    flux = "FF_coal"
  )

## 07.2 Oil -----------------------------------------------------------------
#
# Annual oil consumption is combined with the crude-oil elemental composition
# derived in Section 04.
#
# The uncertainty distribution used for each element follows the information
# available during the oil-composition reconstruction:
#
#   - uniform: concentration sampled between oil-specific minimum and maximum;
#   - normal_cv: uncertainty estimated from the reference trace-element CV;
#   - normal_ucc: UCC fallback concentration and uncertainty;
#   - no estimated uncertainty: central concentration retained unchanged.

oil_mass <- data_oil %>%
  select(
    Year,
    Oil_Tg
  )


### Oil elemental-composition draws -------------------------------------------

oil_composition_draws <- el_oil_final %>%
  crossing(
    sim = seq_len(N_SIM)
  ) %>%
  
  group_by(
    Element
  ) %>%
  
  mutate(
    comp_sim = case_when(
      # Missing central estimates remain missing.
      is.na(med_frac) ~
        NA_real_,
      
      # Elements represented by an oil-specific concentration range are
      # sampled from the corresponding uniform distribution.
      distribution == "uniform" &
        !is.na(range_min_frac) &
        !is.na(range_max_frac) ~
        runif(
          n(),
          min = range_min_frac[1],
          max = range_max_frac[1]
        ),
      
      # Single-value estimates completed using the reference CV, and UCC
      # fallback values, are sampled from normal distributions.
      distribution %in% c(
        "normal_cv",
        "normal_ucc"
      ) &
        !is.na(sd_frac) ~
        rnorm(
          n(),
          mean = med_frac[1],
          sd = sd_frac[1]
        ),
      
      # If no uncertainty distribution is available, retain the central value.
      !is.na(med_frac) ~
        med_frac[1],
      
      TRUE ~
        NA_real_
    ),
    
    # Constrain elemental mass fractions to their physically meaningful range.
    comp_sim = clip(
      comp_sim,
      lower = 0,
      upper = 1
    )
  ) %>%
  
  ungroup() %>%
  
  select(
    Element,
    sim,
    comp_sim,
    source,
    distribution
  )


### Oil-related elemental mobilization ----------------------------------------

# Tg oil yr-1 × kg element/kg oil = Tg element yr-1.

oil_sims <- oil_mass %>%
  crossing(
    sim = seq_len(N_SIM)
  ) %>%
  
  left_join(
    oil_composition_draws,
    by = "sim"
  ) %>%
  
  transmute(
    Year,
    Element,
    sim,
    
    mob_sim =
      Oil_Tg *
      comp_sim,
    
    flux =
      "FF_oil"
  )


oil_mobilization_summary <- oil_sims %>%
  group_by(
    Year,
    Element
  ) %>%
  
  summarise_mc() %>%
  
  mutate(
    unit = "Tg yr-1",
    flux = "FF_oil"
  )

### Export fossil-fuel elemental mobilization ---------------------------------

write_csv(
  coal_mobilization_summary,
  file.path(
    TABULAR_RESULTS_DIR,
    "coal_mobilization_summary.csv"
  )
)

write_csv(
  oil_mobilization_summary,
  file.path(
    TABULAR_RESULTS_DIR,
    "oil_mobilization_summary.csv"
  )
)

# 08. Total anthropogenic elemental mobilization ---------------------------
#
# Total anthropogenic elemental mobilization is calculated by combining four
# fluxes:
#
#   FM1      = mining or synthetic elemental production;
#   FF_coal  = coal-related elemental mobilization;
#   FF_oil   = oil-related elemental mobilization;
#   FC       = construction-related elemental mobilization.
#
# Coal, oil, and construction uncertainty is represented by Monte Carlo
# simulations generated in the preceding sections.
#
# FM1 is currently treated as deterministic and is therefore repeated across
# all Monte Carlo iterations before the fluxes are combined.
#
# All fluxes must be expressed in the same mass unit before aggregation.

## 08.1 Prepare fluxes in a common structure --------------------------------

# Repeat deterministic FM1 values across all Monte Carlo iterations so that
# mining can be combined directly with the simulated fossil-fuel and
# construction fluxes.

fm1_sims <- mob_FM1 %>%
  rename(
    mob_sim = FM1
  ) %>%
  crossing(
    sim = seq_len(N_SIM)
  ) %>%
  transmute(
    Element,
    Year,
    sim,
    mob_sim,
    flux = "FM1"
  )


coal_flux_sims <- coal_sims %>%
  select(
    Element,
    Year,
    sim,
    mob_sim
  ) %>%
  mutate(
    flux = "FF_coal"
  )


oil_flux_sims <- oil_sims %>%
  select(
    Element,
    Year,
    sim,
    mob_sim
  ) %>%
  mutate(
    flux = "FF_oil"
  )


construction_flux_sims <- mob_construction_sims %>%
  select(
    Element,
    Year,
    sim,
    mob_sim
  ) %>%
  mutate(
    flux = "FC"
  )


flux_sims <- bind_rows(
  fm1_sims,
  coal_flux_sims,
  oil_flux_sims,
  construction_flux_sims
)

## 08.2 Aggregate fluxes within each simulation -----------------------------
#
# For each element, year, and Monte Carlo iteration, the available
# anthropogenic fluxes are summed.
#
# `sum_or_na()` implements the following rule:
#
#   - if at least one flux is available, sum the available fluxes;
#   - if all fluxes are missing, return NA rather than zero.
#
# This preserves the original partial-sum approach and avoids interpreting a
# complete absence of information as zero mobilization.

mob_total_sims <- flux_sims %>%
  group_by(
    Element,
    Year,
    sim
  ) %>%
  summarise(
    total_sim =
      sum_or_na(
        mob_sim
      ),
    
    # Record how many of the four flux categories contribute to each total.
    n_flux_observed =
      sum(
        !is.na(mob_sim)
      ),
    
    .groups = "drop"
  )

## 08.3 Annual summaries ----------------------------------------------------

mob_total_summary <- mob_total_sims %>%
  rename(
    mob_sim = total_sim
  ) %>%
  group_by(
    Element,
    Year
  ) %>%
  summarise_mc() %>%
  rename(
    total_mean = mob_mean,
    total_median = mob_median,
    total_sd = mob_sd,
    total_q25 = mob_q25,
    total_q75 = mob_q75,
    total_IQR = mob_IQR
  )

# Summarize the number of contributing anthropogenic fluxes.

mob_total_coverage <- mob_total_sims %>%
  group_by(
    Element,
    Year
  ) %>%
  summarise(
    min_fluxes_observed = min(
      n_flux_observed,
      na.rm = TRUE
    ),
    max_fluxes_observed = max(
      n_flux_observed,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

mob_total <- mob_total_summary %>%
  left_join(
    mob_total_coverage,
    by = c(
      "Element",
      "Year"
    )
  )

mob_total <- mob_total %>%
  mutate(
    unit = "Tg yr-1"
  )

mob_total_sims <- mob_total_sims %>%
  mutate(
    unit = "Tg yr-1"
  )

## 08.4 Log-ratio comparisons -----------------------------------------------
#
# Changes in total anthropogenic mobilization between selected years are
# expressed as natural-log response ratios:
#
#   log(M_t2 / M_t1)
#
# where:
#   > 0 indicates an increase in mobilization;
#   = 0 indicates no change;
#   < 0 indicates a decrease in mobilization.
#
# Log ratios are calculated only when mobilization is positive at both
# endpoints.

comparison_years <- c(
  1850,
  1950,
  1960,
  1980,
  2000,
  2020
)

mob_ratio_sims <- mob_total_sims %>%
  filter(
    Year %in% comparison_years
  ) %>%
  
  select(
    Element,
    Year,
    sim,
    total_sim
  ) %>%
  
  pivot_wider(
    names_from = Year,
    values_from = total_sim
  ) %>%
  
  transmute(
    Element,
    sim,
    
    `1850–2020` = if_else(
      `1850` > 0 &
        `2020` > 0,
      log(
        `2020` /
          `1850`
      ),
      NA_real_
    ),
    
    `1950–2020` = if_else(
      `1950` > 0 &
        `2020` > 0,
      log(
        `2020` /
          `1950`
      ),
      NA_real_
    ),
    
    `1960–1980` = if_else(
      `1960` > 0 &
        `1980` > 0,
      log(
        `1980` /
          `1960`
      ),
      NA_real_
    ),
    
    `1980–2000` = if_else(
      `1980` > 0 &
        `2000` > 0,
      log(
        `2000` /
          `1980`
      ),
      NA_real_
    ),
    
    `2000–2020` = if_else(
      `2000` > 0 &
        `2020` > 0,
      log(
        `2020` /
          `2000`
      ),
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
    
    # Re and He are excluded from the comparison analysis.
    !Element %in% c(
      "Re",
      "He"
    )
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
  group_by(
    Element,
    comparison
  ) %>%
  
  summarise(
    q25 =
      quantile(
        log_ratio,
        0.25,
        na.rm = TRUE
      ),
    
    median =
      median(
        log_ratio,
        na.rm = TRUE
      ),
    
    q75 =
      quantile(
        log_ratio,
        0.75,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

## 08.5 Export results -------------------------------------------------------

write_csv(
  mob_total,
  file.path(
    TABULAR_RESULTS_DIR,
    "total_mobilization_summary.csv"
  )
)

write_csv(
  mob_total_sims,
  file.path(
    DRAW_RESULTS_DIR,
    "total_mobilization_draws.csv"
  )
)

write_csv(
  mob_ratio_summary,
  file.path(
    TABULAR_RESULTS_DIR,
    "mobilization_log_ratio_summary.csv"
  )
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














# 10. Selected-year mobilization tables -------------------------------------
#
# Create publication-ready summary tables for selected years for:
#
#   1. mining-derived mobilization (FM1);
#   2. fossil-fuel mobilization (coal + oil); and
#   3. construction-related mobilization.
#
# For each element and selected year, the tables report:
#   - mean annual elemental mobilization; and
#   - standard deviation across Monte Carlo simulations.
#
# All mobilization values are expressed in Tg yr-1.
#
# Mining-derived mobilization (FM1) is treated as deterministic in the current
# workflow. Because the same FM1 value is repeated across simulations, its
# resulting SD is expected to be zero.


## 10.1 Selected years -------------------------------------------------------

selected_years <- c(
  1850,
  1950,
  1960,
  1980,
  2000,
  2020
)

## 10.2 Helper function ------------------------------------------------------

# Summarize Monte Carlo mobilization draws for the selected years.
#
# Every element present in the input dataset is retained for every selected
# year. Missing element-year combinations are therefore explicitly represented
# by NA rather than being omitted from the final tables.

summarise_selected_years <- function(
    draws,
    selected_years
) {
  
  all_elements <- sort(
    unique(
      draws$Element
    )
  )
  
  draws %>%
    filter(
      Year %in% selected_years
    ) %>%
    
    group_by(
      Element,
      Year
    ) %>%
    
    summarise(
      estimate_Tg = mean(
        mob_sim,
        na.rm = TRUE
      ),
      
      sd_Tg = sd(
        mob_sim,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    
    # mean(..., na.rm = TRUE) returns NaN when all values are missing.
    mutate(
      estimate_Tg = if_else(
        is.nan(estimate_Tg),
        NA_real_,
        estimate_Tg
      ),
      
      sd_Tg = if_else(
        is.nan(sd_Tg),
        NA_real_,
        sd_Tg
      )
    ) %>%
    
    complete(
      Element = all_elements,
      Year = selected_years
    ) %>%
    
    arrange(
      Element,
      Year
    )
}

## 10.3 Mining mobilization --------------------------------------------------
#
# FM1 is deterministic in the current workflow. Multiple mining records for
# the same element and year are summed before the selected-year statistics are
# calculated.

mining_selected_draws <- fm1_sims %>%
  group_by(
    Element,
    Year,
    sim
  ) %>%
  summarise(
    mob_sim = sum_or_na(
      mob_sim
    ),
    .groups = "drop"
  )


table_mining_long <- summarise_selected_years(
  draws = mining_selected_draws,
  selected_years = selected_years
) %>%
  
  mutate(
    mobilization = "Mining",
    unit = "Tg yr-1"
  ) %>%
  
  select(
    mobilization,
    Element,
    Year,
    estimate_Tg,
    sd_Tg,
    unit
  )


## 10.4 Fossil-fuel mobilization --------------------------------------------
#
# Coal- and oil-related mobilization are added within each Monte Carlo
# iteration so that their combined uncertainty is propagated into the
# fossil-fuel total.

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
  
  group_by(
    Element,
    Year,
    sim
  ) %>%
  
  summarise(
    mob_sim = sum_or_na(
      mob_sim
    ),
    .groups = "drop"
  )


table_fossil_fuels_long <- summarise_selected_years(
  draws = fossil_fuel_selected_draws,
  selected_years = selected_years
) %>%
  
  mutate(
    mobilization = "Fossil fuels",
    unit = "Tg yr-1"
  ) %>%
  
  select(
    mobilization,
    Element,
    Year,
    estimate_Tg,
    sd_Tg,
    unit
  )

## 10.5 Construction mobilization -------------------------------------------
#
# The total construction-related Monte Carlo draws generated in Section 05
# already combine sand and gravel, carbonate rocks, and clay.

construction_selected_draws <- mob_construction_sims %>%
  select(
    Element,
    Year,
    sim,
    mob_sim
  )


table_construction_long <- summarise_selected_years(
  draws = construction_selected_draws,
  selected_years = selected_years
) %>%
  
  mutate(
    mobilization = "Construction",
    unit = "Tg yr-1"
  ) %>%
  
  select(
    mobilization,
    Element,
    Year,
    estimate_Tg,
    sd_Tg,
    unit
  )

## 10.6 Convert selected-year tables to wide format -------------------------

# Arrange each selected year as paired estimate and SD columns.

make_wide_mobilization_table <- function(
    summary_table
) {
  
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
    
    arrange(
      Element
    )
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

## 10.7 Export selected-year tables -----------------------------------------

write_csv(
  table_mining_wide,
  file.path(
    TABULAR_RESULTS_DIR,
    "selected_years_mining.csv"
  )
)

write_csv(
  table_fossil_fuels_wide,
  file.path(
    TABULAR_RESULTS_DIR,
    "selected_years_fossil_fuels.csv"
  )
)

write_csv(
  table_construction_wide,
  file.path(
    TABULAR_RESULTS_DIR,
    "selected_years_construction.csv"
  )
)


# 11. Final selected-year mobilization table --------------------------------
#
# Create a single wide table containing:
#   - total anthropogenic mobilization;
#   - coal-related mobilization;
#   - oil-related mobilization;
#   - mining-derived mobilization (FM1); and
#   - construction-related mobilization.
#
# Values are reported as mean ± SD when uncertainty is available.
# Deterministic mining values are reported as single values.
# Missing values are represented by "-".
#
# Units: Tg yr-1.


## 11.1 Years included -------------------------------------------------------

final_table_years <- c(
  1850,
  1950,
  2020
)


## Formatting helpers --------------------------------------------------------

# Format individual values using scientific notation only for small,
# non-zero values.

format_number <- function(
    x,
    digits = 2
) {
  
  vapply(
    x,
    function(value) {
      
      if (is.na(value)) {
        return(NA_character_)
      }
      
      if (value != 0 && abs(value) < 0.001) {
        
        format(
          signif(value, digits),
          scientific = TRUE,
          trim = TRUE
        )
        
      } else {
        
        format(
          signif(value, digits),
          scientific = FALSE,
          trim = TRUE
        )
      }
    },
    character(1)
  )
}


# Combine the estimated mean and SD as "mean ± SD".
# Deterministic values (SD = 0) are reported as a single value,
# and missing estimates as "-".

format_mean_sd <- function(
    mean_value,
    sd_value,
    digits = 2
) {
  
  mean_text <- format_number(
    mean_value,
    digits = digits
  )
  
  sd_text <- format_number(
    sd_value,
    digits = digits
  )
  
  case_when(
    is.na(mean_value) ~ "-",
    
    is.na(sd_value) | sd_value == 0 ~
      mean_text,
    
    TRUE ~
      paste0(
        mean_text,
        " ± ",
        sd_text
      )
  )
}
## 11.3 Total anthropogenic mobilization ------------------------------------

table_total <- mob_total %>%
  filter(
    Year %in% final_table_years
  ) %>%
  
  transmute(
    Element,
    Year,
    
    value = format_mean_sd(
      total_mean,
      total_sd,
      digits = 2
    )
  ) %>%
  
  pivot_wider(
    names_from = Year,
    values_from = value,
    names_glue = "M{Year}"
  )

## 11.4 Coal mobilization ----------------------------------------------------

table_coal <- coal_mobilization_summary %>%
  filter(
    Year %in% final_table_years
  ) %>%
  
  transmute(
    Element,
    Year,
    
    value = format_mean_sd(
      mob_mean,
      mob_sd,
      digits = 2
    )
  ) %>%
  
  pivot_wider(
    names_from = Year,
    values_from = value,
    names_glue = "FFF,coal,{Year}"
  )

## 11.5 Oil mobilization -----------------------------------------------------

table_oil <- oil_mobilization_summary %>%
  filter(
    Year %in% final_table_years
  ) %>%
  
  transmute(
    Element,
    Year,
    
    value = format_mean_sd(
      mob_mean,
      mob_sd,
      digits = 2
    )
  ) %>%
  
  complete(
    Element,
    Year = final_table_years
  ) %>%
  
  mutate(
    value = coalesce(
      value,
      "-"
    )
  ) %>%
  
  pivot_wider(
    names_from = Year,
    values_from = value,
    names_glue = "FFF,oil,{Year}"
  )

## 11.6 Mining mobilization --------------------------------------------------

table_mining <- mob_FM1 %>%
  filter(
    Year %in% final_table_years
  ) %>%
  
  group_by(
    Element,
    Year
  ) %>%
  
  summarise(
    FM1 = sum_or_na(
      FM1
    ),
    .groups = "drop"
  ) %>%
  
  transmute(
    Element,
    Year,
    
    value = if_else(
      is.na(FM1),
      "-",
      format_number(
        FM1,
        digits = 2
      )
    )
  ) %>%
  
  complete(
    Element,
    Year = final_table_years,
    fill = list(
      value = "-"
    )
  ) %>%
  
  pivot_wider(
    names_from = Year,
    values_from = value,
    names_glue = "FM,{Year}"
  )

## 11.7 Construction mobilization -------------------------------------------

table_construction <- elem_mob_aggregate %>%
  select(
    Element,
    Year
  )

table_construction <- mob_construction_total %>%
  filter(
    Year %in% final_table_years
  ) %>%
  
  transmute(
    Element,
    Year,
    
    value = format_mean_sd(
      mob_mean,
      mob_sd,
      digits = 2
    )
  ) %>%
  
  pivot_wider(
    names_from = Year,
    values_from = value,
    names_glue = "FC,{Year}"
  )

## 11.8 Combine all mobilization components ---------------------------------

table_mobilization_final <- table_total %>%
  full_join(
    table_coal,
    by = "Element"
  ) %>%
  
  full_join(
    table_oil,
    by = "Element"
  ) %>%
  
  full_join(
    table_mining,
    by = "Element"
  ) %>%
  
  full_join(
    table_construction,
    by = "Element"
  )

table_mobilization_final <- table_mobilization_final %>%
  mutate(
    across(
      -Element,
      ~ replace_na(
        .x,
        "-"
      )
    )
  )

table_mobilization_final <- table_mobilization_final %>%
  select(
    Element,
    
    M1850,
    M1950,
    M2020,
    
    `FFF,coal,1850`,
    `FFF,coal,1950`,
    `FFF,coal,2020`,
    
    `FFF,oil,1950`,
    `FFF,oil,2020`,
    
    `FM,1850`,
    `FM,1950`,
    `FM,2020`,
    
    `FC,1850`,
    `FC,1950`,
    `FC,2020`
  )

element_order <- mob_total %>%
  filter(
    Year == 2020
  ) %>%
  arrange(
    desc(total_mean)
  ) %>%
  pull(
    Element
  )


table_mobilization_final <- table_mobilization_final %>%
  mutate(
    Element = factor(
      Element,
      levels = element_order
    )
  ) %>%
  
  arrange(
    Element
  )

write_csv(
  table_mobilization_final,
  file.path(
    TABULAR_RESULTS_DIR,
    "selected_year_mobilization_table.csv"
  )
)
