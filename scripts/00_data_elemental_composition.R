# ==============================================================================
# Project: GCB — Anthropogenic Effects on the Elementome
#
# Script: 00_Prepare_Elemental_Composition.R
#
# Purpose:
#   Prepare elemental-composition summaries for limestone and dolomite from the
#   whole-rock geochemical database compiled by Gard et al. (2019).
#
# Analytical framework:
#   The script:
#     1. extracts limestone and dolomite records;
#     2. removes redundant or unsuitable geochemical variables;
#     3. retains samples with CaO concentrations >= 10 wt%;
#     4. removes negative and non-finite values;
#     5. calculates element-specific summary statistics; and
#     6. converts major-oxide concentrations to elemental concentrations.
#
#   These material-specific summaries are subsequently used in the main
#   elemental-mobilization analysis. Where fewer than 50 observations are
#   available for an element, the downstream analysis uses Upper Continental
#   Crust concentrations as fallback estimates.
#
# Units:
#   Source concentrations retain their original reported units during
#   preprocessing. Major oxides are converted to elemental concentrations
#   using stoichiometric mass fractions. Final conversion to kg/kg is performed
#   in `01_Elemental_Mobilization.R`.
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
# Last updated:
#   2026-10-07
#
# R version:
#   R 4.5.1
#
# Main dependencies:
#   tidyverse, readr, here
#
# Inputs:
#   data/data_rocks/complete/complete.csv
#
# Outputs:
#   el_limestone.csv
#   el_dolomite.csv
#
# Downstream use:
#   Limestone and dolomite summaries are combined in
#   `01_Elemental_Mobilization.R` assuming that carbonate-rock extraction is
#   represented by 85% limestone and 15% dolomite.
#
# Reproducibility:
#   Run the script from the root directory of the RStudio project. File paths
#   are constructed relative to the project root using here::here().
#
# Data availability:
#   The Gard et al. (2019) source database is third-party material and is not
#   redistributed with this repository unless permitted by the original source.
#   Derived limestone and dolomite summaries used by the analysis are provided
#   as part of the harmonized elemental-composition inputs.
#
# License:
#   Code: MIT License
#   Data and documentation: CC BY 4.0
#
# Citation:
#   See repository citation information.
# ==============================================================================


# 00. Setup ----------------------------------------------------------------

library(tidyverse)
library(here)


# Project paths -------------------------------------------------------------

ROCK_DATA_FILE <- here(
  "data",
  "data_rocks",
  "complete",
  "complete.csv"
)

ELEMENT_CONCENTRATION_DIR <- here(
  "data",
  "elemental_concentration"
)

INTERMEDIATE_DIR <- file.path(
  ELEMENT_CONCENTRATION_DIR,
  "intermediate"
)

dir.create(
  INTERMEDIATE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# Carbonate-rock assumptions ------------------------------------------------

# Minimum CaO concentration retained in the limestone and dolomite datasets.
#
# This compositional filter was part of the original workflow and should be
# reported explicitly in the Supplementary Methods.
MIN_CAO_WT_PCT <- 10


# Relative contribution of limestone and dolomite to global carbonate-rock
# extraction. These values are applied later in the main analysis.
LIMESTONE_SHARE <- 0.85
DOLOMITE_SHARE <- 0.15

# Minimum number of observations required for a material-specific carbonate
# estimate before falling back to the UCC composition in the main analysis.
MIN_CARBONATE_N <- 50L


# Oxide-to-element conversion factors --------------------------------------

# Stoichiometric mass fraction of the target element in each major oxide.
# These factors convert oxide concentrations in wt% to elemental concentrations
# in wt%.

oxide_to_element <- c(
  sio2      = 0.4674367,
  tio2      = 0.5993489,
  al2o3     = 0.5292612,
  cr2o3     = 0.6842078,
  fe2o3_tot = 0.6994308,
  mgo       = 0.6030419,
  cao       = 0.7146959,
  mno       = 0.7744618,
  k2o       = 0.8301513,
  na2o      = 0.7418625,
  p2o5      = 0.4364271,
  co2       = 0.2729214,
  so3       = 0.4005021
)


# Element associated with each major-oxide variable.

element_name <- c(
  sio2      = "Si",
  tio2      = "Ti",
  al2o3     = "Al",
  cr2o3     = "Cr",
  fe2o3_tot = "Fe",
  mgo       = "Mg",
  cao       = "Ca",
  mno       = "Mn",
  k2o       = "K",
  na2o      = "Na",
  p2o5      = "P",
  co2       = "C",
  so3       = "S"
)


# Geochemical variables excluded from the summaries -------------------------

cols_to_remove <- c(
  
  # Alternative or redundant Fe representations.
  "fe2o3",
  "feo",
  "feo_tot",
  
  # Water and loss-on-ignition variables.
  "h2o_plus",
  "h2o_minus",
  "h2o_tot",
  "loi",
  
  # Elemental ppm measurements duplicating major elements represented by
  # oxide variables.
  "al_ppm",
  "ca_ppm",
  "mg_ppm",
  "na_ppm",
  "k_ppm",
  "si_ppm",
  "ti_ppm",
  "fe_ppm",
  "mn_ppm",
  
  # Rare radionuclides excluded from the present analysis.
  "pa_ppm",
  "pm_ppm"
)


# 01. Load whole-rock geochemical database ---------------------------------

# Whole-rock geochemical database compiled by Gard et al. (2019).

rock_composition_db <- read_csv(
  ROCK_DATA_FILE,
  show_col_types = FALSE
)


# Optional diagnostic: number of records classified as each rock type.
#
# Expected values in the source database:
#   limestone = 5756
#   dolomite  = 1302

carbonate_sample_counts <- rock_composition_db %>%
  filter(
    rock_name %in% c(
      "limestone",
      "dolomite"
    )
  ) %>%
  count(
    rock_name,
    name = "n_samples"
  )

carbonate_sample_counts


# 02. Limestone elemental composition --------------------------------------

# Extract samples classified as limestone in the source database.
#
# Columns 77:180 correspond to the geochemical-composition block used in the
# original analysis. Positional selection is retained here to reproduce the
# published workflow exactly.

limestone_db <- rock_composition_db %>%
  filter(
    rock_name == "limestone"
  ) %>%
  select(
    77:180
  )


# Remove redundant or unsuitable geochemical variables.

limestone_clean <- limestone_db %>%
  select(
    -any_of(cols_to_remove)
  )


# Retain observations with CaO >= 10 wt%.
#
# Negative, NaN, and infinite values are treated as missing independently for
# each geochemical variable.

limestone_clean <- limestone_clean %>%
  filter(
    !is.na(cao),
    cao >= MIN_CAO_WT_PCT
  ) %>%
  mutate(
    across(
      everything(),
      ~ ifelse(
        . < 0 |
          is.nan(.) |
          is.infinite(.),
        NA,
        .
      )
    )
  )


# Calculate summary statistics independently for each retained geochemical
# variable.
#
# The number of observations therefore differs among elements according to
# element-specific data availability.

stats_limestone <- data.frame(
  variable = names(limestone_clean),
  
  mean = sapply(
    limestone_clean,
    function(x) mean(
      x,
      na.rm = TRUE
    )
  ),
  
  median = sapply(
    limestone_clean,
    function(x) median(
      x,
      na.rm = TRUE
    )
  ),
  
  q25 = sapply(
    limestone_clean,
    function(x) quantile(
      x,
      0.25,
      na.rm = TRUE
    )
  ),
  
  q75 = sapply(
    limestone_clean,
    function(x) quantile(
      x,
      0.75,
      na.rm = TRUE
    )
  ),
  
  sd = sapply(
    limestone_clean,
    function(x) sd(
      x,
      na.rm = TRUE
    )
  ),
  
  n = sapply(
    limestone_clean,
    function(x) sum(
      !is.na(x)
    )
  )
) %>%
  
  mutate(
    
    # Variables ending in "_ppm" are treated as elemental concentrations
    # reported in ppm. Other retained major-element variables are initially
    # treated as oxide concentrations in wt%.
    unit = if_else(
      str_detect(
        variable,
        regex(
          "_ppm$",
          ignore_case = TRUE
        )
      ),
      "ppm",
      "wt%"
    ),
    
    # Assign elemental symbols.
    #
    # Major oxides are mapped using `element_name`.
    # Trace-element variables ending in "_ppm" are converted from, for example,
    # "zn_ppm" to "Zn".
    element = case_when(
      variable %in% names(element_name) ~
        unname(
          element_name[variable]
        ),
      
      str_detect(
        variable,
        regex(
          "_ppm$",
          ignore_case = TRUE
        )
      ) ~
        variable %>%
        str_remove(
          regex(
            "_ppm$",
            ignore_case = TRUE
          )
        ) %>%
        str_to_lower() %>%
        str_to_title(),
      
      TRUE ~ NA_character_
    ),
    
    # Retrieve oxide-to-element conversion factor when available.
    factor = if_else(
      variable %in% names(oxide_to_element),
      unname(
        oxide_to_element[variable]
      ),
      NA_real_
    )
  ) %>%
  
  # Convert major-oxide statistics from oxide wt% to elemental wt%.
  #
  # Trace-element concentrations reported in ppm are left unchanged.
  mutate(
    mean = if_else(
      !is.na(factor) &
        unit == "wt%",
      mean * factor,
      mean
    ),
    
    median = if_else(
      !is.na(factor) &
        unit == "wt%",
      median * factor,
      median
    ),
    
    sd = if_else(
      !is.na(factor) &
        unit == "wt%",
      sd * factor,
      sd
    ),
    
    q25 = if_else(
      !is.na(factor) &
        unit == "wt%",
      q25 * factor,
      q25
    ),
    
    q75 = if_else(
      !is.na(factor) &
        unit == "wt%",
      q75 * factor,
      q75
    ),
    
    unit = if_else(
      !is.na(factor) &
        unit == "wt%",
      "wt%_element",
      unit
    )
  ) %>%
  
  # Retain only geochemical variables that can be assigned to an element.
  select(
    element,
    mean,
    sd,
    median,
    q25,
    q75,
    unit,
    n
  ) %>%
  drop_na(
    element
  )


# 03. Dolomite elemental composition ---------------------------------------

# Repeat the same workflow for samples classified as dolomite.

dolomite_db <- rock_composition_db %>%
  filter(
    rock_name == "dolomite"
  ) %>%
  select(
    77:180
  )


# Remove the same redundant or unsuitable geochemical variables used for
# limestone.

dolomite_clean <- dolomite_db %>%
  select(
    -any_of(cols_to_remove)
  )


# Apply the same CaO and data-quality filters.

dolomite_clean <- dolomite_clean %>%
  filter(
    !is.na(cao),
    cao >= MIN_CAO_WT_PCT
  ) %>%
  mutate(
    across(
      everything(),
      ~ ifelse(
        . < 0 |
          is.nan(.) |
          is.infinite(.),
        NA,
        .
      )
    )
  )


# Calculate element-wise summary statistics.

stats_dolomite <- data.frame(
  variable = names(dolomite_clean),
  
  mean = sapply(
    dolomite_clean,
    function(x) mean(
      x,
      na.rm = TRUE
    )
  ),
  
  median = sapply(
    dolomite_clean,
    function(x) median(
      x,
      na.rm = TRUE
    )
  ),
  
  q25 = sapply(
    dolomite_clean,
    function(x) quantile(
      x,
      0.25,
      na.rm = TRUE
    )
  ),
  
  q75 = sapply(
    dolomite_clean,
    function(x) quantile(
      x,
      0.75,
      na.rm = TRUE
    )
  ),
  
  sd = sapply(
    dolomite_clean,
    function(x) sd(
      x,
      na.rm = TRUE
    )
  ),
  
  n = sapply(
    dolomite_clean,
    function(x) sum(
      !is.na(x)
    )
  )
) %>%
  
  mutate(
    
    # Assign original units.
    unit = if_else(
      str_detect(
        variable,
        regex(
          "_ppm$",
          ignore_case = TRUE
        )
      ),
      "ppm",
      "wt%"
    ),
    
    # Assign elemental symbols.
    element = case_when(
      variable %in% names(element_name) ~
        unname(
          element_name[variable]
        ),
      
      str_detect(
        variable,
        regex(
          "_ppm$",
          ignore_case = TRUE
        )
      ) ~
        variable %>%
        str_remove(
          regex(
            "_ppm$",
            ignore_case = TRUE
          )
        ) %>%
        str_to_lower() %>%
        str_to_title(),
      
      TRUE ~ NA_character_
    ),
    
    # Retrieve oxide-to-element conversion factor.
    factor = if_else(
      variable %in% names(oxide_to_element),
      unname(
        oxide_to_element[variable]
      ),
      NA_real_
    )
  ) %>%
  
  # Convert major-oxide statistics to elemental wt%.
  mutate(
    mean = if_else(
      !is.na(factor) &
        unit == "wt%",
      mean * factor,
      mean
    ),
    
    median = if_else(
      !is.na(factor) &
        unit == "wt%",
      median * factor,
      median
    ),
    
    sd = if_else(
      !is.na(factor) &
        unit == "wt%",
      sd * factor,
      sd
    ),
    
    q25 = if_else(
      !is.na(factor) &
        unit == "wt%",
      q25 * factor,
      q25
    ),
    
    q75 = if_else(
      !is.na(factor) &
        unit == "wt%",
      q75 * factor,
      q75
    ),
    
    unit = if_else(
      !is.na(factor) &
        unit == "wt%",
      "wt%_element",
      unit
    )
  ) %>%
  
  select(
    element,
    mean,
    sd,
    median,
    q25,
    q75,
    unit,
    n
  ) %>%
  drop_na(
    element
  )


# 04. Validation checks -----------------------------------------------------

# Check that the cleaned workflow reproduces known values from the original
# analysis. For example, elemental Si in limestone should be based on
# 2161 observations with the current source database.

stats_limestone %>%
  filter(
    element == "Si"
  )


# Compare the number of observations retained after the CaO filter.

carbonate_filtered_counts <- tibble(
  rock_type = c(
    "limestone",
    "dolomite"
  ),
  
  n_after_filter = c(
    nrow(limestone_clean),
    nrow(dolomite_clean)
  )
)

carbonate_filtered_counts


# 05. Export intermediate carbonate summaries ------------------------------

write_csv(
  stats_limestone,
  file.path(
    INTERMEDIATE_DIR,
    "el_limestone.csv"
  )
)

write_csv(
  stats_dolomite,
  file.path(
    INTERMEDIATE_DIR,
    "el_dolomite.csv"
  )
)
