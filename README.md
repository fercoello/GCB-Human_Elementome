# Anthropogenic Effects on the Elementome

Code associated with the study:

Peñuelas, J., Coello, F., de la Casa, J., Nogué, S., Fernandez-Martinez, M., Ogaya, R., & Sardans, J. (in press). *Anthropogenic effects on elementomes
across organisms, ecosystems, and the biosphere*. Global Change Biology.

This repository contains the R workflow used to reconstruct anthropogenic elemental mobilization associated with mining, construction materials, coal consumption, and oil consumption (Figure 3).

The datasets required to run the analyses are archived separately on Zenodo.

## Repository contents

The workflow consists of three R scripts:

- `00_Prepare_Construction_Data.R`  
  Reconstructs annual global mass series for construction materials, including
  aggregates, clay, and carbonate rocks. The workflow integrates historical
  material-production data, population data, railway statistics, and
  literature-derived engineering assumptions.

- `00_Prepare_Elemental_Composition.R`  
  Prepares limestone and dolomite elemental-composition summaries from the
  whole-rock geochemical database of Gard et al. (2019), including filtering,
  summary-statistic calculation, and oxide-to-element conversion.

- `01_Elemental_Mobilization.R`  
  Main analysis script. Reconstructs annual anthropogenic elemental
  mobilization from mining, construction materials, coal, and oil and propagates
  uncertainty through Monte Carlo simulation.

## Analytical framework

Total anthropogenic elemental mobilization is reconstructed from four main sources:

- mining-derived elemental production (`FM1`);
- coal-related mobilization (`FF_coal`);
- oil-related mobilization (`FF_oil`); and
- construction-related mobilization (`FC`).

Mining-derived mobilization is treated as deterministic. Uncertainty in
fossil-fuel elemental composition and construction-material mass and
composition is propagated using Monte Carlo simulation.

Annual elemental mobilization is expressed in **Tg yr⁻¹**. Elemental
compositions are standardized internally as mass fractions (**kg kg⁻¹**).

## Data availability

The datasets supporting this workflow, including the datasets underlying Figure 3 of the associated publication, are archived on Zenodo:

**Zenodo:**  
https://doi.org/10.5281/zenodo.23214154

The deposited materials include:

- harmonized historical mining and ore-production data;
- elemental-composition data for construction materials and fossil fuels; and
- supporting railway-length and railway-gauge data used in the construction-material reconstruction.

Source references, temporal coverage, assumptions, and reconstruction procedures are documented in the Supplementary Methods of the associated publication and in the README documentation included with the
Zenodo datasets.

## Expected data structure

After downloading the archived datasets, the files should be arranged so that the scripts can access them using the relative paths specified in the code.

## Reproducibility
All scripts use relative file paths based on here::here(). The repository should therefore be opened and run from its project root.
A fixed random-number seed is used for Monte Carlo simulations to facilitate reproducibility.
Methodological details, including source selection, interpolation and extrapolation procedures, elemental-concentration assumptions, conversion
factors, and uncertainty propagation, are described in the Supplementary Methods of the associated publication.

## Citation
If using this code or the associated datasets, please cite the associated
publication:
Peñuelas, J., Coello, F., de la Casa, J., Nogué, S., Fernandez-Martinez, M., Ogaya, R., & Sardans, J. (in press). Anthropogenic effects on elementomesacross organisms, ecosystems, and the biosphere (in press). Global Change Biology.

The archived datasets can additionally be cited using the Zenodo DOI:
https://doi.org/10.5281/zenodo.23214154

## Contact
Fernando Coello
CREAF
f.coello@creaf.uab.cat
