# Mapping priority areas for Wolbachia-based malaria control across Sub-Saharan Africa

## Overview

This repository contains all R code and input data for the spatial analysis described in the paper. 

We develop a geospatial framework to map the eligibility of sub-Saharan African locations for _Wolbachia_-based malaria vector control across two demonstration trials and three programmatic use cases, covering five _Anopheles_ species.

## Repository contents

### R scripts

The /R folder contains six Shiny applications. Each app includes a Map tab, an Impact and Coverage tab, and a Sensitivity tab allowing users to explore how results change under different threshold assumptions.

* Demo_trial_1a.R - Demonstration trial 1a: entomological feasibility
* Demo_trial_1b.R - Demonstration trial 1b: epidemiological evaluation
* Use_case_2.R - Use case 2: routine implementation
* Use_case_3.R - Use case 3: elimination support
* Use_case_4.R - Use case 4: prevention of An. stephensi re-emergence
* Overlay_map.R - Spatial overlap between Trial 1a and Trial 1b eligibility maps

### Data files

The /data folder contains all input raster files used in the analysis. See the Data sources section below for details of each file and its source.

## Data sources 

| File | Source | Description |
| --- | --- | --- |
| An. arabiensis.tif | [Malaria Atlas Project](https://malariaatlas.org/) | _An. arabiensis_ occurrence probability (5km) |
| An. coluzzii.tif | [Malaria Atlas Project](https://malariaatlas.org/) | _An. coluzzii_ occurrence probability (5km) |
| An. moucheti.tif | [Malaria Atlas Project](https://malariaatlas.org/) | _An. moucheti_ occurrence probability (5km) |
| An. gambiae.tif | [Malaria Atlas Project](https://malariaatlas.org/) | _An. gambiae s.s._ occurrence probability (5km) |
| An. stephensi.tif | [Malaria Atlas Project](https://malariaatlas.org/) | _An. stephensi_ occurrence probability (5km) |
| PfPR_mean.tif | [Malaria Atlas Project](https://malariaatlas.org/) | Mean PfPR 2000-2020 (5km) |
| ITN_use_rate.tif | [Malaria Atlas Project](https://malariaatlas.org/) | ITN use rate (5km) |
| Pf_incidence_mean_2000.tif | [Malaria Atlas Project](https://malariaatlas.org/) | Mean Pf incidence rate 2000 (5km) |
| POP_MEAN_2000_2020_5km.tif | [World Pop](https://hub.worldpop.org/geodata/listing?id=64) | Mean population count 2000–2020 aggregated to 5km |
| GHS_SMOD_E2025_GLOBE_R2023A_54009_1000_V2_0.tif | [Global Human Settlement Layer](https://human-settlement.emergency.copernicus.eu/ghs_smod2023.php) | Urban/rural typology 2025 (1km) |

## Analytical framework

Each pixel across sub-Saharan Africa is assigned a priority tier based on an additive scoring system. Each parameter is scored 1 (optimal) to 3 (suboptimal) and scores are summed to determine the final tier (Tier 1 = highest priority).

| Use case | Parameter scored | Tier 1 creiteria |
| --- | --- | --- |
| []{rowspan=3} 'Trial 1a' | PfPR | Low |
|          | Dominance | High |
|          | ITN | Low |
