# Which district conditions move with adult obesity in Korea? Walking, stress, drinking, and education over 13 survey years

Analysis code for the manuscript

> Kim S. *Which district conditions move with adult obesity in Korea? Walking, stress, drinking, and education over 13 survey years* Manuscript submitted for publication.

Archived version: Zenodo, https://doi.org/10.5281/zenodo.24518235 (v1.1.0). All versions: https://doi.org/10.5281/zenodo.24390860.

## Contents

Code only. The processed district-year data (master table, analytic panel, survey outcomes, result tables, adjacency graph and boundary file) are not included; they are available from the author on reasonable request. All source statistics are public (see Data sources).

| Path | What it is |
|---|---|
| `01_스크립트/` | Analysis scripts (below) |
| `04_분석/OBES/` | Covariates with expected directions, controls, negative controls and model settings (`sets_OBES_std3.json`) and the model chain (`run_fh_chain.sh`) |
| `04_분석/OBES/MS/build/` | `numbers_obes.py` (reads the result tables and computes every reported number), `figs_obes.py` (figures), `build_refs_obes.py` and `refs_obes.json` (reference metadata from PubMed and Crossref), and document helpers (`docx_engine.py`, `booktabs.py`, `tablefit.py`). `figs_obes_v2.py` draws Figure 2 at print size. The manuscript text generator is not included. |
| `04_분석/OBES/results_plus1_check/` | Code only for the +1 offset sensitivity of the green-zone covariate: `make_green_variants.py` (alternative transforms), `refit_green_plus1.R` (main model refitted), `compare.py` (comparison) and `run.sh`. Result files are not included. |

Folder names are in Korean because the source tables are Korean; each script has a short header describing what it does.

## Data sources

All inputs are public, aggregated statistics. No individual-level data are included.

- Outcome: Korea Community Health Survey (Korea Disease Control and Prevention Agency), district prevalence of obesity (body mass index of 25 kg/m² or higher from self-reported height and weight, adults aged 19 years or older), crude and age-sex-standardized (2005 projected population) with standard errors and numbers of respondents, 2011-2024 (KOSIS table DT_H_OBE_OBE; 2019 not surveyed). Measured obesity, 2018-2019 (DT_117075_MEA_OBE), for comparison.
- Covariates: Community Health Survey (walking, moderate-or-vigorous physical activity, high-risk drinking, current smoking, perceived stress), national public sports facility statistics, village sports facility statistics, urban planning statistics (green zone area), motor vehicle registration, the national local-licensing register (general restaurants and instant-food manufacturing-retail outlets), population and housing census (tertiary education), basic livelihood security statistics (invalid 2021-2023 values in 11 cities set to missing and interpolated), and KOSIS local indicators (controls and negative controls).

## Pipeline

Set `ST_DATA_ROOT` to local copies of the public source files and `ST_GRAPH` to the district contact graph.

| Step | Script | Output |
|---|---|---|
| 1 | `99c_extra_canonical_pool_obes.py` | Basic livelihood security recipients, tertiary education, restaurant and instant-food outlet densities, by district code |
| 2 | `99O_obes_fh_prep.py` | Survey outcomes on the logit scale with delta-method sampling variances, effective sample sizes, and measured obesity |
| 3 | `80_build_disease_panel_v.py OBES` (`ST_VER=std3`) | Transformed and standardized covariates, filling and variance inflation factors |
| 4 | `81fh_obes.R OBES` | Spatiotemporal Fay-Herriot models (R-INLA): screening indices, set comparison, main model, negative controls, windows, subset models, Knorr-Held II-IV, between/within decomposition, lags, leave-one-province-out, alternative likelihoods, graph and prior, measured-obesity comparison, restricted spatial regression, regional slopes, influential observations, checks of associations opposite to the expected direction |
| 5 | `99O4_obes_maps_v2.py` (`99O4_obes_maps.py` earlier layout) | Figure 1 district maps at print size |
| 6 | `04_분석/OBES/MS/build/figs_obes_v2.py` | Figure 2 |
| 7 | `04_분석/OBES/results_plus1_check/run.sh` | +1 offset sensitivity of the green-zone covariate |

## Software

R 4.6.0 with R-INLA 25.10.19, sf and spdep; Python 3.9.6 with pandas, numpy, geopandas and matplotlib.

## Versions

- v1.1.0: print-size figure scripts (Figure 1 maps with scale bar, Figure 2 with an estimate column), +1 offset sensitivity code for the green-zone covariate, updated reference metadata, document helper renamed to `docx_engine.py`.
- v1.0.0: first release.

## License

MIT (see `LICENSE`).
