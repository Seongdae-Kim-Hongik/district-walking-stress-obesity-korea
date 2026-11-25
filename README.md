# Do walking, stress, drinking, and education distinguish districts with more adult obesity?

Analysis code for the manuscript

> Kim S. *Do walking, stress, drinking, and education distinguish districts with more adult obesity?* Manuscript submitted for publication.

Archived version: Zenodo, https://doi.org/10.5281/zenodo.24535579 (v1.2.1). All versions: https://doi.org/10.5281/zenodo.24390860.

## Contents

Code only. The processed district-year data (master table, district-year SQLite, analytic panel, survey outcomes, result tables, adjacency graph, district code tables and boundary file) are not included; they are available from the author on reasonable request. All source statistics are public (see Data sources).

| Path | What it is |
|---|---|
| `01_스크립트/` | Data preparation, panel building, model chain and figure scripts (see Pipeline) |
| `04_분석/OBES/panel_gate/` | Pre-model data gate and panel correction: cell-by-cell comparison with the district-year SQLite, year-jump and zero scans, village sports facility source-table checks and corrections, corrected panel build, denominator and boundary-file audit, zero/year/fill audit, general data and join gates (`data_gate.py`, `panel_gate.py`, `sigungu_key.py`, `70_year_anomaly_gate.py`) and KOSIS population fetchers |
| `04_분석/OBES/run/` | `make_jobs.py` (job lists for the full chain) and `runner.sh` (job runner with load limit, per-job exit codes and time limits) |
| `04_분석/OBES/extra/` | Post-processing of the additional analyses (`02_G_numeric.py`, `03_summarize.py`, `04_tables_md.py`) and post hoc tables (`05_posthoc.py`); no new fits |
| `04_분석/OBES/extra2/` | Round-2 additional analyses: inputs and checks without fits (`e2_00_prep.py`), post hoc standard sets s577, s2 and s8 (`e2_05_make_sets.py`), data-gate threshold alternatives (`e2_10_scan_thresholds.py`, `e2_11_village_rule_thresholds.py`, `e2_12_prep_a205.py`), fits J1-J7 (`e2_fits.R`, run with `make_jobs.py` and `runner.sh`), selection-corrected 2018 split test (`e2_20_split_perm.py`), maps (`e2_30_maps_figs.py`) and the combination of all sensitivity specifications (`e2_40_combine.py`) |
| `04_분석/OBES/sets_OBES_std3.json` | Covariates with expected directions, controls, negative controls and model settings |
| `04_분석/OBES/MS/build/` | `numbers_obes.py` (reads the main-chain result tables and writes `numbers_obes.json`), `figs_common.py` (numbers, half-up rounding and labels shared by the figure scripts), `figs_obes_print.py` (Figure 2 at print size), `figs_si_print.py` (Supplementary Figures S1-S4 from the `extra2/` tables), `build_refs_obes.py` and `refs_obes.json` (reference metadata from PubMed and Crossref; `REFS_SOURCE` = manuscript source with the citation keys), and document helpers (`docx_engine.py`, `booktabs.py`, `tablefit.py`). The manuscript text generator is not included. Numbers of the main chain come from `numbers_obes.json`; numbers of the additional analyses (`extra/`, `extra2/`) are read directly from the tables those scripts write. |
| `04_분석/OBES/results_plus1_check/` | Code only for the +1 offset sensitivity of the green-zone covariate: `make_green_variants.py` (alternative transforms), `compare.py` (comparison) and `run.sh`. Result files are not included. |

Scripts whose first line says `SUPERSEDED in v1.2.0` (`99O4_obes_maps.py`, `figs_obes.py`, `run_fh_chain.sh`, `results_plus1_check/refit_green_plus1.R`) are kept for reference only; the replacement is named in that line.

Folder names are in Korean because the source tables are Korean; each script has a short header describing what it does.

## Data sources

All inputs are public, aggregated statistics. No individual-level data are included.

- Outcome: Korea Community Health Survey (Korea Disease Control and Prevention Agency), district prevalence of obesity (body mass index of 25 kg/m² or higher from self-reported height and weight, adults aged 19 years or older), crude and age-sex-standardized (2005 projected population) with standard errors and numbers of respondents, 2011-2024 (KOSIS table DT_H_OBE_OBE; 2019 not surveyed). Measured obesity, 2018-2019 (DT_117075_MEA_OBE), for comparison.
- Covariates: Community Health Survey (walking, moderate-or-vigorous physical activity, high-risk drinking, current smoking, perceived stress), national public sports facility statistics, village sports facility statistics (Ministry of Culture, Sports and Tourism yearly tables), urban planning statistics (green zone area), motor vehicle registration, the national local-licensing register (general restaurants and instant-food manufacturing-retail outlets), population and housing census (tertiary education), basic livelihood security statistics (official recipient counts and population), and KOSIS local indicators (controls and negative controls).
- Denominator check: KOSIS resident registration population by district (DT_1B04005N), fetched with `kosis_fetch_pop.py` (set `KOSIS_API_KEY`).

## Pipeline

Environment variables: `ST_DATA_ROOT` (local copies of the public source files), `ST_SQLITE` (district-year SQLite built from them, opened read-only), `PANEL_WORKDIR` (work folder of the data gate), `SIGUNGU_KEY_DIR` (district code tables), `ST_SHP` (boundary file, default `data/boundary/final.shp`), `ST_GRAPH` (district contact graph).

| Step | Script | Output |
|---|---|---|
| 1 | `01_스크립트/99c_extra_canonical_pool_obes.py` | Tertiary education, restaurant and instant-food outlet densities, by district code |
| 2 | `01_스크립트/99O_obes_fh_prep.py` | Survey outcomes on the logit scale with delta-method sampling variances, effective sample sizes, and measured obesity |
| 3 | `01_스크립트/80_build_disease_panel_v.py OBES` (`ST_VER=std3`) | First-release panel (transformed and standardized covariates, filling and variance inflation factors); copied to `PANEL_WORKDIR/input_first/` |
| 4 | `panel_gate/01_sqlite_compare.py`, `02_year_jump_scan.py` | Cell-by-cell comparison with the SQLite; year-jump and zero-cell flags |
| 5 | `panel_gate/03_village_rowshift_check.py`, `04_village_permutation.py`, `05_village_count_swap.py`, `06_village_corrections.py` | Village sports facility source-table defects (row shifts, rotated rows, swapped counts, double counting) and the correction record |
| 6 | `panel_gate/07_build_panel.py` | Corrected panel `out/panel_OBES_std3.csv` with `panel_META`, `panel_INFO` and the fill ledger `fill_cells_OBES_std3.csv`; the uncorrected rebuild is checked against the first-release panel cell by cell first |
| 7 | `panel_gate/08_denominator_shape_audit.py`, `09_zero_year_fill_audit.py`, `data_gate.py`, `panel_gate.py` | Denominator, boundary-file, zero, year-anomaly and fill audits of the corrected panel; the gates must pass before any model is fitted |
| 8 | copy `out/panel_OBES_std3.csv`, `panel_META_OBES_std3.csv`, `panel_INFO_OBES_std3.json` to `04_분석/OBES/` of the run root | Model input |
| 9 | `01_스크립트/82_extra_inputs.py <run root> <fill ledger> 04_분석/OBES/results_extra`; `results_plus1_check/make_green_variants.py <run root>` | Inputs of the additional analyses and of the +1 sensitivity |
| 10 | `run/make_jobs.py <run root> <boundary file> <graph> <jobs folder>`, then `RSLOT=01_스크립트/rslot.sh run/runner.sh <log folder> A.jobs B.jobs C.jobs` | Phase A: screening, set comparison and main model (`81fh_obes.R` S0, S2, S3) and the S7 specification list. Phase B: S7 multiverse in four blocks; S4-S6, S8-S14, S16, S18, S19 (negative controls, windows, subsets, Knorr-Held II-IV, between/within, lags, leave-one-province-out, alternative likelihoods, graph and prior, measured obesity, restricted spatial regression, influence, 2018 step); S20-S29, S31 (k-nearest-neighbour graphs, leave-one-year-out, period restrictions, lag-matched samples, Spatial+, mixing-parameter priors, residual Moran's I by specification, survey-covariate checks, forced-control removal, time-trend and holdout checks, untransformed covariates); additional analyses A-F (`83_extra_fits.R`); +1 sensitivity (`84_refit_green_plus1.R`). Phase C: S7 summary, S17 regional slopes, S44 checks of associations opposite to the expected direction |
| 11 | `extra/02_G_numeric.py`, `03_summarize.py`, `04_tables_md.py`, `05_posthoc.py`; `results_plus1_check/compare.py` | Summaries and post hoc tables (no new fits) |
| 11b | `extra2/e2_00_prep.py`, `e2_05_make_sets.py`, `e2_10_scan_thresholds.py`, `e2_11_village_rule_thresholds.py`, `e2_12_prep_a205.py`; `extra2/make_jobs.py <run root> <boundary file> <graph> <jobs folder>`, then `RSLOT=01_스크립트/rslot.sh extra2/runner.sh <log folder> E2.jobs`; `e2_20_split_perm.py`, `e2_30_maps_figs.py`, `e2_40_combine.py` | Round-2 additional analyses (post hoc), written to `04_분석/OBES/results_extra2` |
| 12 | `04_분석/OBES/MS/build/numbers_obes.py`, `01_스크립트/99O4_obes_maps_print.py`, `figs_obes_print.py`, `figs_si_print.py` | Reported numbers of the main chain, Figure 1 district maps, Figure 2, Supplementary Figures S1-S4 |

All fits use `ST_THREADS=1` (`num.threads = 1:1`); at most four R jobs run at once.

Names in the R code: prior lists `PRI_BYM2` (PC priors of the BYM2 district term), `PRI_IID` (PC precision prior of the iid, random-walk, AR1, Besag and generic terms) and `PRI_FIX` (fixed effects); loop variables at the top level of a script have two or more letters.

## Software

R 4.6.0 with R-INLA 25.10.19, sf 1.1.1 and spdep 1.4.2; Python 3.9.6 with pandas, numpy, scipy, geopandas and matplotlib. Newer R-INLA versions may shift DIC and WAIC by a few points without changing any OR.

## Versions

- v1.2.1: code of the final manuscript build: figure scripts no longer need the manuscript text module (`figs_common.py`), Supplementary Figure script `figs_si_print.py`, round-2 additional analyses (`extra2/`), updated document helpers and reference metadata (`refs_obes.json` identical to the build copy); R prior lists renamed (`PRI_BYM2`, `PRI_IID`, `PRI_FIX`) and top-level loop variables lengthened (results unchanged); boundary-file default `data/boundary/final.shp`; data-gate work-folder and column names without release tags (`input_first/`, `vil_n_corrected`, `T07b_마을체육시설_계_정정_전칸.csv`). Run from this copy, `numbers_obes.py` and all figure scripts reproduced the manuscript build byte for byte, the data gate rebuilt the corrected panel with identical values, and the main model (S0, S2, S3) reproduced the main odds ratios.
- v1.2.0: new manuscript title; pre-model data gate and panel correction code (40 panel cells of village sports facilities corrected for source-table row shifts and double counting; basic livelihood security rates replaced with official values restored from recipient counts and population); full chain refitted on the corrected panel with a job runner; additional sensitivity analyses (S20-S29, S31, analyses A-F, post hoc tables); `make_green_variants.py` now writes the `green_p1_raw` column needed by fit C_plus1_drop16; print-size figure scripts renamed (`99O4_obes_maps_print.py`, `figs_obes_print.py`); scripts that were replaced are marked SUPERSEDED.
- v1.1.0: print-size figure scripts (Figure 1 maps with scale bar, Figure 2 with an estimate column), +1 offset sensitivity code for the green-zone covariate, updated reference metadata, document helper renamed to `docx_engine.py`.
- v1.0.0: first release.

## License

MIT (see `LICENSE`).
