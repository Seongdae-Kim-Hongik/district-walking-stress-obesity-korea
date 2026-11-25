#!/bin/bash
# Runs one set of +1 sensitivity fits (comma-separated names, e.g. orig_rerun,A_halfmin). X_ROOT = project root, ST_SHP = boundary file, ST_GRAPH = district contact graph.
# Run make_green_variants.py first; results go to this folder.
D=$(cd "$(dirname "$0")" && pwd); X_ROOT=${X_ROOT:-$(cd "$D/../../.." && pwd)}
cd "$X_ROOT/01_스크립트" && PC_RES="$D" ST_SHP="${ST_SHP:-data/boundary/final.shp}" ST_GRAPH="${ST_GRAPH:-$X_ROOT/04_분석/OBES/graph_593.graph}" ST_THREADS=${ST_THREADS:-1} PC_ONLY=$1 Rscript 84_refit_green_plus1.R OBES
