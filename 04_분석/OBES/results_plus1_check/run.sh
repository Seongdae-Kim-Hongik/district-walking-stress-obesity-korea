#!/bin/bash
# Runs one set of +1 sensitivity fits (comma-separated names, e.g. orig_rerun,A_halfmin). X_ROOT = project root, ST_SHP = boundary file, ST_GRAPH = district contact graph.
D=$(cd "$(dirname "$0")" && pwd); X_ROOT=${X_ROOT:-$(cd "$D/../../.." && pwd)}
cd "$D" && X_ROOT="$X_ROOT" ST_SHP="${ST_SHP:-final.shp}" ST_GRAPH="${ST_GRAPH:-$X_ROOT/04_분석/OBES/graph_593.graph}" ST_THREADS=${ST_THREADS:-1} PC_ONLY=$1 Rscript refit_green_plus1.R OBES
