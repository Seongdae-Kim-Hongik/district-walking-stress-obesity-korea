#!/bin/bash
# Runs the model chain.
S=${ST_ROOT:-.}/01_스크립트; cd $S; T=OBES; L=${ST_ROOT:-.}/04_분석/$T
export ST_GRAPH=$L/graph_593.graph ST_THREADS=1 ST_VER=std3
ST_FH_ONLY=S2,S3 ./rslot.sh Rscript 81fh_obes.R $T > $L/run_fh_main.log 2>&1
ST_FH_ONLY=S7 ST_MV_STOP=1 Rscript 81fh_obes.R $T >> $L/run_fh_main.log 2>&1
ST_FH_ONLY=S4,S5,S6,S8,S9 ./rslot.sh Rscript 81fh_obes.R $T > $L/run_fh_sens.log 2>&1 &
ST_FH_ONLY=S10,S11,S12,S13,S14,S16,S18 ./rslot.sh Rscript 81fh_obes.R $T > $L/run_fh_sens2.log 2>&1 &
for R in 1-50 51-100 101-150 151-200; do ST_FH_ONLY=S7 ST_MV_RANGE=$R ./rslot.sh Rscript 81fh_obes.R $T > $L/run_fh_s7_$R.log 2>&1 & done
wait
ST_FH_ONLY=S7,S17,S44 ./rslot.sh Rscript 81fh_obes.R $T >> $L/run_fh_main.log 2>&1
echo "chain done $(date)" >> $L/run_fh_main.log; echo ALLDONE
