#!/bin/bash

OUT_DIR=/PATH/TO/OUT/DIR #set output directory
K=1000 #set carrying capacity

cd $OUT_DIR
for SEED in {1..3} #set number of replicates
do
slim -d seed=$SEED -d K=$K GEMnWF.slim &> GEMnWF.seed${SEED}.oe &
done