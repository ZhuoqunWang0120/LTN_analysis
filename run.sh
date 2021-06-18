#!/bin/bash

cov_nsim=$1
mkdir $WORK_DIR/LTN_analysis/results
mkdir $WORK_DIR/LTN_analysis/cache
$WORK_DIR/LTN_analysis/src/data/process_data_otu50.R $WORK_DIR
$WORK_DIR/LTN_analysis/src/data/process_data_simulation.R $WORK_DIR
$WORK_DIR/LTN_analysis/src/data/process_data_application.R $WORK_DIR
$WORK_DIR/LTN_analysis/src/data/dtm_data.R $WORK_DIR
$WORK_DIR/LTN_analysis/src/experiments/covariance_estimation/dtm_clrcov.R $WORK_DIR 1000000
for i in `seq 1 $cov_nsim`; do
    $WORK_DIR/LTN_analysis/src/experiments/covariance_estimation/dtm_sim.R $WORK_DIR $i
done