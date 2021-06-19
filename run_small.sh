#!/bin/bash

cov_nsim=$1
mkdir $WORK_DIR/LTN_analysis/results
mkdir $WORK_DIR/LTN_analysis/cache
mkdir $WORK_DIR/LTN_analysis/results/covariance_estimation
mkdir $WORK_DIR/LTN_analysis/results/covariance_estimation/dtm

$WORK_DIR/LTN_analysis/src/data/process_data_otu50.R $WORK_DIR
$WORK_DIR/LTN_analysis/src/data/process_data_simulation.R $WORK_DIR
$WORK_DIR/LTN_analysis/src/data/process_data_application.R $WORK_DIR
$WORK_DIR/LTN_analysis/src/data/dtm_data.R $WORK_DIR
$WORK_DIR/LTN_analysis/src/experiments/covariance_estimation/dtm_clrcov.R $WORK_DIR 10

declare -a lamvec=(0 0.1 1 10)
for lambda in "${lamvec[@]}";do
for i in `seq 1 $cov_nsim`; do
    $WORK_DIR/LTN_analysis/src/experiments/covariance_estimation/dtm_sim.R $WORK_DIR $i
    $WORK_DIR/LTN_analysis/src/experiments/covariance_estimation/dtm_fit_parse.R --i $i --niter 3 --nmc 10 --lambda $lambda --WORK_DIR $WORK_DIR
done
done
