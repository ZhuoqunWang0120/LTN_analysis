#!/bin/bash

cov_nsim=$1
mkdir $WORK_DIR/results
mkdir $WORK_DIR/cache
mkdir -p $WORK_DIR/results/covariance_estimation/dtm
mkdir -p $WORK_DIR/results/covariance_estimation/ln

# make data
$WORK_DIR/src/data/process_data_otu50.R $WORK_DIR
$WORK_DIR/src/data/process_data_simulation.R $WORK_DIR
$WORK_DIR/src/data/process_data_application.R $WORK_DIR
$WORK_DIR/src/data/dtm_data.R $WORK_DIR
$WORK_DIR/src/experiments/covariance_estimation/dtm_clrcov.R $WORK_DIR 1000000

# covariance estimation
for i in `seq 1 $cov_nsim`; do
    $WORK_DIR/src/experiments/covariance_estimation/dtm_sim.R $WORK_DIR $i
    for j in `seq 1 3`;do
        $WORK_DIR/src/experiments/covariance_estimation/ln_sim.R $WORK_DIR $i $j
    done
done
declare -a lamvec=(0 0.1 1 10)
for lambda in "${lamvec[@]}";do
for i in `seq 1 $cov_nsim`; do
    $WORK_DIR/src/experiments/covariance_estimation/dtm_fit_parse.R --i $i --niter 10000 --nmc 1000000 --lambda $lambda --WORK_DIR $WORK_DIR
done
done

