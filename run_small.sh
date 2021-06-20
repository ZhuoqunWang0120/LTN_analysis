#!/bin/bash

nsim1=$1
nsim2=$2
mkdir $WORK_DIR/cache
mkdir -p $WORK_DIR/results/covariance_estimation/dtm
mkdir -p $WORK_DIR/results/covariance_estimation/ln

# make data
echo "processing data"
$WORK_DIR/src/data/process_data_otu50.R $WORK_DIR
$WORK_DIR/src/data/process_data_simulation.R $WORK_DIR
$WORK_DIR/src/data/process_data_application.R $WORK_DIR
$WORK_DIR/src/data/dtm_data.R $WORK_DIR
$WORK_DIR/src/experiments/covariance_estimation/dtm_clrcov.R $WORK_DIR 5

# covariance estimation
echo "generating datasets for covariance estimation"
for i in `seq 1 $nsim1`; do
    $WORK_DIR/src/experiments/covariance_estimation/dtm_sim.R $WORK_DIR $i
    for j in `seq 1 3`;do
        $WORK_DIR/src/experiments/covariance_estimation/ln_sim.R $WORK_DIR $i $j
    done
done

declare -a lamvec=(0 10)
for lambda in "${lamvec[@]}";do
for i in `seq 1 $nsim1`; do
    echo "fitting LTN on covariance simulation "$i
    $WORK_DIR/src/experiments/covariance_estimation/dtm_fit_parse.R --i $i --niter 5 --nmc 5 --lambda $lambda --WORK_DIR $WORK_DIR
for j in `seq 1 3`;do
    $WORK_DIR/src/experiments/covariance_estimation/ln_fit_parse.R --SEED $i --modelCov $j --niter 5 --nmc 5 --lambda $lambda --WORK_DIR $WORK_DIR
done
done
$WORK_DIR/src/experiments/covariance_estimation/collect_results.R $WORK_DIR $lambda $nsim1
done

# cross-group comparison
echo "generating datasets for cross-group comparison"
for i in `seq 1 $nsim2`; do
for h in `seq 0 1`;do
    $WORK_DIR/src/experiments/cross_group_comparison/single_otu_sim.R $WORK_DIR $i $h
    $WORK_DIR/src/experiments/cross_group_comparison/multi_otu_sim.R $WORK_DIR $i $h
done
done

declare -a lamvec=(0 10)
for i in `seq 1 $nsim2`; do
echo "fitting LTN on cross-group comparison simulation "$i
for h in `seq 0 1`;do
$WORK_DIR/src/experiments/cross_group_comparison/fit.R --h $h --niter 5 --reff --reffcov 1 --pi_only --lambda 0 --i $i --s single_otu --WORK_DIR $WORK_DIR
$WORK_DIR/src/experiments/cross_group_comparison/fit.R --h $h --niter 5 --reff --reffcov 1 --pi_only --lambda 0 --i $i --s multi_otu --WORK_DIR $WORK_DIR
for lambda in "${lamvec[@]}";do
$WORK_DIR/src/experiments/cross_group_comparison/fit.R --h $h --niter 5 --reff --reffcov 2 --pi_only --lambda $lambda --i $i --s single_otu --WORK_DIR $WORK_DIR
$WORK_DIR/src/experiments/cross_group_comparison/fit.R --h $h --niter 5 --reff --reffcov 2 --pi_only --lambda $lambda --i $i --s multi_otu --WORK_DIR $WORK_DIR
done
done
done

for lambda in "${lamvec[@]}";do
$WORK_DIR/src/experiments/cross_group_comparison/collect_results.R $WORK_DIR $lambda
done

# application
declare -a lamvec=(0.1 10)
for lambda in "${lamvec[@]}";do
for i in `seq 1 9`; do
$WORK_DIR/src/experiments/application/application.R $WORK_DIR 5 $i 2 $lambda
done
$WORK_DIR/src/experiments/application/collect_results.R $WORK_DIR $lambda
done



