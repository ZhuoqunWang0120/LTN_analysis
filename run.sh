#!/bin/bash

mkdir $WORK_DIR/LTN_analysis/results
mkdir $WORK_DIR/LTN_analysis/cache
$WORK_DIR/LTN_analysis/src/data/process_data_otu50.R $WORK_DIR
$WORK_DIR/LTN_analysis/src/data/process_data_simulation.R $WORK_DIR
$WORK_DIR/LTN_analysis/src/data/process_data_application.R $WORK_DIR