# LTN_analysis

Please set the working directory before running `run.sh`.
```bash
export WORK_DIR=.../LTN_analysis/
$WORK_DIR/run.sh 100 1000
```  
Results from LTN for the numerical examples of covariance estimation are stored in the folder `results/covariance_estimation/risk`. The folder `results/cross_group_comparison/single_otu/summary` and `results/cross_group_comparison/multi_otu/summary` contain results in the cross-group comparison examples, including the PJAPs and labels for making the ROC curves as well as the ROC plots. We fitted COAT and DirFactor to simulated datasets with code from https://github.com/yuanpeicao/COAT.git and https://github.com/boyuren158/DirFactor-fix respectively. 
The plots and tables shown in the application section are stored in `results/cross_group_comparison/application/`. 



