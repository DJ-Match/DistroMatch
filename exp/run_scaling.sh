#!/bin/bash

pwd > exp/util/local_dir.txt

cd exp/scaling
rm -r run_scaling
 
echo "--- Generating experiment run files ---"
/bin/bash gen_exp_scaling.sh
 
echo "--- Running experiment ---"
/bin/bash ../util/run.sh run_scaling 1 # Change this to run experiments parallel if there are enough cores

echo "--- Extracting and Plotting results ---"
/bin/bash ../util/extract_result.sh run_scaling
mv run_scaling/results.txt results_run_scaling.txt

python3 plot_efficiency.py
