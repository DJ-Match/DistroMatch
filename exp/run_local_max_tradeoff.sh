#!/bin/bash

pwd > exp/util/local_dir.txt

cd exp/local_max_tradeoff
rm -r run
 
# echo "--- Generating experiment run files ---"
/bin/bash gen_exp_local_max_tradeoff.sh
 
echo "--- Running experiment ---"
/bin/bash ../util/run.sh run 1 # Change this to run experiments parallel if there are enough cores

echo "--- Extracting and Plotting results ---"
/bin/bash ../util/extract_result.sh run
mv run/results.txt results_run_local_max_tradeoff.txt
 
# rm -r run_local_max_tradeoff
 

python3 plot_local_max_tradeoff.py
