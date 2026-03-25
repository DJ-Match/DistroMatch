#!/bin/bash

pwd > exp/util/local_dir.txt

cd exp/rmat
rm -r run
 
# echo "--- Generating experiment run files ---"
/bin/bash gen_exp_rmat.sh
 
echo "--- Running experiment ---"
/bin/bash ../util/run.sh run 1 # Change this to run experiments parallel if there are enough cores

echo "--- Extracting and Plotting results ---"
/bin/bash ../util/extract_run_comp.sh run
mv run/results.csv results_run_rmat.csv

python3 plot_rmat.py
