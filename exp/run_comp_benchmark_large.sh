#!/bin/bash

pwd > exp/util/local_dir.txt

cd exp/comp_benchmark
rm -r run_large
 
echo "--- Generating experiment run files ---"
/bin/bash gen_exp_comp_large.sh
 
echo "--- Running experiment (large) ---"
/bin/bash ../util/run.sh run_large 1 # Change this to run experiments parallel if there are enough cores


echo "--- Extracting and Plotting results ---"
/bin/bash ../util/extract_run_comp.sh run_large
mv run_large/results.csv results_run_seq_comp_large.csv

python3 plot_comp.py --size large
