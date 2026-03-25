#!/bin/bash

pwd > exp/util/local_dir.txt

cd exp/comp_benchmark
rm -r run_small
 
echo "--- Generating experiment run files ---"
/bin/bash gen_exp_comp.sh
 
echo "--- Running experiment (small) ---"
/bin/bash ../util/run.sh run_small 1 # Change this to run experiments parallel if there are enough cores

echo "--- Extracting and Plotting results ---"
/bin/bash ../util/extract_run_comp.sh run_small
mv run_small/results.csv results_run_seq_comp.csv

python3 plot_comp.py --size small
