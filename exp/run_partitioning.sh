pwd > exp/util/local_dir.txt

cd exp/partitioning
rm -r run
 
echo "--- Generating experiment run files ---"
/bin/bash gen_exp
_partitioning.sh
 
echo "--- Running experiment ---"
/bin/bash ../util/run.sh run 1 # Change this to run experiments parallel if there are enough cores

echo "--- Extracting and Plotting results ---"
/bin/bash ../util/crawl_partitioning_times.sh

/bin/bash ../util/extract_result.sh run
mv run/results.txt results_run_partitioning.txt
cp results_run_partitioning.txt results_run_scaling.txt

#python plot_partitioning.py
#python plot_efficiency.py
