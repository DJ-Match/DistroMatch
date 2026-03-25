#!/bin/bash

local_dir=`cat ../util/local_dir.txt`
working_dir="."
program="${local_dir}/`cat ../util/program.txt`"
database="${local_dir}/`cat ../util/database.txt`"

nc_program="${local_dir}/../DJ-Match/build/Release/DJMatch"
stk_program="${local_dir}/../GStream-dev/build/apps/kstmatch"

declare -a graphs=("GAP-urand" "GAP-kron" "MOLIERE_2016" 
                  "com-Friendster_uni" "AGATHA_2015_uni" "mycielskian20_uni"
                  "com-Friendster_exp" "AGATHA_2015_exp" "mycielskian20_exp" )


output_base="${working_dir}/run_large"
mpi="mpirun --mca btl tcp,self --mca btl_tcp_if_include lo"

echo $database
algos=("stk" "kCS_128")


for graph in "${graphs[@]}"
do
  echo "${graph}"

  # Loop over each seed and value of k
  for seed in 0 1 2 ; do
    for k in 2 8 32 64; do
      # Define output directory
      address="${output_base}/${graph}"

      # Loop over each algorithm and create the appropriate command
      for algo in "${algos[@]}"; do
        # Construct the address and ensure directories exist
        algo_dir="${address}/${algo}/${k}/${seed}"
        mkdir -p "${algo_dir}"

        # Start generating the script for this experiment
        gen_script="${algo_dir}/gen_map.sh"
        echo "#!/bin/bash" > "${gen_script}"
        echo "now=\"\$(date +\"%T\")\"; echo [\$now] Running ${algo_dir}" >> "${gen_script}"

        # Initialize the specific command based on the algorithm
        case "$algo" in
          "kCS_128")
            echo "$mpi -np 128 ${program} -i ${database}/${graph}.mtx --algo=kCS --k=${k} --seed=${seed} --partitioning=random > ${algo_dir}/result.txt" >> "${gen_script}"
            ;;
          "stk")
            echo "${stk_program} -f ${database}/${graph}.mtx -a stk -k ${k} --sim > ${algo_dir}/result.txt" >> "${gen_script}"
            ;;
        esac
      done
    done
  done
done

chmod +x -R "${output_base}"

