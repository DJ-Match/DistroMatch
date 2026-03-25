#!/bin/bash

local_dir=`cat ../util/local_dir.txt`
working_dir="."
program="${local_dir}/`cat ../util/program.txt`"
database="${local_dir}/`cat ../util/database.txt`"

nc_program="${local_dir}/../DJ-Match/build/Release/DJMatch"
stk_program="${local_dir}/../GStream-dev/build/apps/kstmatch"


output_base="${working_dir}/run_small"

declare -a graphs=("fbA" "fbB" "fbC" 
                  "rmat_s20_ef48_0.25_0.25_0.25_0.25_exp"
                  "rmat_s20_ef48_0.25_0.25_0.25_0.25_uni"
                  "rmat_s20_ef48_0.45_0.15_0.15_0.25_exp"
                  "rmat_s20_ef48_0.45_0.15_0.15_0.25_uni"
                  "rmat_s20_ef48_0.55_0.15_0.15_0.15_exp"
                  "rmat_s20_ef48_0.55_0.15_0.15_0.15_uni"
                )

# Algorithms to test and their specific commands
algos=("kCS_1" "kCS_32" "kMM_1" "kMM_32" "MRepS_1" "MRepS_32" "MRepLM_1" "MRepLM_32" "kEC" "stk-dp" "stk")

mpi="mpirun --mca btl tcp,self --mca btl_tcp_if_include lo"

mkdir -p $output_base

# Loop over each graph file in the database directory
for graph in "${graphs[@]}"
do
  echo "${graph}"

  # Loop over each seed and value of k
  for seed in 0 1 2; do
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
          "kCS_1")
            echo "$mpi -np 1 ${program} -i ${database}/${graph}.mtx --algo=kCS --k=${k} --seed=${seed} --partitioning=random > ${algo_dir}/result.txt" >> "${gen_script}"
            ;;
          "kCS_32")
            echo "$mpi -np 32 ${program} -i ${database}/${graph}.mtx --algo=kCS --k=${k} --seed=${seed} --partitioning=random > ${algo_dir}/result.txt" >> "${gen_script}"
            ;;
          "kMM_1")
            echo "$mpi -np 1 ${program} -i ${database}/${graph}.mtx --algo=kMM --k=${k} --seed=${seed} --partitioning=random > ${algo_dir}/result.txt" >> "${gen_script}"
            ;;
          "kMM_32")
            echo "$mpi -np 32 ${program} -i ${database}/${graph}.mtx --algo=kMM --k=${k} --seed=${seed} --partitioning=random > ${algo_dir}/result.txt" >> "${gen_script}"
            ;;
          "MRepS_1")
            echo "$mpi -np 1 ${program} -i ${database}/${graph}.mtx --algo=MRepS --k=${k} --seed=${seed} --partitioning=random > ${algo_dir}/result.txt" >> "${gen_script}"
            ;;
          "MRepS_32")
            echo "$mpi -np 32 ${program} -i ${database}/${graph}.mtx --algo=MRepS --k=${k} --seed=${seed} --partitioning=random > ${algo_dir}/result.txt" >> "${gen_script}"
            ;;
          "MRepLM_1")
            echo "$mpi -np 1 ${program} -i ${database}/${graph}.mtx --algo=MRepLM --k=${k} --seed=${seed} --partitioning=random > ${algo_dir}/result.txt" >> "${gen_script}"
            ;;
          "MRepLM_32")
            echo "$mpi -np 32 ${program} -i ${database}/${graph}.mtx --algo=MRepLM --k=${k} --seed=${seed} --partitioning=random > ${algo_dir}/result.txt" >> "${gen_script}"
            ;;
          "kEC")
            echo "${nc_program} --b=${k} -a k-ec ${database}/${graph}.edge_stream > ${algo_dir}/result.txt" >> "${gen_script}"
            ;;
          "stk")
              echo "${stk_program} -f ${database}/${graph}.mtx -a stk -k ${k} --sim > ${algo_dir}/result.txt" >> "${gen_script}"
              ;;
          "stk-dp")
              echo "${stk_program} -f ${database}/${graph}.mtx -a stk -k ${k} --dp --sim > ${algo_dir}/result.txt" >> "${gen_script}"
              ;;
        esac
      done
    done
  done
done

chmod +x -R "${output_base}"

