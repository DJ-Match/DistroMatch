#!/bin/bash

local_dir=`cat ../util/local_dir.txt`
working_dir="."
program="${local_dir}/`cat ../util/program.txt`"
database="${local_dir}/`cat ../util/database.txt`"

address="${working_dir}/run" #_$(date +"%s")"

declare -a graphs=("fbA" "fbB" "fbC" 
                  "rmat_s20_ef48_0.25_0.25_0.25_0.25_exp"
                  "rmat_s20_ef48_0.25_0.25_0.25_0.25_uni"
                  "rmat_s20_ef48_0.45_0.15_0.15_0.25_exp"
                  "rmat_s20_ef48_0.45_0.15_0.15_0.25_uni"
                  "rmat_s20_ef48_0.55_0.15_0.15_0.15_exp"
                  "rmat_s20_ef48_0.55_0.15_0.15_0.15_uni"
                )

mpi="mpirun --mca btl tcp,self --mca btl_tcp_if_include lo"

# for  in ${database}/*.mtx
for graph in "${graphs[@]}"
do
  echo ${graph}
  for seed in 0 1 2
  do
    for k in 2 8 16
    do
      for cores in 1 8 32
      do
        for epsilon in 0 0.0001 0.001 0.01 0.1 0.5 1
        do
          for algo in kMM
          do
            address2=${address}/${graph}/${algo}/${epsilon}/${k}/${cores}/${seed}

            mkdir -p ${address2}

            echo "now=\"\$(date +\"%T\")\"; echo [\$now] Running ${address2}" >> ${address2}/gen_map.sh
            echo "$mpi -np ${cores}  ${program} -i ${database}/${graph}.mtx --algo=${algo} --k=${k} --CPal_epsilon=${epsilon} --seed=$seed --partitioning=random > ${address2}/result.txt" >> ${address2}/gen_map.sh

          done
        done
      done
    done
  done
done


chmod +x -R ${address}
