#!/bin/bash


local_dir=`cat ../util/local_dir.txt`
working_dir="."
program="${local_dir}/`cat ../util/program.txt`"
database="${local_dir}/`cat ../util/database.txt`"

address="${working_dir}/run_scaling" #_$(date +"%s")"


mkdir $address

declare -a graphs=("fbA" "fbC" 
                  "rmat_s20_ef48_0.25_0.25_0.25_0.25_exp"
                  "rmat_s20_ef48_0.25_0.25_0.25_0.25_uni"
                  "rmat_s20_ef48_0.45_0.15_0.15_0.25_exp"
                  "rmat_s20_ef48_0.45_0.15_0.15_0.25_uni"
                  "rmat_s20_ef48_0.55_0.15_0.15_0.15_exp"
                  "rmat_s20_ef48_0.55_0.15_0.15_0.15_uni"
                )
mpi="mpirun --mca btl tcp,self --mca btl_tcp_if_include lo"

for graph in "${graphs[@]}"
do
  echo "${graph}"
  for seed in 0 1 2 #3 4 5 6 7 8 9
  do
    for cores in 1 4 8 16 32
    do
      for k in 4 8 16 32
      do
        for algo in kMM kCS MRepS MRepLM
        do
          for p in random
          do
            if [ "$cores" -eq 1 ] && [ "$p" != random ]; then
              continue
            fi
            address2=${address}/${graph}/${algo}/${p}/${k}/${cores}/${seed}
            parameters=''

            mkdir -p ${address2}
            echo "now=\"\$(date +\"%T\")\"; echo [\$now] $address2" >> ${address2}/gen_map.sh

            ext='mtx'

            echo "$mpi -np ${cores}  ${program} -i ${database}/${graph}.${ext} --algo=${algo} --k=${k} --seed=${seed} --partitioning=$p  > ${address2}/result.txt" >> ${address2}/gen_map.sh

          done
        done
      done
    done
  done
done

chmod +x -R ${address}
