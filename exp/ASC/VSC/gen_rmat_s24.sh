#!/bin/bash


local_dir=`$HOME/exp/gkm_scaling`
working_dir="."
program="$HOME/DistroMatch/build/Release/DistroMatch"
database="$HOME/data"

address="${working_dir}/run_s24_p$1" #_only_scaling" #_$(date +"%s")"


mkdir $address

declare -a graphs=(
                  "rmat_s24_ef24_0.25_0.25_0.25_0.25_exp"
                  "rmat_s24_ef24_0.25_0.25_0.25_0.25_uni"
                  "rmat_s24_ef24_0.45_0.15_0.15_0.25_exp"
                  "rmat_s24_ef24_0.45_0.15_0.15_0.25_uni"
                  "rmat_s24_ef24_0.55_0.15_0.15_0.15_exp"
                  "rmat_s24_ef24_0.55_0.15_0.15_0.15_uni"
                )

declare -A partitioning_methods
partitioning_methods[GKM]="random"

for graph in "${graphs[@]}"
do
  echo "${graph}"
  for seed in 0 1 # 2 3 4 5 6 7 8 9
  do
    for cores in $1 # 2 4 8 16 32
    do
      for k in  8 # 16 32
      do
        for algo in GKM # CPal kM MRep
        do
          for p in ${partitioning_methods[$algo]}
          do
            if [ "$cores" -eq 1 ] && [ "$p" != random ]; then
              continue
            fi
            address2=${address}/${graph}/${algo}/${p}/${k}/${cores}/${seed}
            parameters=''

            mkdir -p ${address2}
            echo "now=\"\$(date +\"%T\")\"; echo [\$now] $address2" >> ${address2}/gen_map.sh

            ext='mtx'

            echo "mpirun -np ${cores}  ${program} -i ${database}/${graph}.${ext} --algo=${algo} --k=${k} --seed=${seed} --partitioning=$p  > ${address2}/result.txt" >> ${address2}/gen_map.sh

          done
        done
      done
    done
  done
done

chmod +x -R ${address}
