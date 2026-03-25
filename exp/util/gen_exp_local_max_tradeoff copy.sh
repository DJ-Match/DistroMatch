#!/bin/bash

home=`cat home.txt`
working_dir="${home}/`cat rel_path.txt`"
program="${home}/`cat program.txt`"
database="${home}/data/ready"
parameters=""

address="${working_dir}/run_local_max_tradeoff" #_$(date +"%s")"

for graph in ${database}/*.mtx
do
  graph="${graph##*/}"
  graph="${graph%.*}"
  if [[ "${graph}" == *_u ]]; then
    continue
  fi
  #if [[ "${graph}" != rmat_s20*uni* ]]; then
  #  continue
  #fi


  echo ${graph}
  for seed in 0 1 2
  do
    for k in 2 8 16
    do
      for cores in 1 8 32
      do
        for epsilon in 0 0.0001 0.001 0.01 0.1 0.5 1
        do
          for algo in kM
          do
            address2=${address}/${graph}/${algo}/${epsilon}/${k}/${cores}/${seed}

            mkdir -p ${address2}

            echo "now=\"\$(date +\"%T\")\"; echo [\$now] $address2" >> ${address2}/gen_map.sh
            echo "mpirun -np ${cores}  ${program} -i ${database}/${graph}.mtx --algo=${algo} --k=${k} --CPal_epsilon=${epsilon} --seed=$seed --partitioning=random > ${address2}/result.txt" >> ${address2}/gen_map.sh

          done
        done
      done
    done
  done
done

chmod +x -R ${address}
