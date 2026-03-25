#!/bin/bash

local_dir=`cat ../util/local_dir.txt`
working_dir="."
program="${local_dir}/`cat ../util/program.txt`"
database="${local_dir}/`cat ../util/database.txt`"

address="${working_dir}/run" #_$(date +"%s")"

declare -a graphs=("fbA" "fbC" 
                  "rmat_s20_ef48_0.25_0.25_0.25_0.25_exp"
                  "rmat_s20_ef48_0.25_0.25_0.25_0.25_uni"
                  "rmat_s20_ef48_0.45_0.15_0.15_0.25_exp"
                  "rmat_s20_ef48_0.45_0.15_0.15_0.25_uni"
                  "rmat_s20_ef48_0.55_0.15_0.15_0.15_exp"
                  "rmat_s20_ef48_0.55_0.15_0.15_0.15_uni"
                )


# fuction to run partitioning commands
get_partitioning_file_name() {
  local seed=$1
  local cores=$2
  local p=$3
  local file=""

  case $p in
    metis)
      file="${database}/part/${graph}.graph.part.$cores.metis.s$seed"
      ;;
    umetis)
      file="${database}/part/${graph}_u.graph.part.$cores.metis.s$seed"
      ;;
    uwmetis)
      file="${database}/part/${graph}_uw.graph.part.$cores.metis.s$seed"
      ;;
    ukahip)
      file="${database}/part/${graph}_u.graph.part.${cores}.kahip.s${seed}"
      ;;
    none | random | diststream)
      return 0  # Always true
      ;;
    *)
      echo "Unknown partitioning method: $p"
      return 1
      ;;
  esac

  # Check if the file exists
  if [[ -f "$file" ]]; then
    return 0  # True: file exists
  else
    return 1  # False: file does not exist
  fi
}



for graph in "${graphs[@]}"
do
  echo "${graph}"

  for seed in 0 1 2 #3 4 5 6 7 8 9
  do
    for cores in 1 4 8 16 32
    do
      for p in uwmetis metis umetis ukahip diststream random
      do
        for k in 4 8 16 32
        do
          for algo in kMM kCS MRepS MRepLM
          do
            if [ "$cores" -eq 1 ] && [ "$p" != random ]; then
              continue
            fi
            if ! get_partitioning_file_name $seed $cores $p; then
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
