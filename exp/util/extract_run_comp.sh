#!/bin/bash

declare -A algora_des=( ["NC"]="node_centered-SUM\s*\+\s*threshold\s*0\.2"
                        ["NCb"]="node_centered-B_SUM\s*\+\s*threshold\s*0\.2"
                        ["kEC"]="k-Edge Coloring \(CC, RL\)")

# Loop through each address
for address in $1; do
  # Set the output CSV file name
  output_file="$address/results.csv"
  # Write the CSV header
  echo "graph,algo,k,seed,n,m,read_time,running_time,matching_size,matching_weight" > "$output_file"
  # Find each result.txt file in the given address
  for result_file in $(find "$address" -iname result.txt); do
    # Parse the path to get graph, algo, k, and seed
    path_components=$(echo "$result_file" | sed 's|/| |g')
    read -r _ graph algo k seed _ <<<"$path_components"

    # Initialize variables
    n="" m="" read_time="" running_time="" matching_size="" matching_weight=""

    # Parse based on the algorithm type
    if [[ $algo == stk ]]; then
      # Extract relevant data for 'stk' algorithm
      line=$(tail -n 2 "$result_file" | head -n 1)
      read -r n m matching_weight matching_size _ read_time proc_time read_proc_time post_proc_time _ <<<"$line"
      running_time=$(echo "$read_proc_time $post_proc_time" | awk '{printf "%f", $1 + $2}')

    elif [[ $algo == stk-dp ]]; then
      # Extract relevant data for 'stk' algorithm
      line=$(tail -n $((k+3)) "$result_file" | head -n 1)
      #echo "first:" $((k+3)) $line
      read -r _ _ _ _ _ read_time proc_time read_proc_time _ _ <<<"$line"
      line=$(tail -n 2 "$result_file" | head -n 1)
      #echo "second: " 2 $line
      read -r n m matching_weight matching_size _ _ _ post_proc_time _ dp_time <<<"$line"
      running_time=$(echo "$read_proc_time $post_proc_time $dp_time" | awk '{printf "%f", $1 + $2 + $3}')
      #echo $read_proc_time $post_proc_time $dp_time $running_time

    elif [[ $algo == CCM_* || $algo == MRep_* || $algo == kM_* ]]; then
     # echo $result_file 
     # Extract n, m, and matching_weight for 'CCM'
      n=$(grep -oP '^vertices: \K\d+' "$result_file")
      local_edges=$(grep -oP 'local_edges: \K\d+' "$result_file")
      cross_edges=$(grep -oP 'cross_edges: \K\d+' "$result_file")
      m=$((local_edges + cross_edges))

      read_time=$(grep -oP 'read time: \K[\d\.]+' "$result_file")
      running_time=$(grep -oP 'running_time: \K[\d\.]+' "$result_file")
      matching_size=$(grep -oP 'total_size: \K\d+' "$result_file")
      matching_weight=$(grep -oP 'total_weight: \K[\d\.]+' "$result_file")

    #echo "$graph,$algo,$k,$seed,$n,$m"

  elif [[ "$algo" == "NC" || $algo == NCb || $algo == kEC ]]; then
      # Extract n, m, and matching_weight for 'NC'
      n=$(grep -oP '^%n,m \K\d+' "$result_file")
      m=$(grep -oP '^%n,m \d+,\K\d+' "$result_file")
      name=${algora_des["$algo"]}
#      echo $name
      read_time=$(grep -oP 'Input I/O took \K[\d\.]+' "$result_file")
#      echo "grep -oP '\|\s*'"$name"'\s*\|\s*\K[\d,]+' "$result_file" | tr -d ',')"
      matching_weight=$(grep -oP '\|\s*'"$name"'\s*\|\s*\K[\d,]+' "$result_file" | tr -d ',')
      matching_size=0
      running_time=$(grep -oP '\|\s*'"$name"'\s*\|\s*[\d,]+\s*\|\s*\K[\d\.]+' "$result_file")
fi

    # Append extracted data to CSV file
    echo "$graph,$algo,$k,$seed,$n,$m,$read_time,$running_time,$matching_size,$matching_weight" >> "$output_file"
  done
done

