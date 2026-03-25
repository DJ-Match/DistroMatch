graphs=("fb_clusterA_full" "fb_clusterC_full" "rgg_s20_ef48_0.25_0.25_0.25_0.25_exp" "rgg_s20_ef48_0.45_0.15_0.15_0.25_uni" "rgg_s20_ef48_0.25_0.25_0.25_0.25_uni" "rgg_s20_ef48_0.55_0.15_0.15_0.15_exp" "rgg_s20_ef48_0.45_0.15_0.15_0.25_exp" "rgg_s20_ef48_0.55_0.15_0.15_0.15_uni")
# "rgg_s20_ef32_0.25_0.25_0.25_0.25_uni" "rgg_s20_ef48_0.55_0.15_0.15_0.15_exp")

cores_list=(4 8 16 32)  # Adjust the cores list based on the cores you want to use.
seeds=(0 1 2 ) # 3 4 5 6 7 8 9)       # Adjust the seeds as per your requirements.

# CSV Output file
output_file="graph_partition_times.csv"

# Write header to CSV file
echo "graph,cores,seed,metis,umetis,uwmetis,ukahip,random,diststream" > $output_file

# Function to extract timing from a file
extract_time() {
    local file_path=$1
    local pattern=$2
    local field=$3
    local result=$(grep "$pattern" "$file_path" | awk -v var="$field" '{print $var}')
    echo $result
}

# Loop over graphs, cores, and seeds
for graph in "${graphs[@]}"; do
    for cores in "${cores_list[@]}"; do
        for seed in "${seeds[@]}"; do
            # Paths to the different files for extraction
            metis_file="$HOME/data/ready/part/${graph}.graph.part.${cores}.metis.s${seed}.metis_out.txt"
            umetis_file="$HOME/data/ready/part/${graph}_u.graph.part.${cores}.metis.s${seed}.metis_out.txt"
            uwmetis_file="$HOME/data/ready/part/${graph}_uw.graph.part.${cores}.metis.s${seed}.metis_out.txt"
            ukahip_file="$HOME/data/ready/part/${graph}_u.graph.part.${cores}.kahip.s${seed}.kahip_out.txt"
            random_file="$HOME/exp/run_partitioning/${graph}/CCM/random/16/${cores}/${seed}/result.txt"
            stream_file="$HOME/exp/run_partitioning/${graph}/CCM/diststream/16/${cores}/${seed}/result.txt"
#          echo "grep 'Partitioning took' $stream_file | awk -v var=3 '{print \$var}'"

           # Extract the timings using the correct patterns
            metis_time=$(extract_time "$metis_file" 'Partitioning:' 2)
            umetis_time=$(extract_time "$umetis_file" 'Partitioning:' 2)
            uwmetis_time=$(extract_time "$uwmetis_file" 'Partitioning:' 2)
            ukahip_time=$(extract_time "$ukahip_file" 'time spent for partitioning' 5)
            random_time=$(extract_time "$random_file" 'Random generation took' 4)
            stream_time=$(extract_time "$stream_file" 'Partitioning took' 3)

#            echo HI $umetis_time
            # If any of the values is empty, set it to "N/A"
            metis_time=${metis_time:-"N/A"}
            umetis_time=${umetis_time:-"N/A"}
            uwmetis_time=${uwmetis_time:-"N/A"}
            ukahip_time=${ukahip_time:-"N/A"}
            random_time=${random_time:-"N/A"}
            stream_time=${stream_time:-"N/A"}

            # Append the data to the CSV file
            echo "${graph},${cores},${seed},${metis_time},${umetis_time},${uwmetis_time},${ukahip_time},${random_time},${stream_time}" >> $output_file
        done
    done
done

echo "Data collection complete! Output saved to $output_file"
