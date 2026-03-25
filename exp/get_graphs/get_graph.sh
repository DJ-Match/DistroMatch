#!/bin/bash

# Set to true to download Facebook traffic graphs for the experiments
get_fb_graphs=true
# Set to true to generate synthetic R-MAT graphs for the experiments
generate_rmat_graphs=true
# Set to true to generate only R-MAT graphs with 2^20 vertices; false to generate graphs of all sizes; works only with generate_rmat_graphs=true
generate_only_x20_rmat_graphs=false
# Set to true to format the generated graphs for use with offline reference algorithms
format_for_offline_reference_algos=true
# Set to true to download and prepare large real-world graphs for scalability testing
get_large_graphs=true

working_dir='exp/get_graphs'

mkdir 'data'

if $get_fb_graphs; then

echo '--- Downloading facebook graphs ---'
wget -O 'data/fbA.graph' 'https://nextcloud.inet.tu-berlin.de/public.php/dav/files/YLDnX4zCe4PcgGW'
wget -O 'data/fbB.graph' 'https://nextcloud.inet.tu-berlin.de/public.php/dav/files/JkBeGEnXCoM7pHB'
wget -O 'data/fbC.graph' 'https://nextcloud.inet.tu-berlin.de/public.php/dav/files/zNa5sjH3ECTn5ei'

fi

if $generate_rmat_graphs; then

echo '--- Generating rmat graphs ---'
program="python3 $working_dir/generate_rmat.py"

declare -a init_matrices=(
  "0.25 0.25 0.25 0.25"
  "0.55 0.15 0.15 0.15"
  "0.45 0.15 0.15 0.25"
)

declare -a scale_edgefactor=(
  "14 48"
  "15 48"
  "16 48"
  "17 48"
  "18 48"
  "19 48"
  "20 48"
  "21 24"
  "22 24"
  "23 24"
  "24 24"
)

if $generate_only_x20_rmat_graphs; then
  scale_edgefactor=(
    "20 48"
  )
fi


for s_ef in "${scale_edgefactor[@]}"; do
  read -r s ef <<< "$s_ef"
  for edge_weights in uni exp; do
    for init_m in "${init_matrices[@]}"; do
      read -r p_a p_b p_c p_d <<< "$init_m"
      
      echo "scale: $s, edgefactor: $ef, p_a: $p_a, p_b: $p_b, p_c: $p_c, p_d: $p_d, distribution: $edge_weights"
      
      $program --scale "$s" --edgefactor "$ef" \
        --p_a "$p_a" --p_b "$p_b" --p_c "$p_c" --p_d "$p_d" \
        --distribution "$edge_weights"

      mv "rmat_s${s}_ef${ef}_${p_a}_${p_b}_${p_c}_${p_d}_${edge_weights}.graph" "data/rmat_s${s}_ef${ef}_${p_a}_${p_b}_${p_c}_${p_d}_${edge_weights}.graph"


    done
  done
done

fi

echo '--- Converting into mtx format ---'
python3 exp/get_graphs/convert_metis_to_mtx.py 'data'

if $format_for_offline_reference_algos; then

echo '--- Converting into edge stream format... (for benchmark algorithms) ---'
python3 exp/get_graphs/convert_mtx_to_edge_stream.py 'data'

fi

if $get_large_graphs; then

echo '--- Downloading large graphs ---'
link="https://suitesparse-collection-website.herokuapp.com/MM/"

data=("SNAP/com-Friendster.tar.gz" "Sybrandt/MOLIERE_2016.tar.gz" "Sybrandt/AGATHA_2015.tar.gz" "GAP/GAP-urand.tar.gz" "GAP/GAP-kron.tar.gz" "Mycielski/mycielskian20.tar.gz")

for dataset in "${data[@]}"; do
    data_name=$(basename "$dataset" .tar.gz)

    wget "${link}${dataset}"
 
    tar -xvf "$data_name.tar.gz"
    mv "$data_name/$data_name.mtx" "data/$data_name.mtx"

    rm -r "$data_name"
    rm "$data_name.tar.gz"

    python3 exp/get_graphs/mtx_add_edge_weights.py "$data_name"

done

fi


