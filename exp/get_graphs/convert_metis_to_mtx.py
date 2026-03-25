import argparse
import os
import glob

def process_metis_file(input_file):
    try:
        # Prepare the output .mtx file
        output_file = os.path.splitext(input_file)[0] + ".mtx"
        
        # Skip if output file already exists
        if os.path.exists(output_file):
            print(f"Skipping {input_file}, output {output_file} already exists.")
            return
        
        with open(input_file, 'r') as f:
            # Read the header
            header = f.readline().strip()
            num_vertices, num_edges, has_weights = map(int, header.split())

        with open(output_file, 'w') as out_f:
            # Write Matrix Market header
            out_f.write("%%MatrixMarket matrix coordinate real symmetric\n")
            out_f.write(f"{num_vertices} {num_vertices} {num_edges}\n")

            # Process each line of the .graph file and write directly
            with open(input_file, 'r') as f:
                f.readline()  # Skip the header line
                for i, line in enumerate(f):
                    if not line.strip():
                        continue
                    elements = line.strip().split()

                    if has_weights == 1:
                        for j in range(0, len(elements), 2):
                            v = int(elements[j])
                            w = float(elements[j+1])
                            if i + 1 < v:  # Only write if i+1 < v for symmetric graphs
                                out_f.write(f"{v} {i + 1} {w}\n")
                    else:
                        for v in elements:
                            v = int(v)
                            if i + 1 < v:  # Only write if i+1 < v for symmetric graphs
                                out_f.write(f"{v} {i + 1} 1.0\n")

        print(f"Conversion to MTX format complete for {input_file}. Output saved in {output_file}")

    except Exception as e:
        print(f"Failed to process {input_file}: {e}")

def process_directory(directory):
    # Find all .graph files in the directory
    graph_files = glob.glob(os.path.join(directory, "*.graph"))

    if not graph_files:
        print(f"No .graph files found in {directory}")
        return

    print(f"Found {len(graph_files)} .graph files in {directory}. Converting...")

    # Process each .graph file
    for graph_file in graph_files:
        process_metis_file(graph_file)

if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="Convert all METIS (.graph) files in a directory to Matrix Market (.mtx) format.")
    parser.add_argument("directory", help="Path to the directory containing .graph files.")

    args = parser.parse_args()

    # Process all .graph files in the directory
    process_directory(args.directory)
