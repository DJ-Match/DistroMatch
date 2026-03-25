import os
import sys
import glob

def process_mtx_file(input_file, output_file):
    """
    Convert a Matrix Market (.mtx) file to edge stream format.
    
    Arguments:
    - input_file: Path to the .mtx file.
    - output_file: Path to save the converted .edge_stream file.
    """
    print(f"Converting {input_file} to {output_file}...")

    try:
        with open(input_file, 'r') as f, open(output_file, 'w') as out_f:
            # Skip comments/header lines
            line = f.readline()
            while line.startswith('%'):
                line = f.readline()
            
            # Read matrix dimensions
            num_vertices, num_cols, num_edges = map(int, line.split())
            
            for line in f:
                if line.startswith('%'):
                    continue
                
                data = line.split()
                if len(data) == 3:
                    i, j, w = data
                    w = int(float(w))  # Convert weight to integer if necessary
                else:
                    i, j = data
                    w = 1  # Default weight if none provided
                
                i, j = int(i), int(j)
                
                # Avoid self-loops and duplicate edges by enforcing i < j
                if i != j and i < j:
                    out_f.write(f"{i} {j} {w} 0\n")
        
    except Exception as e:
        print(f"Error processing {input_file}: {e}")


def convert_mtx_to_edge_stream(directory):
    """
    Convert all .mtx files in the directory to edge stream format.
    
    Arguments:
    - directory: Path to the directory containing .mtx files.
    """
    mtx_files = glob.glob(os.path.join(directory, "*.mtx"))
    if not mtx_files:
        print(f"No .mtx files found in {directory}")
        return
    
    for mtx_file in mtx_files:
        output_file = mtx_file.replace(".mtx", ".edge_stream")

        # Skip for rmat graphs with other scale than 20
        if 'rmat' in mtx_file and 's20' not in mtx_file:
            continue

        # Skip if output file already exists
        if os.path.exists(output_file):
            # print(f"Skipping {mtx_file}, output {output_file} already exists.")
            continue
        process_mtx_file(mtx_file, output_file)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 convert_mtx_to_edge_stream.py <input_directory>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a valid directory.")
        sys.exit(1)
    
    convert_mtx_to_edge_stream(input_dir)
