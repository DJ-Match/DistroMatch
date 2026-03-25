import os
import sys

def is_weighted(first_line):
    """
    Determine if the graph is weighted based on the first line of the file.
    The METIS format includes an optional third parameter that indicates weighting.
    """
    parts = first_line.strip().split()
    if len(parts) == 3 and int(parts[2]) == 1:
        return True
    return False

def convert_to_unweighted_metis(input_file, output_file):
    """
    Converts a METIS graph file to an unweighted version, writing the result to output_file.
    If the graph is already unweighted, it simply copies the file.
    Skips lines starting with '%' (comments).
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Read first line, skipping comments
        first_line = infile.readline().strip()
        while first_line.startswith('%'):
            first_line = infile.readline().strip()

        # Check if the graph is weighted
        weighted = is_weighted(first_line)

        # Write the first line without the third value (if weighted)
        if weighted:
            parts = first_line.split()
            # Rewrite first line without the weighting flag
            outfile.write(f"{parts[0]} {parts[1]}\n")
        else:
            outfile.write(first_line + "\n")

        # Process the adjacency list, skipping comment lines and removing weights if necessary
        for line in infile:
            line = line.strip()
            if line.startswith('%'):
                # Skip comment lines
                continue

            parts = line.split()
            if weighted:
                # If weighted, every edge has a pair (vertex, weight), so we take every other element
                vertices = parts[::2]  # Extract only the vertices, ignore weights
                outfile.write(" ".join(vertices) + "\n")
            else:
                # If already unweighted, just write the line as is
                outfile.write(line + "\n")

def process_graphs(input_folder, output_folder):
    """
    Processes all .graph files in the input_folder, converts them to unweighted METIS format,
    and saves them in the output_folder.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Iterate over all files in the input folder
    for filename in os.listdir(input_folder):
        if filename.endswith(".graph"):
            # Skip for rmat graphs with other scale than 20
            if 'rmat' in filename and 's20' not in filename:
                continue
            # Skip for not rmat or fb graphs
            if not ('rmat' in filename or 'fb' in filename) or 'fbB' in filename:
                continue

            input_path = os.path.join(input_folder, filename)
            output_filename = filename.replace(".graph", "_u.graph")
            output_path = os.path.join(output_folder, output_filename)

            print(f"Processing {filename}...")

            if filename.endswith("_u.graph") or filename.endswith("_uw.graph"):
                print("Graph is alrady unweighted")
                continue

            if os.path.isfile(output_path):
                print(f"Unweighted graph {output_filename} alrady exists")
                continue

            # Convert to unweighted METIS format
            convert_to_unweighted_metis(input_path, output_path)

            print(f"Saved unweighted graph to {output_filename}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 metis_create_unweighted.py <input_folder> <output_folder>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    process_graphs(input_folder, output_folder)
