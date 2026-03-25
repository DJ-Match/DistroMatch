import os
import sys

def process_large_graph(input_file, output_file):
    """
    Process the input graph file line by line and write the modified graph
    with node degrees as vertex weights to the output file.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Read the header line
        header = infile.readline().strip().split()
        num_vertices = int(header[0])
        num_edges = int(header[1])
        fmt = int(header[2]) if len(header) > 2 else 0

        # Check if the graph already has vertex weights
        has_vertex_weights = fmt == 10 or fmt == 11
        has_edge_weights = fmt == 1 or fmt == 11

        if has_vertex_weights:
            print(f"Skipping {input_file} (already has node weights)")
            return

        # Initialize node degrees and updated format
        node_degrees = [0] * num_vertices
        new_fmt = 11 if has_edge_weights else 10  # Vertex weights added, keep edge weights if present

        # First write the header to the output file
        # num_vertices num_edges new_fmt
        outfile.write(f"{num_vertices} {num_edges} {new_fmt}\n")

        # Process each line of the graph, calculate degrees, and write the updated graph
        for vertex_idx in range(num_vertices):
            line = infile.readline().strip()
            edges = list(map(int, line.split()))

            # Calculate the degree (number of neighbors)
            node_degrees[vertex_idx] = len(edges) // (2 if has_edge_weights else 1)

            # Write vertex degree and edges to the output file
            outfile.write(f"{node_degrees[vertex_idx]} {line}\n")


def main(folder_path):
    """
    Find all files with '_u.graph' suffix, process them, and save the result
    with '_uw.graph' suffix.
    """
    # Iterate over all files in the directory
    for file_name in os.listdir(folder_path):
        if file_name.endswith("_u.graph"):
            input_file = os.path.join(folder_path, file_name)
            output_file = input_file.replace("_u.graph", "_uw.graph")

            if os.path.isfile(output_file):
                print(f"Target graph {output_file} alrady exists")
                continue

            print(f"Processing {input_file} -> {output_file}")
            process_large_graph(input_file, output_file)
            print(f"Finished processing {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 create_degree_node_weight_graphs.py <folder_path>")
        sys.exit(1)

    folder_path = sys.argv[1]
    if not os.path.isdir(folder_path):
        print(f"Error: {folder_path} is not a valid directory.")
        sys.exit(1)

    main(folder_path)

