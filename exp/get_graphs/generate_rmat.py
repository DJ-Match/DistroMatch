import argparse
import networkit as nk
import numpy as np

def generate_rmat_graph(scale, edgefactor, p_a, p_b, p_c, p_d, distribution):

    # Generate the R-MAT graph using NetworKit's RmatGenerator
    rmat_gen = nk.generators.RmatGenerator(scale, edgefactor, p_a, p_b, p_c,
            p_d, weighted = True)
    G = rmat_gen.generate()
    # Remove self-loops and multi-edges
    G.removeSelfLoops()
    G.removeMultiEdges()

    # Assign weights to edges based on the selected distribution
    weights = []
    if distribution == 'uni':
        weights = np.random.uniform(size=G.numberOfEdges(),low=1,high=500000)

    elif distribution == 'exp':
        weights = np.clip(np.random.exponential(scale=32000, size=G.numberOfEdges()), 1, 500000)

    weights = weights.astype(int)
    for idx, (u, v) in enumerate(G.iterEdges()):
        G.setWeight(u, v, weights[idx])

    return G

def save_graph_metis(G, filename):
    nk.graphio.writeGraph(G, filename, nk.graphio.Format.METIS)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate an R-MAT graph in METIS format.")
    parser.add_argument("--scale", type=int, required=True, help="Scale (log2 number of vertices)")
    parser.add_argument("--edgefactor", type=int, required=True, help="Edge factor (average edges per vertex)")
    parser.add_argument("--p_a", type=float, required=True, help="R-MAT parameter p_a")
    parser.add_argument("--p_b", type=float, required=True, help="R-MAT parameter p_b")
    parser.add_argument("--p_c", type=float, required=True, help="R-MAT parameter p_c")
    parser.add_argument("--p_d", type=float, required=True, help="R-MAT parameter p_d")
    parser.add_argument("--distribution", choices=['uni', 'exp'], required=True, help="Weight distribution ('uni' for uniform, 'exp' for exponential)")

    args = parser.parse_args()

    # Generate the R-MAT graph
    G = generate_rmat_graph(args.scale, args.edgefactor, args.p_a, args.p_b, args.p_c, args.p_d, args.distribution)

    # Prepare the filename according to the format 'rmat_#scale_#edgefactor_(#a_#b_#c_#d)_#distribution.graph'
    filename = str(f"rmat_s{args.scale}_ef{args.edgefactor}_{args.p_a}_{args.p_b}_{args.p_c}_{args.p_d}_{args.distribution}.graph")

    # Save the graph in METIS format
    save_graph_metis(G, filename)

