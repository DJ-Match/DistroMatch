/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef MTX_GRAPH_READER
#define MTX_GRAPH_READER

#include <sstream>
#include <string>
#include <vector>

#include <datastructure/parallel_edge_graph.h>
#include <io/reader.h>

namespace graph_io {

template <class Graph, typename Graph::NodeIDType (*adjust_id)(
                           const typename Graph::NodeIDType) =
                           parkec::identity<typename Graph::NodeIDType>>
class MtxGraphReader : public GraphReader<Graph> {
    typedef typename Graph::WeightType WeightType;
    typedef typename Graph::NodeIDType NodeIDType;
    typedef typename Graph::EdgeIDType EdgeIDType;
    typedef typename Graph::Edge Edge;

  public:
    void read_header(std::ifstream &in, NodeIDType &n) {
        NodeIDType m;
        read_header_with_m(in, n, m);
    }
    void read_header_with_m(std::ifstream &in, NodeIDType &n, NodeIDType &m) {
        std::string line;
        std::getline(in, line);

        // read type code
        const std::string weightedUndirectedHeader =
            "%%MatrixMarket matrix coordinate integer symmetric";
        const std::string unweightedUndirectedHeader =
            "%%MatrixMarket matrix coordinate pattern symmetric";
        const std::string weightedRealUndirectedHeader =
            "%%MatrixMarket matrix coordinate real symmetric";

        m_weighted = line.compare(unweightedUndirectedHeader) != 0;

        if (!m_weighted && line.compare(unweightedUndirectedHeader) != 0 &&
            line.compare(weightedRealUndirectedHeader) != 0)
            throw std::invalid_argument("Illegal argument: Header does not "
                                        "match expected MatrixMarket formats.");

        // ignore comment line
        do {
            std::getline(in, line);
        } while (line[0] == '%');

        std::stringstream strs(line);
        // read header (format rows, columns, entries). #rows=#cols as the
        // matrix is symmetric
        int cols;
        strs >> n >> cols >> m;
    }

    void read_graph(std::ifstream &in, Graph &g, bool use_partitioning,
                    NodeIDType num_vertices, unsigned int proc_id,
                    unsigned int procs,
                    std::vector<NodeIDType> &first_global_vertex_of_proc,
                    std::vector<NodeIDType> &vertex_mapping) {
        NodeIDType src_vertex, tar_vertex;

        std::string line;
        unsigned int self_loops = 0;

        while (std::getline(in, line)) {
            if (line[0] == '%') {
                // ignore comment line
                continue;
            }

            std::stringstream strs(line);

            strs >> src_vertex >> tar_vertex;

            // Ids should start with zero
            --src_vertex;
            --tar_vertex;

            double weight = 1;
            if (m_weighted) {
                strs >> weight;
            }

            if (use_partitioning) {
                src_vertex = vertex_mapping[src_vertex];
                tar_vertex = vertex_mapping[tar_vertex];
            }

            auto src_is_local = g.global_id_is_local(src_vertex);
            auto tar_is_local = g.global_id_is_local(tar_vertex);

            if (!(src_is_local || tar_is_local)) {
                continue;
            }

            if (src_vertex == tar_vertex) {
                // ignore self edges
                ++self_loops;
                continue;
            }

            auto e = Edge(src_vertex, tar_vertex, weight);
            if (src_is_local && tar_is_local) {
                g.add_local_edge(e);
            } else {
                g.add_cross_edge(e, this->proc_id_of_global_vertex(
                                        src_is_local ? tar_vertex : src_vertex,
                                        num_vertices,
                                        first_global_vertex_of_proc));
            }
        }

        if (self_loops)
            printf("[%d] Ignored %d self loops\n", proc_id, self_loops);
    }

  private:
    int determine_ownership(std::vector<unsigned int> &vertex_ownership,
                            std::vector<unsigned int> &edge_counts,
                            NodeIDType u, NodeIDType v, unsigned int proc_id,
                            unsigned int procs) {

        auto u_owned = vertex_ownership[u] != no_owner;
        auto v_owned = vertex_ownership[v] != no_owner;
        if (u_owned) {
            ++edge_counts[vertex_ownership[u]];
        }
        if (v_owned) {
            ++edge_counts[vertex_ownership[u]];
        }

        if (u_owned && v_owned) {
            // No new node for me -> return wether its ny edge
            return 0;
        }

        if (!u_owned | !v_owned) {

            // proc with fewest edges + lowest id takes it!

            unsigned int best_proc = 0;
            unsigned int min_edges = edge_counts[0];

            for (unsigned int p = 1; p < procs; ++p) {
                if (edge_counts[p] < min_edges) {
                    best_proc = p;
                    min_edges = edge_counts[p];
                }
            }

            vertex_ownership[u] = best_proc;
            vertex_ownership[v] = best_proc;
            edge_counts[best_proc] += 2;

            return 2;
        }

        auto owned_vertex = u_owned ? u : v;
        auto unowned_vertex = u_owned ? v : u;

        auto owned_is_mine = vertex_ownership[owned_vertex] == proc_id;

        if (owned_is_mine /* && balance is okay */) {
            // Own unowned_vertex if balance is ok
        }

        if (true /* && my balance is off */) {
            // own vertex anyway
        }

        return 1;
    }

    unsigned int no_owner = -1;

  protected:
    bool m_weighted = 0;
};

} // namespace graph_io

#endif
