/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef BIN_GRAPH_READER
#define BIN_GRAPH_READER

#include <sstream>
#include <string>
#include <vector>

#include <datastructure/parallel_edge_graph.h>
#include <io/reader.h>

namespace graph_io {

template <class Graph, typename Graph::NodeIDType (*adjust_id)(
                           const typename Graph::NodeIDType) =
                           parkec::identity<typename Graph::NodeIDType>>
class BinGraphReader : public GraphReader<Graph> {
    typedef typename Graph::WeightType WeightType;
    typedef typename Graph::NodeIDType NodeIDType;
    typedef typename Graph::EdgeIDType EdgeIDType;
    typedef typename Graph::Edge Edge;

  public:
    void read_header(std::ifstream &in, NodeIDType &n) {
        int numEdges;

        // Read the header (n and m)
        in.read(reinterpret_cast<char *>(&n), sizeof(int));
        in.read(reinterpret_cast<char *>(&numEdges), sizeof(int));
    }

    void read_graph(std::ifstream &in, Graph &g, bool use_partitioning,
                    NodeIDType num_vertices, unsigned int proc_id,
                    unsigned int procs,
                    std::vector<NodeIDType> &first_global_vertex_of_proc,
                    std::vector<NodeIDType> &vertex_mapping) {
        NodeIDType src_vertex, tar_vertex;
        double weight;

        unsigned int self_loops = 0;

        while (in.read(reinterpret_cast<char *>(&src_vertex), sizeof(int))) {
            in.read(reinterpret_cast<char *>(&tar_vertex), sizeof(int));
            in.read(reinterpret_cast<char *>(&weight), sizeof(double));

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

  protected:
    bool m_weighted = 0;
};

} // namespace graph_io

#endif
