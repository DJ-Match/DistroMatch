/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef METIS_GRAPH_READER
#define METIS_GRAPH_READER

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <mpi.h>

#include <config.h>
#include <datastructure/parallel_edge_graph.h>
#include <io/reader.h>
#include <util/mpi_tags.h>
#include <util/random.h>

namespace graph_io {

template <class Graph, typename Graph::NodeIDType (*adjust_id)(
                           const typename Graph::NodeIDType) =
                           parkec::identity<typename Graph::NodeIDType>>
class MetisGraphReader : public GraphReader<Graph> {
    typedef typename Graph::WeightType WeightType;
    typedef typename Graph::NodeIDType NodeIDType;
    typedef typename Graph::EdgeIDType EdgeIDType;
    typedef typename Graph::Edge Edge;

    void read_header(std::ifstream &in, NodeIDType &n) {
        size_t m = 0;
        n = 0;

        std::string line;
        std::getline(in, line);

        // ignore comment line
        while (line[0] == '%') {
            std::getline(in, line);
        }

        std::stringstream strs(line);
        // read header
        strs >> n >> m >> m_fmt;
    }

    void read_graph(std::ifstream &in, Graph &g, bool use_partitioning,
                    NodeIDType num_vertices, unsigned int proc_id,
                    unsigned int procs,
                    std::vector<NodeIDType> &first_global_vertex_of_proc,
                    std::vector<NodeIDType> &vertex_mapping) {
        NodeIDType src_vertex_idx = 0;

        bool has_edge_weights = m_fmt == 1;
        std::string line;

        while (std::getline(in, line)) {
            if (line[0] == '%') {
                // ignore comment line
                continue;
            }
            auto src_vertex = use_partitioning ? vertex_mapping[src_vertex_idx]
                                               : src_vertex_idx;

            auto src_is_local =
                is_local(src_vertex, proc_id, first_global_vertex_of_proc);

            if (!src_is_local) {
                ++src_vertex_idx;
                continue;
            }

            // read edges
            std::stringstream strs(line);
            NodeIDType tar_vertex;
            while (strs >> tar_vertex) {
                double weight = 1;
                if (has_edge_weights) {
                    strs >> weight;
                }

                --tar_vertex; // start node count at 0!
                if (use_partitioning)
                    tar_vertex = vertex_mapping[tar_vertex];

                auto tar_is_local =
                    is_local(tar_vertex, proc_id, first_global_vertex_of_proc);

                auto e = Edge(src_vertex, tar_vertex, weight);

                if (tar_is_local) {
                    if (tar_vertex < src_vertex) {
                        continue;
                    }
                    g.add_local_edge(e);
                } else {
                    g.add_cross_edge(e, this->proc_id_of_global_vertex(
                                            tar_vertex, num_vertices,
                                            first_global_vertex_of_proc));
                }
            }
            ++src_vertex_idx;
        }
    }

  private:
    bool is_local(NodeIDType id, int proc_id,
                  const std::vector<NodeIDType> &first_global_vertex_of_proc) {
        return id >= first_global_vertex_of_proc[proc_id] &&
               id < first_global_vertex_of_proc[proc_id + 1];
    };

    int m_fmt = 0;
};

} // namespace graph_io

#endif
