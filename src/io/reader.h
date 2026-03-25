/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef GRAPH_READER
#define GRAPH_READER

#include <iostream>
#include <limits>
#include <list>
#include <utility>
#include <vector>

#include <mpi.h>

#include <util/hash_functions.h>
#include <util/math_funcs.h>
#include <util/mpi_tags.h>

namespace graph_io {

template <class Graph, typename Graph::NodeIDType (*adjust_id)(
                           const typename Graph::NodeIDType) =
                           parkec::identity<typename Graph::NodeIDType>>
class GraphReader {
    typedef typename Graph::NodeIDType NodeIDType;
    typedef typename Graph::Edge Edge;

  public:
    virtual ~GraphReader() = default;

    virtual void read_header(std::ifstream &in, NodeIDType &n) {}

    virtual void
    read_graph(std::ifstream &in, Graph &g, bool use_partitioning,
               NodeIDType num_vertices, unsigned int proc_id,
               unsigned int procs,
               std::vector<NodeIDType> &first_global_vertex_of_proc,
               std::vector<NodeIDType> &vertex_mapping) {}

    virtual void set_seed(int seed) {}

    int proc_id_of_global_vertex(
        NodeIDType global_id, NodeIDType num_vertices,
        const std::vector<NodeIDType> &first_global_vertex_of_proc) {
        unsigned int procs = first_global_vertex_of_proc.size() - 1;
        // guess proc_id based on equal distribution
        unsigned int proc_id = std::floor(global_id / (num_vertices / procs));
        proc_id = std::min(proc_id, procs - 1);

        // correct guessed proc_id
        while (global_id < first_global_vertex_of_proc[proc_id])
            --proc_id;
        while (global_id >= first_global_vertex_of_proc[proc_id + 1])
            ++proc_id;

        return proc_id;
    }
};

} // namespace graph_io

#endif
