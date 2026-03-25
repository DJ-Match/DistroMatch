/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef K_EDGE_COLORING_H
#define K_EDGE_COLORING_H

#include <iostream>
#include <limits>
#include <list>
#include <utility>
#include <vector>

#include <mpi.h>

#include <util/hash_functions.h>
#include <util/math_funcs.h>
#include <util/mpi_tags.h>

// #define LOGGING

namespace parkec {

static std::vector<char> DEFAULT_VECTOR;

template <class Graph, typename Graph::NodeIDType (*adjust_id)(
                           const typename Graph::NodeIDType) =
                           parkec::identity<typename Graph::NodeIDType>>
class KEdgeColoring {
    typedef typename Graph::Edge Edge;

  public:
    virtual void compute_k_edge_coloring(
        std::vector<std::vector<Edge>> &matchings,
        int &total_global_matching_size, int &total_number_of_rounds,
        double &total_global_matching_weight, Graph &g, unsigned int k) {}

    void set_return_local(bool return_local_edges) {
        m_return_local_edges = return_local_edges;
    }

#ifdef TESTING
    void set_test_param_1(bool test_param_1) { m_test_param_1 = test_param_1; }
    void set_test_param_2(bool test_param_2) { m_test_param_2 = test_param_2; }
#endif

  protected:
#ifdef TESTING
    bool m_test_param_1;
    bool m_test_param_2;
#endif
    bool m_return_local_edges;
};

} // namespace parkec

#endif
