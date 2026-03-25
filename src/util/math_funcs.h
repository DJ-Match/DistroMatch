/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef _MATH_FUNCS_H_
#define _MATH_FUNCS_H_

#include <iostream>
#include <limits>
#include <stdlib.h>
#include <time.h>

namespace parkec {

template <typename T> inline T min_val() {
    if (std::numeric_limits<T>::is_signed) {
        return -std::numeric_limits<T>::max();
    }

    return 0;
}

double random_weight() { return (double)rand() / (double)RAND_MAX; }

template <typename WeightType> WeightType const_weight() { return 1; }

template <typename T> inline T identity(const T v) { return v; }

/*!
 * \brief Returns true if the rating of e1 is larger than the rating of e2
 *
 * @param e1	first input edge
 * @param e2	second input edge
 * @return True if rating of e1 is larger than rating of e1
 */
template <class Graph>
bool vertex_tie_breaking(const Graph &g, const typename Graph::Edge &e1,
                         const typename Graph::Edge &e2) {
    /// \todo think about the type of the hash, also the return type
    unsigned int e1_n1 = hash(g.local_vertex_id_to_global_id(e1.n1));
    unsigned int e1_n2 = hash(g.local_vertex_id_to_global_id(e1.n2));
    unsigned int e1_max, e1_min;
    if (e1_n1 < e1_n2) {
        e1_max = e1_n2;
        e1_min = e1_n1;
    } else {
        e1_max = e1_n1;
        e1_min = e1_n2;
    }

    unsigned int e2_n1 = hash2(g.local_vertex_id_to_global_id(e2.n1));
    unsigned int e2_n2 = hash2(g.local_vertex_id_to_global_id(e2.n2));
    unsigned int e2_max, e2_min;
    if (e2_n1 < e2_n2) {
        e2_max = e2_n2;
        e2_min = e2_n1;
    } else {
        e2_max = e2_n1;
        e2_min = e2_n2;
    }

    return (e1.weight > e2.weight) ||
           (e1.weight == e2.weight &&
            (e1_max > e2_max || (e1_max == e2_max && e1_min > e2_min)));
}

} // namespace parkec

template <typename NodeIDType>
int proc_id_of_vertex(NodeIDType vid, NodeIDType num_vertices, int procs) {
    int proc_id;

    NodeIDType base_number_of_proc_vertices = num_vertices / procs;
    NodeIDType remainder = num_vertices % procs;

    // Processes with proc-ranges smaller than threshold_id, have
    // local vertex count of base_number_of_proc_vertices+1
    NodeIDType threshold_id = remainder * (base_number_of_proc_vertices + 1);

    if (vid < threshold_id) {
        proc_id = vid / (base_number_of_proc_vertices + 1);
    } else {
        proc_id =
            remainder + ((vid - threshold_id) / base_number_of_proc_vertices);
    }

    return proc_id;
}

template <typename NodeIDType>
int proc_id_of_vertex(
    NodeIDType vid,
    std::vector<std::pair<NodeIDType, NodeIDType>> &proc_ranges) {
    for (NodeIDType p = 0; p < proc_ranges.size(); p++) {
        if (proc_ranges[p].first <= vid && vid < proc_ranges[p].second) {
            return p;
        }
    }

    return -1;
}

#endif
