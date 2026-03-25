/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef _MIN_PARALLEL_EDGE_GRAPH_
#define _MIN_PARALLEL_EDGE_GRAPH_

#include <algorithm>
#include <iostream>
#include <limits>
#include <list>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"

#include <mpi.h>

#include "edge.h"
#include <datastructure/parallel_edge_graph.h>
#include <util/math_funcs.h>
#include <util/mpi_tags.h>

namespace parkec {

/// _NodeIDType should be at least unsigned long if there are more than 2^32 - 1
/// nodes _EdgeIDType should be at least unsigned long if there are more than
/// 2^32 - 1 edges
template <typename _WeightType = double, typename _NodeIDType = unsigned int,
          typename _EdgeIDType = long long>
class MinParallelEdgeGraph : public ParallelEdgeGraph<> {

  public:
    typedef _WeightType WeightType;
    typedef _NodeIDType NodeIDType;
    typedef _EdgeIDType EdgeIDType;

    typedef parkec::Edge<WeightType, NodeIDType> Edge;

    typedef typename std::vector<Edge>::iterator edge_iterator;
    typedef typename std::vector<Edge>::const_iterator const_edge_iterator;

    typedef std::list<NodeIDType> ConnectedComponent;
    typedef std::list<ConnectedComponent> ConnectedComponents;

  public:
    MinParallelEdgeGraph() {
        m_initialized = false;
        m_dummy_vertex = std::numeric_limits<NodeIDType>::max();
        m_dummy_edge.n1 = std::numeric_limits<NodeIDType>::max();
        m_dummy_edge.n2 = std::numeric_limits<NodeIDType>::max();
        m_dummy_edge.weight = min_val<WeightType>();

        m_num_local_edges = 0;
        m_num_cross_edges = 0;

        MPI_Comm_size(MPI_COMM_WORLD, &m_procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &m_proc_id);
    }

    inline edge_iterator begin_active_edges() override {
        // printf("get begin_active_edges: ")
        return m_local_edges.begin() + begin_edges();
    }
    inline edge_iterator end_active_edges() override {
        return m_local_edges.begin() + m_end_active_local_edges;
    }

    inline edge_iterator deactivate_edge(const edge_iterator edge_it) override {
        EdgeIDType edge_ref = edge_it - m_local_edges.begin();

        if (edge_ref >= m_end_active_local_edges) {
            return edge_it;
        }

        --m_end_active_local_edges;

        swap(m_end_active_local_edges, edge_ref, m_local_edges);

        return m_local_edges.begin() + m_end_active_local_edges;
    }

    inline edge_iterator
    edge_is_candidate(const edge_iterator edge_it) override {
        EdgeIDType edge_ref = edge_it - m_local_edges.begin();

        if (edge_ref < m_end_active_candidates) {
            return edge_it;
        }

        if (m_end_active_candidates < edge_ref)
            swap(m_end_active_candidates, edge_ref, m_local_edges);

        ++m_end_active_candidates;

        return m_local_edges.begin() + (m_end_active_candidates - 1);
    }

    inline edge_iterator end_active_candidates() override {
        return m_local_edges.begin() + m_end_active_candidates;
    }
    inline void reset_candidates() override {
        m_end_active_candidates = begin_edges();
    }

    inline edge_iterator
    deactivate_candidate(const edge_iterator edge_it) override {
        EdgeIDType edge_ref = edge_it - m_local_edges.begin();

        if (edge_ref >= m_end_active_candidates) {
            return edge_it;
        }
        --m_end_active_candidates;

        swap(m_end_active_candidates, edge_ref, m_local_edges);

        return m_local_edges.begin() + m_end_active_candidates;
    }

    inline NodeIDType get_ghost_id(const NodeIDType n_id) const override {

        if (m_ghost_global_to_local_hash.contains(n_id))
            return global_vertex_id_to_local_id_of_ghost_vertex(n_id);
        else {
            return -1;
        }
    }

    inline NodeIDType
    global_vertex_id_to_local_id(const NodeIDType n_id) const override {
        if (global_id_is_local(n_id)) {
            return global_vertex_id_to_local_id_of_local_vertex(n_id);
        } else {
            return get_ghost_id(n_id);
        }
    }

    inline EdgeIDType get_edge_id(const edge_iterator &edge_it) const override {
        return edge_it - m_local_edges.begin();
    }
    inline EdgeIDType
    get_edge_id(const const_edge_iterator &edge_it) const override {
        return edge_it - m_local_edges.begin();
    }

    inline void add_local_edge(Edge &e) override {
        ++m_num_local_edges;
        add_edge(e);
    }
    inline void add_cross_edge(Edge &e, int proc_of_ghost) override {
        ++m_num_cross_edges;

        if (m_active_cross_edges_of_proc[proc_of_ghost] == 0) {
            ++m_count_active_partners;
        }

        ++m_active_cross_edges_of_proc[proc_of_ghost];

        add_edge(e, proc_of_ghost);
    }

    inline void add_edge(Edge &e, int proc_of_ghost = -1) {
        if (m_initialized) {
            printf("Graph is already final!");
            exit(1);
        }
        if (!global_id_is_local(e.n1)) {
            if (!global_id_is_local(e.n2)) {
                // This is not my edge!
                printf("[%d] Edge (%d, %d) does not belong the assigned vertex "
                       "range.",
                       m_proc_id, e.n1, e.n2);
                throw std::logic_error("Initialization error");
            }

            // swap
            auto tmp = e.n1;
            e.n1 = e.n2;
            e.n2 = tmp;
        }

        if (!global_id_is_local(e.n2)) {

            auto local_ghost_id = -1;

            if (m_ghost_global_to_local_hash.contains(e.n2)) {
                local_ghost_id = m_ghost_global_to_local_hash.at(e.n2);

            } else {
                // New ghost vertex
                local_ghost_id = m_num_local_vertices + m_num_ghost_vertices;
                m_ghost_global_to_local_hash.insert(
                    std::make_pair(e.n2, local_ghost_id));
                m_local_ghost_to_global.push_back(e.n2);
                m_proc_of_ghost_vertex.push_back(proc_of_ghost);

                ++m_num_ghost_vertices;
            }

            e.n2 = local_ghost_id;

        } else {
            e.n2 -= m_first_global_vertex;
        }

        // Adjust local ids
        e.n1 -= m_first_global_vertex;

        if (e.weight > m_max_weight) {
            m_max_weight = e.weight;
        }

        m_local_edges.push_back(e);
    }

    inline void finalize_initialization() override {
        m_initialized = true;

        m_end_active_local_edges = m_local_edges.size();
        m_end_active_cross_edges = 0;
        m_end_active_candidates = 0;

        m_start_low_weight_edges = 0;

        // add one dummy node
        m_local_ghost_to_global.push_back(m_dummy_vertex);

        m_local_edges.push_back(m_dummy_edge);
        m_cross_edges.push_back(m_dummy_edge);
    }

    inline void activate_edges() override {
        m_end_active_local_edges = m_local_edges.size() - 1;
    }

    void virtual reset() override {
        activate_edges();

        // deactivate colored edges
        auto e_it = begin_active_cross_edges();
        while (e_it != end_active_cross_edges()) {
            if (e_it->color < color_set::not_colorable) {
                deactivate_cross_edge(e_it);
            } else {
                ++e_it;
            }
        }

        e_it = begin_active_local_edges();
        while (e_it < end_active_local_edges()) {
            if (e_it->color < color_set::not_colorable) {
                deactivate_local_edge(e_it);
            } else {
                ++e_it;
            }
        }

        e_it = begin_active_edges();
        while (e_it != end_active_edges()) {
            NodeIDType g_n1 = local_vertex_id_to_global_id(e_it->n1);
            NodeIDType g_n2 = local_vertex_id_to_global_id(e_it->n2);

            auto hash_n1 = hash(g_n1);
            auto hash_n2 = hash(g_n2);

            if ((!is_local(e_it->n2)) && (hash_n1 > hash_n2)) {
                deactivate_edge(e_it);
                continue;
            }
            ++e_it;
        }
    }

    inline virtual bool has_active_edges() const {
        return m_end_active_local_edges - begin_edges();
    }

    inline void sort_high_low(WeightType weight) {

        m_theshold_weight = weight;
        auto e_it = begin_active_edges();
        while (e_it < end_active_edges()) {
            if (e_it->weight <= weight) {
                ++m_num_low_edges;
                deactivate_edge(e_it);
                continue;
            }
            ++e_it;
        }

        m_start_low_weight_edges = m_end_active_local_edges;
    }
    inline void swap_to_low() {

        use_low = true;
        m_end_active_local_edges = m_start_low_weight_edges + m_num_low_edges;
    }
    inline WeightType get_max_weight() const { return m_max_weight; }

    inline void print() const override {
        for (EdgeIDType pos = 0; pos < m_num_local_edges + m_num_cross_edges;
             pos++) {
            const auto &e = m_local_edges[pos];
            printf("(%d<->%d, %f),", local_vertex_id_to_global_id(e.n1),
                   local_vertex_id_to_global_id(e.n2), e.weight);
        }
        printf("\n");

        printf("use_low%d, theshold_weight: %.2f,max_weight: %.2f, "
               "low_start:%lld, low_num:%lld, end_active: %lld\n",
               use_low, m_theshold_weight, m_max_weight,
               m_start_low_weight_edges, m_num_low_edges,
               m_end_active_local_edges);
    }

    inline EdgeIDType getNumHighWeightEdges() const {
        return m_start_low_weight_edges;
    }

  private:
    inline EdgeIDType begin_edges() const {
        return use_low ? m_start_low_weight_edges : 0;
    }

    EdgeIDType m_start_low_weight_edges = 0;
    EdgeIDType m_num_low_edges = 0;
    bool use_low = 0;
    WeightType m_max_weight = 0;
    WeightType m_theshold_weight = 0;
};

} // namespace parkec

#endif
