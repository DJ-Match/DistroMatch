/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef GRAPH_H
#define GRAPH_H

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
#include <util/color_set.h>
#include <util/math_funcs.h>
#include <util/mpi_tags.h>

namespace parkec {
class BSuitor;

/// _NodeIDType should be at least unsigned long if there are more than 2^32 - 1
/// nodes _EdgeIDType should be at least unsigned long if there are more than
/// 2^32 - 1 edges
template <typename _WeightType = double, typename _NodeIDType = unsigned int,
          typename _EdgeIDType = long long>
class IGraph {
    friend class BSuitor;

  public:
    typedef _WeightType WeightType;
    typedef _NodeIDType NodeIDType;
    typedef _EdgeIDType EdgeIDType;

    typedef parkec::Edge<WeightType, NodeIDType> Edge;

    typedef typename std::vector<Edge>::iterator edge_iterator;
    typedef typename std::vector<Edge>::const_iterator const_edge_iterator;

  public:
    IGraph() {
        MPI_Comm_size(MPI_COMM_WORLD, &m_procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &m_proc_id);
        m_initialized = false;
    }

    virtual ~IGraph() {}

    inline NodeIDType num_all_vertices() const {
        return m_num_local_vertices + m_num_ghost_vertices;
    }

    inline NodeIDType num_local_vertices() const {
        return m_num_local_vertices;
    }

    inline NodeIDType num_ghost_vertices() const {
        return m_num_ghost_vertices;
    }

    inline virtual EdgeIDType getNumAllEdges() const {
        return m_num_local_edges + m_num_cross_edges;
    }
    inline virtual EdgeIDType getNumLocalEdges() const {
        return m_num_local_edges;
    }
    inline virtual EdgeIDType getNumCrossEdges() const {
        return m_num_cross_edges;
    }

    inline bool is_ghost(const NodeIDType n_id) const {
        return n_id - m_num_local_vertices < m_num_ghost_vertices;
    }

    inline NodeIDType num_global_vertices() const {
        return m_num_global_vertices;
    }

    inline NodeIDType first_global_vertex() const {
        return m_first_global_vertex;
    }

    inline NodeIDType end_global_vertex() const { return m_end_global_vertex; }

    inline NodeIDType
    global_vertex_id_to_local_id_of_local_vertex(const NodeIDType n_id) const {
        return n_id - m_first_global_vertex;
    }

    inline NodeIDType
    global_vertex_id_to_local_id_of_ghost_vertex(const NodeIDType n_id) const {
        return m_ghost_global_to_local_hash.at(n_id);
    }

    inline bool global_id_is_local(const NodeIDType n_id) const {
        return n_id >= m_first_global_vertex && n_id < m_end_global_vertex;
    }

    inline bool is_local(const NodeIDType n_id) const {
        return n_id < m_num_local_vertices;
    }

    inline virtual NodeIDType
    global_vertex_id_to_local_id(const NodeIDType n_id) const {
        if (global_id_is_local(n_id)) {
            return global_vertex_id_to_local_id_of_local_vertex(n_id);
        } else {
            return global_vertex_id_to_local_id_of_ghost_vertex(n_id);
        }
    }

    inline int get_proc_of_ghost_vertex(const NodeIDType n_id) const {
        return m_proc_of_ghost_vertex[n_id - m_num_local_vertices];
    }

    inline NodeIDType
    local_vertex_id_to_global_id(const NodeIDType n_id) const {

        return is_local(n_id)
                   ? m_first_global_vertex + n_id
                   : m_local_ghost_to_global[n_id - m_num_local_vertices];
    }

    inline int get_active_partner_count() const {
        return m_count_active_partners;
    }

    inline virtual void
    start_initialization(std::vector<NodeIDType> &first_global_vertex_of_proc) {

        if (m_initialized) {
            printf("Graph is already final!");
            exit(1);
        }

        m_num_global_vertices = first_global_vertex_of_proc[m_procs];
        m_first_global_vertex = first_global_vertex_of_proc[m_proc_id];
        m_end_global_vertex = first_global_vertex_of_proc[m_proc_id + 1];

        m_num_local_vertices = m_end_global_vertex - m_first_global_vertex;
        m_num_ghost_vertices = 0;

        m_active_cross_edges_of_proc.resize(m_procs, 0);
        m_count_active_partners = 0;
    }

    inline virtual bool is_active_partner(const int proc_id) const {
        return m_active_cross_edges_of_proc[proc_id];
    }

    inline void reduce_active_edges_of_proc(int p) {
        // printf("[%d] Reduce NOW: m_active_cross_edges_of_proc[p%d]= %d\n",
        //    m_proc_id, p, m_active_cross_edges_of_proc[p] - 1);

        // .n2 is ghost - we ensure this in the constructor
        --m_active_cross_edges_of_proc[p];

        m_count_active_partners -= (m_active_cross_edges_of_proc[p] == 0);
    }

    inline virtual void add_local_edge(Edge &e) = 0;

    inline virtual void add_cross_edge(Edge &e, int proc_of_ghost) = 0;

    inline virtual void finalize_initialization() { m_initialized = true; }

    bool virtual valid_k_matchings(std::vector<std::vector<Edge>> &matchings,
                                   bool check_maximal = true) = 0;

  protected:
    // Helper to send and receive messages with MPI
    template <typename MSG>
    void send_receive(std::vector<std::pair<std::vector<MSG>, MPI_Request>>
                          &messages_for_proc,
                      MPI_Datatype mpiDatatype, int tag,
                      std::function<void(MSG &)> process_msg) {

        for (int p = 0; p < m_procs; p++) {
            if (p != m_proc_id) {
                if (!messages_for_proc[p].first.empty()) {
                    MPI_Isend(&messages_for_proc[p].first[0],
                              messages_for_proc[p].first.size(), mpiDatatype, p,
                              tag, MPI_COMM_WORLD,
                              &messages_for_proc[p].second);
                } else {
                    MPI_Isend(0, 0, mpiDatatype, p, tag, MPI_COMM_WORLD,
                              &messages_for_proc[p].second);
                }
            }
        }

        // Receive messages
        for (int p = 0; p < m_procs; p++) {
            if (p != m_proc_id) {
                MPI_Status status;
                MPI_Probe(p, tag, MPI_COMM_WORLD, &status);

                int count;
                MPI_Get_count(&status, mpiDatatype, &count);

                if (count > 0) {
                    MSG *messages = new MSG[count];
                    MPI_Recv(messages, count, mpiDatatype, status.MPI_SOURCE,
                             tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    for (int i = 0; i < count; i++) {
                        process_msg(messages[i]);
                    }
                }
            }
        }

        for (unsigned int p = 0; p < messages_for_proc.size(); p++) {
            if (!messages_for_proc[p].first.empty()) {
                MPI_Wait(&messages_for_proc[p].second, MPI_STATUS_IGNORE);
                messages_for_proc[p].first.clear();
            }
        }
    }

    NodeIDType m_first_global_vertex;
    NodeIDType m_end_global_vertex;
    NodeIDType m_num_global_vertices;

    NodeIDType m_num_local_vertices;
    NodeIDType m_num_ghost_vertices;

    EdgeIDType m_num_local_edges;
    EdgeIDType m_num_cross_edges;

    std::vector<int> m_proc_of_ghost_vertex;

    std::vector<NodeIDType> m_local_ghost_to_global;

    absl::flat_hash_map<NodeIDType, NodeIDType> m_ghost_global_to_local_hash;

    std::vector<EdgeIDType> m_active_cross_edges_of_proc;

    int m_count_active_partners;

    int m_proc_id, m_procs;

    bool m_initialized;
};

} // namespace parkec

#endif
