/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef DCCM_H
#define DCCM_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

#include <mpi.h>

#include <datastructure/min_parallel_edge_graph.h>
#include <util/color_set.h>
#include <util/hash_functions.h>
#include <util/math_funcs.h>
#include <util/mpi_tags.h>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"

#include "k_edge_coloring.h"

// #define LOGGING
// #define LOGGING_V
// #define LOGGING_HASH_MAP
// #define LOGGING_L0
// #define DEBUG
// #define USE_HASH_MAP

// #define MONITORING

namespace parkec {

template <class GraphBase> class DirectCCM : public KEdgeColoring<GraphBase> {
    typedef typename parkec::MinParallelEdgeGraph<> Graph;
    typedef typename Graph::WeightType WeightType;
    typedef typename Graph::NodeIDType NodeIDType;
    typedef typename Graph::EdgeIDType EdgeIDType;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::edge_iterator edge_iterator;
    typedef typename Graph::const_edge_iterator const_edge_iterator;

    struct EdgeMinimal {
        int first;
        int second;
        int color;
    };

  public:
    DirectCCM(double ccm_epsilon) : m_edge_weight_threshold(1 - ccm_epsilon) {
        MPI_Comm_size(MPI_COMM_WORLD, &m_procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &m_proc_id);
    }

    void compute_k_edge_coloring(std::vector<std::vector<Edge>> &matchings,
                                 int &total_global_matching_size,
                                 int &total_number_of_rounds,
                                 double &total_global_matching_weight,
                                 GraphBase &g_ref, unsigned int k) {
        Graph &g = static_cast<Graph &>(g_ref);

#ifdef LOGGING_L0
        printf(
            "[%d] Global node range: [%d,%d[, num local edges: %lu, num cross "
            "edges: %lu\n",
            m_proc_id, g.first_global_vertex(), g.end_global_vertex(),
            g.getNumLocalEdges(), g.getNumCrossEdges());

        printf("[%d] num local vertices: %d, num ghost vertices: %d\n",
               m_proc_id, g.num_local_vertices(), g.num_ghost_vertices());

        // printf("[%d] CrossCandidate: size %ld, align %ld; Candidate:
        // size %ld, align
        // "
        //        "%ld\n",
        //        m_proc_id, sizeof(CrossCandidate),
        //        alignof(CrossCandidate), sizeof(Candidate),
        //        alignof(Candidate));

#endif
        int local_size = g.num_local_vertices();
        // Gather sizes from all processes
        std::vector<int> vertices_per_proc(m_procs);
        MPI_Allgather(&local_size, 1, MPI_INT, vertices_per_proc.data(), 1,
                      MPI_INT, MPI_COMM_WORLD);

        // Compute displacements
        std::vector<int> displacements(m_procs);
        m_global_number_of_vertices = vertices_per_proc[0];
        displacements[0] = 0;
        for (int i = 1; i < m_procs; ++i) {
            m_global_number_of_vertices += vertices_per_proc[i];
        }

        m_global_to_local.resize(m_global_number_of_vertices);

        for (NodeIDType i = 0; i < m_global_number_of_vertices; ++i) {
            m_global_to_local[i] = g.global_vertex_id_to_local_id(i);
        }

        // printf("SIZE colorset: %ld, nops: %d\n", sizeof(color_set::bit_type),
        //    color_set::npos);

        // buffers for messages that are used to send candidate requests of
        // ghostnodes
        std::vector<EdgeMinimal> candidates;
        candidates.resize(m_procs);

        matchings.resize(k, std::vector<Edge>(0));
        m_color_available.resize(g.num_all_vertices(), color_set(k));
        m_k = k;

        m_candidate_weight_of_node.resize(g.num_all_vertices(), 0);

        edge_iterator e_it = g.begin_active_edges();

        while (e_it != g.end_active_edges()) {
            NodeIDType g_n1 = g.local_vertex_id_to_global_id(e_it->n1);
            NodeIDType g_n2 = g.local_vertex_id_to_global_id(e_it->n2);

            auto hash_n1 = hash(g_n1);
            auto hash_n2 = hash(g_n2);

            if ((!g.is_local(e_it->n2)) && (hash_n1 > hash_n2)) {
                // printf("[%d] (c) Deactivated edge (v%d (%d), v%d (%d)) [%d >
                // "
                //        "%d]\n",
                //        m_proc_id, g_n1, e_it->n1, g_n2, e_it->n2,
                //        hash_n1, hash_n2);

                deactivate_edge(g, e_it);
                continue;
            }
            ++e_it;
        }

#ifdef USE_HASH_MAP
        while (e_it != g.end_active_edges()) {

            if ((!g.is_local(e_it->n2))) {

                m_cross_edge_id.insert(std::make_pair(
                    std::make_pair(e_it->n1, e_it->n2), g.get_edge_id(e_it)));
            }
            ++e_it;
        }
#endif

        int round = 0;
        m_comm_rounds = 0;
        m_comm_wait_rounds = 0;
#ifdef MONITORING
        double start_round_time = MPI_Wtime();
        double set_weights_time, handle_local_time, handle_cross_time,
            deactivate_time;
#endif

        // Mersenne Twister RNG with the provided seed
        std::mt19937 rng(m_proc_id);

        while (g.has_active_edges()) {

#ifdef LOGGING_L0
            printf("[%d] --- --- Start round %d --- ---\n", m_proc_id, round);
            auto edges_active_local =
                g.end_active_local_edges() - g.begin_active_local_edges();

            printf("[%d] Start round %d: edges_active %ld", m_proc_id, round,
                   edges_active_local);
#endif

#ifdef MONITORING
            m_asked_candidates_per_round = 0;

            auto edges_active_local =
                g.end_active_local_edges() - g.begin_active_local_edges();
            auto edges_active_cross =
                g.end_active_cross_edges() - g.begin_active_cross_edges();

            printf("[%d] Start round %d: edges_active_local %ld, "
                   "edges_active_cross %ld\n",
                   m_proc_id, round, edges_active_local, edges_active_cross);

#endif
            std::vector<NodeIDType> nonzero_candidates;

            set_candidate_weight(g, nonzero_candidates);

// Use in-place MPI_Allgatherv (send and receive buffer are the
// same)
// MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
//                m_candidate_weight_of_node.data(),
//                vertices_per_proc.data(), displacements.data(),
//                MPI_DOUBLE, MPI_COMM_WORLD);

// In-place reduction: Every process gets the element-wise maximum

// MPI_Allreduce(MPI_IN_PLACE, m_candidate_weight_of_node.data(),
//   total_size, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

// double start_time = MPI_Wtime();
#ifdef LOGGING_L0
            auto candidates_communicated =
#endif
                sparse_global_max(g, nonzero_candidates);
            // printf("[%d] Round %d: sparse_global_max took %.6f s\n",
            // m_proc_id,
            //    round, MPI_Wtime() - start_time);

            g.reset_candidates();

            edge_iterator e_it = g.begin_active_edges();
            while (e_it < g.end_active_edges()) {

                if (is_candidate(g, e_it)) {
                    auto cc = find_free_color(g, e_it, matchings, rng);

                    if (cc == color_set::not_colorable) {
                        deactivate_candidate(g, e_it);
                        continue;
                    }

                    // swap edges to front
                    edge_is_candidate(g, e_it);
                }
                ++e_it;
            }

#ifdef LOGGING_L0
            auto candidates_active =
                g.end_active_candidates() - g.begin_active_edges();

            printf(", candidates_communicated %d, candidates_active %ld\n",
                   candidates_communicated, candidates_active);
#endif

#ifdef MONITORING
            set_weights_time = MPI_Wtime();
#endif

            handle_matching_candidates(g, candidates, matchings, rng);

            // reset candidates
            std::fill(m_candidate_weight_of_node.begin(),
                      m_candidate_weight_of_node.end(), 0);

            ++round;
        }

        answer_until_finished(g, candidates, false);

        printf("[%d] Comm rounds: %d + waited %d\n", m_proc_id, m_comm_rounds,
               m_comm_wait_rounds);

        printf("[%d]: round %d, size %d, weight %.2f\n", m_proc_id, round,
               m_local_matching_size, m_local_matching_weight);

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Reduce(&m_local_matching_size, &total_global_matching_size, 1,
                   MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

        MPI_Reduce(&m_local_matching_weight, &total_global_matching_weight, 1,
                   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        MPI_Reduce(&round, &total_number_of_rounds, 1, MPI_INT, MPI_MAX, 0,
                   MPI_COMM_WORLD);
    }

  private:
    bool is_candidate(Graph &g, const_edge_iterator e_it) const {
        return is_candidate(g, e_it->n1, e_it->n2, e_it->weight);
    }

    bool is_candidate(Graph &g, NodeIDType local_n1, NodeIDType local_n2,
                      WeightType weight) const {

        // auto global_n1 = g.local_vertex_id_to_global_id(local_n1);
        // auto global_n2 = g.local_vertex_id_to_global_id(local_n2);

        return weight >= m_edge_weight_threshold *
                             m_candidate_weight_of_node[local_n1] &&
               weight >= m_edge_weight_threshold *
                             m_candidate_weight_of_node[local_n2];
    }

    void set_candidate_weight(Graph &g,
                              std::vector<NodeIDType> &nonzero_candidates) {
        edge_iterator e_it = g.begin_active_edges();
        while (e_it != g.end_active_edges()) {

            // auto global_n1 = g.local_vertex_id_to_global_id(e_it->n1);
            // auto global_n2 = g.local_vertex_id_to_global_id(e_it->n2);

            if (m_candidate_weight_of_node[e_it->n1] < e_it->weight) {
                if (m_candidate_weight_of_node[e_it->n1] == 0) {
                    nonzero_candidates.push_back(e_it->n1);
                }
                m_candidate_weight_of_node[e_it->n1] = e_it->weight;
            }
            if (m_candidate_weight_of_node[e_it->n2] < e_it->weight) {
                if ( //*g.is_local(e_it->n2) &&
                    m_candidate_weight_of_node[e_it->n2] == 0) {
                    nonzero_candidates.push_back(e_it->n2);
                }
                m_candidate_weight_of_node[e_it->n2] = e_it->weight;
            }

            // #ifdef LOGGING
            //             printf(
            //                 "[%d] Consider edge as candidate (v%d(%d),
            //                 v%d(%d), w%.2f)\n ", m_proc_id,
            //                 g.local_vertex_id_to_global_id(e_it->n1),
            //                 e_it->n1,
            //                 g.local_vertex_id_to_global_id(e_it->n2),
            //                 e_it->n2, e_it->weight);
            // #endif
            // #ifdef TESTING
            //             if (this->m_test_param_1 && is_candidate(e_it)) {
            //                 // swap edges to front
            //                 edge_is_candidate(g, e_it);
            //             }
            // #endif

            ++e_it;
        }
    }

    int sparse_global_max(Graph &g, std::vector<NodeIDType> &indices) {

        int local_count = indices.size();
        // Number of nonzero elements from each process
        std::vector<int> recv_counts(m_procs);
        // Gather Sizes
        MPI_Allgather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
                      MPI_COMM_WORLD);

        // Compute Displacements and Total Size
        std::vector<int> displs(m_procs, 0);
        unsigned int total_nonzero = 0;
        for (int i = 0; i < m_procs; ++i) {
            displs[i] = total_nonzero;
            total_nonzero += recv_counts[i];
        }

        if (total_nonzero > 3 * m_global_number_of_vertices) {

            std::vector<double> global_data(m_global_number_of_vertices, 0);

            for (unsigned int i = 0; i < indices.size(); ++i) {
                global_data[g.local_vertex_id_to_global_id(indices[i])] =
                    m_candidate_weight_of_node[indices[i]];
            }

            MPI_Allreduce(MPI_IN_PLACE, global_data.data(),
                          m_global_number_of_vertices, MPI_DOUBLE, MPI_MAX,
                          MPI_COMM_WORLD);

            // Update weights
            for (unsigned int i = 0; i < m_global_number_of_vertices; ++i) {
                auto local_node = m_global_to_local[i];

                if (local_node == (NodeIDType)-1) {
                    continue;
                }

                if (m_candidate_weight_of_node[local_node] < global_data[i]) {
                    m_candidate_weight_of_node[local_node] = global_data[i];
                }
            }
        } else {

            std::vector<double> values;
            values.reserve(indices.size());

            for (unsigned int i = 0; i < indices.size(); ++i) {
                values.push_back(m_candidate_weight_of_node[indices[i]]);
                indices[i] = g.local_vertex_id_to_global_id(indices[i]);
            }

            // Gather All (Global Index, Value) Pairs
            std::vector<int> global_indices(total_nonzero);
            std::vector<double> global_values(total_nonzero);

            MPI_Allgatherv(indices.data(), local_count, MPI_INT,
                           global_indices.data(), recv_counts.data(),
                           displs.data(), MPI_INT, MPI_COMM_WORLD);
            MPI_Allgatherv(values.data(), local_count, MPI_DOUBLE,
                           global_values.data(), recv_counts.data(),
                           displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

            // Update weights
            for (unsigned int i = 0; i < total_nonzero; ++i) {
                auto local_node = m_global_to_local[global_indices[i]];

                if (local_node == (NodeIDType)-1) {
                    continue;
                }

                if (m_candidate_weight_of_node[local_node] < global_values[i]) {
                    m_candidate_weight_of_node[local_node] = global_values[i];
                }
            }
        }

        return total_nonzero > 3 * m_global_number_of_vertices
                   ? m_global_number_of_vertices
                   : total_nonzero;
    }

    void handle_matching_candidates(Graph &g,
                                    std::vector<EdgeMinimal> &candidates,
                                    std::vector<std::vector<Edge>> &matchings,
                                    std::mt19937 &rng) {

        edge_iterator e_it = g.begin_active_edges();
        edge_iterator next;

        // int skipped_count = 0;
        while (e_it < g.end_active_candidates()) {
            next = e_it + 1;

            auto cc = find_free_color(g, e_it, matchings, rng);

            if (cc == color_set::not_colorable) {
                deactivate_candidate(g, e_it);
                // e_it = next;
                // ++skipped_count;
                continue;
            }

            EdgeMinimal edge;
            edge.first = g.local_vertex_id_to_global_id(e_it->n1);
            edge.second = g.local_vertex_id_to_global_id(e_it->n2);
            edge.color = cc;

            m_used_colors_not_own_v.clear();
            m_colored_edges.clear();
            ++m_comm_rounds;
            MPI_Allgather(&edge, 3, MPI_INT, candidates.data(), 3, MPI_INT,
                          MPI_COMM_WORLD);

            char valid_color = 1;
#ifdef LOGGING

#ifdef LOGGING_V
            m_print = false;
            for (int i = 0; i < m_procs; ++i) {
                m_print |= candidates[i].first == m_watched_global_id ||
                           candidates[i].second == m_watched_global_id;
            }
            if (m_print) {
#endif
                printf("[%d] candidates(cr%d):(", m_proc_id,
                       m_comm_rounds + m_comm_wait_rounds);
                for (int i = 0; i < m_procs; ++i) {
                    printf("[%d,%d,%d],", candidates[i].first,
                           candidates[i].second, candidates[i].color);
                }
                printf(")\n");

#ifdef LOGGING_V
            }
#endif
#endif
            for (int i = 0; i < m_procs; ++i) {

                if (i == m_proc_id) {
                    valid_color = true
#ifdef USE_HASH_MAP

                        !already_colored(edge.first, edge.second, false)
#endif
                        ;

// if (cc == color_set::not_colorable) {
// continue;
// }
#ifdef DEBUG

                    if (e_it->n1 !=
                            g.global_vertex_id_to_local_id(edge.first) ||
                        e_it->n2 !=
                            g.global_vertex_id_to_local_id(edge.second)) {
                        printf("[%d] Internal iterator error! (looked for "
                               "(%d,%d) v%d,v%d) found (%d,%d))\n",
                               m_proc_id,
                               g.global_vertex_id_to_local_id(edge.first),
                               g.global_vertex_id_to_local_id(edge.second),
                               edge.first, edge.second, e_it->n1, e_it->n2);

                        throw std::logic_error("Internal iterator error");
                    }
#endif

#ifdef LOGGING
#ifdef LOGGING_V
                    if (m_print)
#endif
                        printf("%d: (%d, %d), already_colored: %d, "
                               "avil n1: %d, avil n2: %d )\n",
                               m_proc_id, e_it->n1, e_it->n2, !valid_color,
                               m_color_available[e_it->n1][cc],
                               m_color_available[e_it->n2][cc]);
#endif
                    if (valid_color && m_color_available[e_it->n1][cc] &&
                        m_color_available[e_it->n2][cc]) {

#ifdef LOGGING
#ifdef LOGGING_V
                        if (m_print)
#endif
                            printf("[%d] Matched (v%d (%d), v%d (%d), "
                                   "color:%d)\n",
                                   m_proc_id, edge.first, e_it->n1, edge.second,
                                   e_it->n2, cc);
#endif
                        // No conflict for matched edge
                        e_it->color = cc;
                        add_to_matching(edge.first, edge.second, e_it->weight,
                                        matchings, cc);

                        m_color_available[e_it->n1].setOff(cc);
                        m_color_available[e_it->n2].setOff(cc);
#ifdef USE_HASH_MAP

                        m_colored_edges.insert(
                            edge.first < edge.second
                                ? std::make_pair(edge.first, edge.second)
                                : std::make_pair(edge.second, edge.first));
#endif
                    } else {
                        valid_color = 0;
                    }
                    continue;
                }
                process_cadidate(g, candidates[i], e_it);
            }

            if (valid_color || m_color_available[e_it->n1].none() ||
                m_color_available[e_it->n2].none()) {

                deactivate_candidate(g, e_it);

                continue;
            }

            e_it = next;
        }

        // printf("[%d] skipped %d candidates\n", m_proc_id, skipped_count);

        answer_until_finished(g, candidates, true);

        // e_it = g.begin_active_edges();
        while (e_it < g.end_active_edges()) {
            if (
                // e_it->color != color_set::npos ||
                m_color_available[e_it->n1].none() ||
                m_color_available[e_it->n2].none()) {
                deactivate_edge(g, e_it);
                continue;
            }

            ++e_it;
        }
    }

  private:
    void process_cadidate(Graph &g, EdgeMinimal candidate,
                          edge_iterator &e_it) {

        // if (candidate.color == color_set::not_colorable) {
        // invalid candidate
        // return;
        // }

        // candidate.first might be a local ghost
        auto local_n2 = candidate.first >= 0
                            ? m_global_to_local[candidate.first]
                            : candidate.first;

        // candidate.second might be a local node or a local
        // ghost
        auto local_n1 = candidate.second >= 0
                            ? m_global_to_local[candidate.second]
                            : candidate.second;
#ifdef LOGGING
#ifdef LOGGING_V
        if (m_print)
#endif
            printf("(v%d,v%d) local: (%d,%d)\n", candidate.second,
                   candidate.first, local_n1, local_n2);
#endif
        if (local_n1 == (NodeIDType)-1 && local_n2 == (NodeIDType)-1) {
            // not relevant -> but store used colors
            add_to_hashmap(candidate.second);
            add_to_hashmap(candidate.first);

            if (m_used_colors_not_own_v[candidate.second][candidate.color] &&
                m_used_colors_not_own_v[candidate.first][candidate.color]
#ifdef USE_HASH_MAP
                && !already_colored(candidate.first, candidate.second)
#endif
            ) {

                m_used_colors_not_own_v[candidate.second].setOff(
                    candidate.color);
                m_used_colors_not_own_v[candidate.first].setOff(
                    candidate.color);
            }
            return;
        }

        if (local_n1 != (NodeIDType)-1 && local_n2 != (NodeIDType)-1) {
            // safety check
            // if (!is_candidate(g, local_n1, local_n2, candidate.weight)) {
            //     printf("(v%d,v%d, w%d) local: (%d,%d)\n", candidate.second,
            //            candidate.first, candidate.weight, local_n1,
            //            local_n2);
            //     throw std::logic_error("No global candidate");
            // }
            // check for conflict
            if (m_color_available[local_n2][candidate.color] &&
                m_color_available[local_n1][candidate.color]) {
                // no conflict

#ifdef LOGGING
#ifdef LOGGING_V
                if (m_print)
#endif
                    printf("color avil c%d\n", candidate.color);
#endif
#ifdef USE_HASH_MAP

                // CAUTION:mutiple processes could
                // try to color the same edge in different
                // colors (!)
                if (already_colored(candidate.first, candidate.second)) {
                    return;
                }
#endif
#ifdef LOGGING
#ifdef LOGGING_V
                if (m_print)
#endif
                    printf("set color off\n");
#endif
                m_color_available[local_n1].setOff(candidate.color);
                m_color_available[local_n2].setOff(candidate.color);

                if (!g.is_local(local_n1)) {
                    // ghost-ghost edge -> to need to
                    // deactivate

#ifdef LOGGING
#ifdef LOGGING_V
                    if (m_print)
#endif
                        printf("ghost-ghost\n");
#endif
                    return;
                }
#ifdef USE_HASH_MAP

                // need to deactivate edge.
                // CAUTION: do not conflict with the current
                // edge pointer... (it only affected if cur
                // e_it is at the end of the active edges)

                auto received_pos =
                    m_cross_edge_id.at(std::make_pair(local_n1, local_n2));
                auto received_e_it = g.begin_active_edges() + received_pos;

                if (received_e_it->n1 != local_n1 ||
                    received_e_it->n2 != local_n2) {
                    printf(
                        "[%d] Internal hash map error! (looked for (% d, % d) "
                        "(v%d,v%d) at id%d, found (% d, % d))\n ",
                        m_proc_id, local_n1, local_n2, candidate.second,
                        candidate.first,
                        m_cross_edge_id.at(std::make_pair(local_n1, local_n2)),
                        received_e_it->n1, e_it->n2);

                    throw std::logic_error("Internal hash map error");
                }

                received_e_it->color = candidate.color;

                if (is_candidate(g, received_e_it)) {

                    // e_it is an active candidate!
                    // It might be swapped away from the end of active
                    // candidates

                    auto e_it_swapped_from_end_candidates =
                        e_it == (g.end_active_candidates() - 1);

                    auto e_it_equal_received_e_it = e_it == received_e_it;

                    deactivate_candidate(g, received_e_it);

                    if (e_it_equal_received_e_it) {

#ifdef LOGGING
#ifdef LOGGING_V
                        if (m_print)
#endif
                            printf("e_it_equal_received_e_it\n");
#endif

                        e_it = g.end_active_edges();
                    } else if (e_it_swapped_from_end_candidates) {
#ifdef LOGGING
#ifdef LOGGING_V
                        if (m_print)
#endif
                            printf("e_it_swapped_from_end_candidates\n");
#endif

                        // e_it was the last active edge and is now swapped
                        // to the previous position of the received_e_it we
                        // correct the postion of e_it
                        e_it = g.begin_active_edges() + received_pos;
                    }

                } else {
                    auto e_it_swapped_from_end_active =
                        e_it == (g.end_active_edges() - 1);

                    deactivate_edge(g, received_e_it);

                    if (e_it_swapped_from_end_active) {

#ifdef LOGGING
#ifdef LOGGING_V
                        if (m_print)
#endif
                            printf("e_it_swapped_from_end_active\n");
#endif

                        // e_it was the last active end and is now
                        // swapped to the previous position of the
                        // received_e_it we correct the postion of e_it
                        e_it = g.begin_active_edges() + received_pos;
                    }
                }
#endif
            }
            return;
        }

        if (local_n1 != (NodeIDType)-1) {
            add_to_hashmap(candidate.first);

            if (m_color_available[local_n1][candidate.color] &&
                m_used_colors_not_own_v[candidate.first][candidate.color]
#ifdef USE_HASH_MAP
                && !already_colored(candidate.first, candidate.second)
#endif
            ) {
                // no conflict
#ifdef LOGGING
#ifdef LOGGING_V
                if (m_print)
#endif
                    printf("local_n1: c%d\n", candidate.color);
#endif
                m_color_available[local_n1].setOff(candidate.color);
                m_used_colors_not_own_v[candidate.first].setOff(
                    candidate.color);
            }
            return;
        }

        if (local_n2 != (NodeIDType)-1) {
            add_to_hashmap(candidate.second);

            if (m_color_available[local_n2][candidate.color] &&
                m_used_colors_not_own_v[candidate.second][candidate.color]
#ifdef USE_HASH_MAP
                && !already_colored(candidate.first, candidate.second)
#endif
            ) {
                // no conflict
#ifdef LOGGING
#ifdef LOGGING_V
                if (m_print)
#endif
                    printf("local_n2: c%d\n", candidate.color);
#endif
                m_color_available[local_n2].setOff(candidate.color);
                m_used_colors_not_own_v[candidate.second].setOff(
                    candidate.color);
            }
            return;
        }
    }

    void answer_until_finished(Graph &g, std::vector<EdgeMinimal> &candidates,
                               bool handleUpdates = true) {

        int default_val = -1; // - (!handleUpdates);

        EdgeMinimal edge;
        edge.first = default_val;
        edge.second = default_val;

        edge_iterator e_it = g.end_active_edges();

        int total_nonzero = 0;
        do {
            if (!handleUpdates) {
                // reset candidates
                std::fill(m_candidate_weight_of_node.begin(),
                          m_candidate_weight_of_node.end(), 0);

                std::vector<NodeIDType> nonzero_candidates;
                total_nonzero = sparse_global_max(g, nonzero_candidates);

#ifdef LOGGING_L0

                printf("[%d] Communicated %d candidates weights\n", m_proc_id,
                       total_nonzero);
#endif

                // MPI_Allreduce(MPI_IN_PLACE,
                // m_candidate_weight_of_node.data(),
                //   m_candidate_weight_of_node.size(), MPI_DOUBLE,
                //   MPI_MAX, MPI_COMM_WORLD);
                //
                // printf("[%d] m_candidate_weight_of_node(cr%d):(", m_proc_id,
                //        m_comm_rounds + m_comm_wait_rounds);
                // for (int i = 0; i < m_candidate_weight_of_node.size(); ++i) {
                //     printf("v%d:w%.2f,", i, m_candidate_weight_of_node[i]);
                // }
                // printf(")\n");
            }
            char finished = 0;

            while (!finished) {
                finished = 1;

                m_used_colors_not_own_v.clear();
                m_colored_edges.clear();
                ++m_comm_wait_rounds;
                MPI_Allgather(&edge, 3, MPI_INT, candidates.data(), 3, MPI_INT,
                              MPI_COMM_WORLD);

#ifdef LOGGING
#ifdef LOGGING_V
                m_print = false;
                for (int i = 0; i < m_procs; ++i) {
                    m_print |= candidates[i].first == m_watched_global_id ||
                               candidates[i].second == m_watched_global_id;
                }
                if (m_print) {
#endif

                    printf("[%d] answer candidates(cr%d):(", m_proc_id,
                           m_comm_rounds + m_comm_wait_rounds);
                    for (int i = 0; i < m_procs; ++i) {
                        printf("[%d,%d,%d]", candidates[i].first,
                               candidates[i].second, candidates[i].color);
                    }
                    printf(")\n");
#ifdef LOGGING_V
                }
#endif
#endif

                for (int i = 0; i < m_procs; ++i) {
                    if (i == m_proc_id) {
                        continue;
                    }
                    if (handleUpdates) {
                        process_cadidate(g, candidates[i], e_it);
                    }
                    if (candidates[i].first <= default_val) {
                        continue;
                    }
                    finished = 0;
                }

#ifdef LOGGING
#ifdef LOGGING_V
                if (m_print)
#endif
                    printf("finished:%d\n", finished);
#endif
            }
#ifdef LOGGING
#ifdef LOGGING_V
            if (m_print)
#endif
                printf("any total_nonzero:%d, !handleUpdates %d\n",
                       total_nonzero, !handleUpdates);

#endif
        } while ((!handleUpdates) && total_nonzero);
    }

    void deactivate_candidate(Graph &g, edge_iterator &e_it) {
#ifdef USE_HASH_MAP
        auto id_before = g.get_edge_id(e_it);
#endif
        auto new_e_it = g.deactivate_candidate(e_it);
#ifdef USE_HASH_MAP

        if (e_it != g.end_active_edges() && !g.is_local(e_it->n2)) {

            m_cross_edge_id.find(std::make_pair(e_it->n1, e_it->n2))->second =
                id_before;
        }

        if (new_e_it != e_it && !g.is_local(new_e_it->n2)) {

            m_cross_edge_id.find(std::make_pair(new_e_it->n1, new_e_it->n2))
                ->second = new_e_it - g.begin_active_local_edges();
        }
#endif
        deactivate_edge(g, new_e_it);
    }
    void deactivate_edge(Graph &g, edge_iterator &e_it) {

#ifdef LOGGING_HASH_MAP
#ifdef LOGGING_V
        if (g.local_vertex_id_to_global_id(e_it->n1) == m_watched_global_id ||
            g.local_vertex_id_to_global_id(e_it->n2) == m_watched_global_id)
#endif
            printf("[%d] (c) Deactivated edge (v%d (%d), v%d (%d)) "
                   "[%d,%d,%d, %d]\n",
                   m_proc_id, g.local_vertex_id_to_global_id(e_it->n1),
                   e_it->n1, g.local_vertex_id_to_global_id(e_it->n2), e_it->n2,
                   e_it->color != color_set::npos,
                   m_color_available[e_it->n1].none(),
                   m_color_available[e_it->n2].none(), e_it->palette.none());
#endif
#ifdef USE_HASH_MAP
        auto id_before = g.get_edge_id(e_it);
#endif
        g.deactivate_edge(e_it);
#ifdef USE_HASH_MAP

        if (e_it != g.end_active_edges() && !g.is_local(e_it->n2)) {

            m_cross_edge_id.find(std::make_pair(e_it->n1, e_it->n2))->second =
                id_before;
        }
#endif
    }

    void edge_is_candidate(Graph &g, edge_iterator &e_it) {

        auto n1 = e_it->n1;
        auto n2 = e_it->n2;

        // #ifdef LOGGING
        //         printf("[%d] Edge id%d is candidate (v%d (%d), v%d (%d))\n",
        //         m_proc_id,
        //                id_before, g.local_vertex_id_to_global_id(e_it->n1),
        //                e_it->n1, g.local_vertex_id_to_global_id(e_it->n2),
        //                e_it->n2);
        // #endif

        g.edge_is_candidate(e_it);

#ifdef LOGGING_HASH_MAP
#ifdef LOGGING_V
        if (g.local_vertex_id_to_global_id(n1) == m_watched_global_id ||
            g.local_vertex_id_to_global_id(n2) == m_watched_global_id)
#endif
            printf("[%d] swapped (v%d (%d), v%d (%d)) - (v%d (%d), v%d (%d))\n",
                   m_proc_id, g.local_vertex_id_to_global_id(e_it->n1),
                   e_it->n1, g.local_vertex_id_to_global_id(e_it->n2), e_it->n2,
                   g.local_vertex_id_to_global_id(n1), n1,
                   g.local_vertex_id_to_global_id(n2), n2);
#endif

        if (n1 == e_it->n1 && n2 == e_it->n2) {
            return;
        }
#ifdef USE_HASH_MAP

        if (e_it != g.end_active_edges() && !g.is_local(e_it->n2)) {
            m_cross_edge_id.find(std::make_pair(e_it->n1, e_it->n2))->second =
                id_before;
#ifdef LOGGING_HASH_MAP
#ifdef LOGGING_V
            if (g.local_vertex_id_to_global_id(n1) == m_watched_global_id ||
                g.local_vertex_id_to_global_id(n2) == m_watched_global_id)
#endif
                printf("[%d] (v%d (%d), v%d (%d)) - id%d\n", m_proc_id,
                       g.local_vertex_id_to_global_id(e_it->n1), e_it->n1,
                       g.local_vertex_id_to_global_id(e_it->n2), e_it->n2,
                       id_before);
#endif
        }
        if (!g.is_local(n2)) {
            m_cross_edge_id.find(std::make_pair(n1, n2))->second =
                g.end_active_candidates() - g.begin_active_local_edges() - 1;

#ifdef LOGGING_HASH_MAP
#ifdef LOGGING_V
            if (g.local_vertex_id_to_global_id(n1) == m_watched_global_id ||
                g.local_vertex_id_to_global_id(n2) == m_watched_global_id)
#endif
                printf("[%d] (v%d (%d), v%d (%d)) - id%ld\n", m_proc_id,
                       g.local_vertex_id_to_global_id(n1), n1,
                       g.local_vertex_id_to_global_id(n2), n2,
                       g.end_active_candidates() -
                           g.begin_active_local_edges());
#endif
        }
#endif
    }

    uint_fast8_t find_free_color(Graph &g, edge_iterator &e_it,
                                 std::vector<std::vector<Edge>> &matchings,
                                 std::mt19937 &rng) {
        return find_free_color(g, e_it->n1, e_it->n2, e_it->weight, matchings,
                               rng);
    };

    uint_fast8_t find_free_color(Graph &g, NodeIDType n1, NodeIDType n2,
                                 WeightType weight,
                                 std::vector<std::vector<Edge>> &matchings,
                                 std::mt19937 &rng) {

        auto cc = color_set::common_colors(m_color_available[n1],
                                           m_color_available[n2])
                      .find_random(rng);

        if (cc == color_set::npos) {
            return color_set::not_colorable;
        }

        return cc;
    }
    void add_to_matching(NodeIDType n1, NodeIDType n2, WeightType weight,
                         std::vector<std::vector<Edge>> &matchings,
                         uint_fast8_t color) {

        // push back edge the since this fuction in only called by local
        // edges or one pf the processes of an cross edge
        matchings[color].push_back(Edge(n1, n2, weight));

        ++m_local_matching_size;
        m_local_matching_weight += weight;
        // printf("Add edge: %d, %d, %.2f\n", n1, n2, weight);
    }

    bool already_colored(NodeIDType first, NodeIDType second,
                         bool insert = true) {
        if (m_colored_edges.contains(first < second
                                         ? std::make_pair(first, second)
                                         : std::make_pair(second, first))) {

#ifdef LOGGING
#ifdef LOGGING_V
            if (m_print)
#endif
                printf("edge already colored\n");
#endif
            return true;
        }
        if (insert) {
            m_colored_edges.insert(first < second
                                       ? std::make_pair(first, second)
                                       : std::make_pair(second, first));
        }
        return false;
    }

    void add_to_hashmap(NodeIDType node) {
        if (!m_used_colors_not_own_v.contains(node)) {
            m_used_colors_not_own_v.insert(
                std::make_pair(node, color_set(m_k)));
        }
    }

#ifdef MONITORING
    int m_asked_candidates_per_round = 0;
#endif

#ifdef LOGGING_V
    NodeIDType m_watched_global_id = 0;
    int m_print = false;
#endif

    int m_procs, m_proc_id;

    int m_blocksize;

    unsigned int m_global_number_of_vertices;
    int m_k;
    int m_comm_rounds, m_comm_wait_rounds;
    unsigned int m_local_matching_size = 0;
    double m_local_matching_weight = 0.0;

    double m_edge_weight_threshold;

    std::vector<NodeIDType> m_global_to_local;

    absl::flat_hash_map<NodeIDType, color_set> m_used_colors_not_own_v;
    absl::flat_hash_set<std::pair<NodeIDType, NodeIDType>> m_colored_edges;

    std::vector<color_set> m_color_available;
    std::vector<WeightType> m_candidate_weight_of_node;

    absl::flat_hash_map<std::pair<NodeIDType, NodeIDType>, EdgeIDType>
        m_cross_edge_id;

    std::vector<unsigned int> m_cross_deg;

    std::vector<std::vector<int>> m_procs_of_cross_candidates;
};

} // namespace parkec

#endif
