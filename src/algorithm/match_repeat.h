/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *
 * The code for a single matching is based on LocalMaxMatching (MIT Licence,
 * Marcel Birn) https://github.com/LocalMaxMatching/LocalMaxMatching
 *
 * M. Birn, V. Osipov, P. Sanders, C. Schulz, and N. Sitchinava,
 * “Efficient Parallel and External Matching,” in
 * Euro-Par 2013 Parallel Processing - 19th International Conference,
 * Aachen, Germany, August 26-30, 2013. Proceedings,
 * in Lecture Notes in Computer Science, vol. 8097. Springer, 2013, pp. 659–670.
 * doi: 10.1007/978-3-642-40047-6_66.
 *
 *****************************************************************************/

#ifndef GREEDY_MATCHING_REPEAT_H
#define GREEDY_MATCHING_REPEAT_H

#include <iostream>
#include <limits>
#include <list>
#include <utility>
#include <vector>

#include <mpi.h>

#include <util/hash_functions.h>
#include <util/math_funcs.h>
#include <util/mpi_tags.h>

#include "datastructure/parallel_edge_graph.h"
#include "k_edge_coloring.h"

// #define LOGGING
// #define CHECK

namespace parkec {

template <class GraphBase>
class GreedyMatchingRepeat : public KEdgeColoring<GraphBase> {
    typedef typename parkec::ParallelEdgeGraph<> Graph;
    typedef typename Graph::WeightType WeightType;
    typedef typename Graph::NodeIDType NodeIDType;
    typedef typename Graph::EdgeIDType EdgeIDType;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::edge_iterator edge_iterator;
    typedef typename Graph::const_edge_iterator const_edge_iterator;

    struct Candidate {
        WeightType weight;
        NodeIDType partner;

        Candidate()
            : weight(min_val<WeightType>()),
              partner(
                  0) /// \todo this only works as long min_val<WeightType>() is
                     /// smaller then any weight that might occur in graphs
        {}

        Candidate(const WeightType weight, const NodeIDType partner)
            : weight(weight), partner(partner) {}

        void set(const WeightType weight, const NodeIDType partner) {
            this->weight = weight;
            this->partner = partner;
        }
    };

    struct Message {
        NodeIDType src;
        NodeIDType candidate;

        Message() : src(-1), candidate(-1) {}

        Message(const NodeIDType src, const NodeIDType candidate)
            : src(src), candidate(candidate) {}
    };

  public:
    GreedyMatchingRepeat() {
        // init MSG sizes
        MPI_Type_contiguous(sizeof(Message), MPI_BYTE, &CANDIDATE_MESSAGE_TYPE);
        MPI_Type_commit(&CANDIDATE_MESSAGE_TYPE);

        MPI_Type_contiguous(sizeof(NodeIDType), MPI_BYTE,
                            &MATCHED_MESSAGE_TYPE);
        MPI_Type_commit(&MATCHED_MESSAGE_TYPE);
    }

    bool better_than_old_partner(const Graph &g, const NodeIDType src_id,
                                 const WeightType old_weight,
                                 const NodeIDType old_partner_id,
                                 const WeightType new_weight,
                                 const NodeIDType new_partner_id) {
        /// \todo think about the type of the hash, also the return type
        unsigned int old_hash = hash(old_partner_id ^ src_id);
        unsigned int new_hash = hash(new_partner_id ^ src_id);

        return (new_weight > old_weight) ||
               (new_weight == old_weight && new_hash < old_hash) ||
               (new_weight == old_weight && new_hash == old_hash &&
                (new_partner_id ^ src_id) < (old_partner_id ^ src_id));
    }

    void compute_k_edge_coloring(std::vector<std::vector<Edge>> &matchings,
                                 int &total_global_matching_size,
                                 int &total_number_of_rounds,
                                 double &total_global_matching_weight,
                                 GraphBase &g_ref, unsigned int k) {

        Graph &g = static_cast<Graph &>(g_ref);
        int proc_id;
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

        bool is_root = proc_id == 0;
        for (unsigned int i = 0; i < k; ++i) {
#ifdef LOGGING
            printf("[%d] Start color %d, active edges left: %d [local: %ld, "
                   "cross: %ld]\n",
                   proc_id, i, g.has_active_edges(),
                   g.end_active_local_edges() - g.begin_active_local_edges(),
                   g.end_active_cross_edges() - g.begin_active_cross_edges());
#endif
            MPI_Barrier(MPI_COMM_WORLD);
            std::vector<Edge> matching;

            std::vector<bool> matched(g.num_all_vertices(), false);

            unsigned int depth = 0;

            compute_weighted_matching(g, matching, depth, matched, i);

            MPI_Barrier(MPI_COMM_WORLD);

            matchings.push_back(matching);

#ifdef CHECK
            bool maximal_matching = g.is_valid_matching(matching, i);
#endif
            unsigned int global_matching_size;
            unsigned int local_matching_size = matching.size();

            MPI_Reduce(&local_matching_size, &global_matching_size, 1,
                       MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

            double global_weight = get_weight(g, matching);

            unsigned int global_depth;
            MPI_Reduce(&depth, &global_depth, 1, MPI_UNSIGNED, MPI_MAX, 0,
                       MPI_COMM_WORLD);

            if (is_root) {
                printf("Color: %d/%d\n", i + 1, k);
                printf("Rounds: %d\n", global_depth);
                printf("Size of Matching: %d\n", global_matching_size);
                printf("Weight of Matching: %.2f\n", global_weight);
#ifdef CHECK
                printf("Is maximal Matching: %d\n", maximal_matching);
#endif
            }

            total_global_matching_size += global_matching_size;
            total_global_matching_weight += global_weight;
            total_number_of_rounds += global_depth;

            if (i == k - 1) {
                // no need to reactivate the graph in the last iteration
                break;
            }

            g.activate_edges();

            // deactivate colored edges
            edge_iterator e_it = g.begin_active_cross_edges();
            while (e_it != g.end_active_cross_edges()) {
                if (e_it->color != Edge::no_color) {
                    g.deactivate_cross_edge(e_it);
                } else {
                    ++e_it;
                }
            }

            e_it = g.begin_active_local_edges();
            while (e_it < g.end_active_local_edges()) {
                if (e_it->color != Edge::no_color) {
                    g.deactivate_local_edge(e_it);
                } else {
                    ++e_it;
                }
            }
        }
    }

    /*!
     * \brief Computes a maximal weighted matching of the given graph.
     *
     * Runtime: ??
     *
     * @param g			The Graph
     * @param matching	Return list of the resulting matching
     */
    void compute_weighted_matching(Graph &g, std::vector<Edge> &matching,
                                   unsigned int &round,
                                   std::vector<bool> &matched,
                                   unsigned int color) {

        int procs, proc_id;
        MPI_Comm_size(MPI_COMM_WORLD, &procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

        // buffers for messages that are used to send candidate requests of
        // ghostnodes
        std::vector<std::pair<std::vector<Message>, MPI_Request>>
            candidate_messages_of_proc;
        // buffers for messages that are used to inform other processes about
        // matched ghostnodes
        std::vector<std::pair<std::vector<NodeIDType>, MPI_Request>>
            matched_messages_of_proc;

        candidate_messages_of_proc.resize(
            procs, make_pair(std::vector<Message>(), MPI_Request()));
        matched_messages_of_proc.resize(
            procs, make_pair(std::vector<NodeIDType>(), MPI_Request()));

        // stores for each vertex the incident edge id, of the edge with the
        // largest weight
        std::vector<Candidate> candidate_of_node(
            g.num_all_vertices(),
            Candidate()); // initialize with dummy candidate

        round = 0;

        while (g.has_active_edges() /*&& round<4*/) {
#ifdef LOGGING
            printf("[%d] Start round %d\n", proc_id, round);
#endif
            set_local_maximal_candidate_of_nodes(g, candidate_of_node, matched);

            set_matched_local_nodes_and_add_edge_to_matching(
                g, candidate_of_node, matched, matching);
            set_matched_ghost_nodes_and_add_edge_to_matching(
                g, candidate_of_node, candidate_messages_of_proc,
                matched_messages_of_proc, color * 1000 + round, matched,
                matching);

            deactivate_edges_incident_to_matched_nodes(g, matched,
                                                       candidate_of_node);

#ifdef LOGGING
            printf("[%d] End round %d, has active edges left: %d [local: %ld, "
                   "cross: %ld]\n",
                   proc_id, round, g.has_active_edges(),
                   g.end_active_local_edges() - g.begin_active_local_edges(),
                   g.end_active_cross_edges() - g.begin_active_cross_edges());
#endif
            round++;
        }
    }

    // only process 0 gets the correct global result
    double get_weight(const Graph &g, std::vector<Edge> &matching) {
        double local_result = 0.;
        double result = 0.;

        for (typename std::vector<Edge>::iterator it = matching.begin();
             it != matching.end(); it++) {
            local_result += it->weight;
        }

        MPI_Reduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);

        return result;
    }

  private:
    void set_local_maximal_candidate_of_nodes(
        Graph &g, std::vector<Candidate> &candidate_of_node,
        const std::vector<bool> &matched) {
        // iterate over each active local edge, to set local max edges
        for (edge_iterator e_it = g.begin_active_local_edges(),
                           it_end = g.end_active_local_edges();
             e_it < it_end; e_it++) {
            const NodeIDType adjusted_e_it_n1 =
                g.local_vertex_id_to_global_id(e_it->n1);
            const NodeIDType adjusted_e_it_n2 =
                g.local_vertex_id_to_global_id(e_it->n2);

            /// \todo maybe store the global vertex IDs in the candidate
            /// objects?
            // check first endpoint of current edge
            if (better_than_old_partner(
                    g, adjusted_e_it_n1, candidate_of_node[e_it->n1].weight,
                    g.local_vertex_id_to_global_id(
                        candidate_of_node[e_it->n1].partner),
                    e_it->weight, adjusted_e_it_n2)) {
                candidate_of_node[e_it->n1].set(e_it->weight, e_it->n2);
            }

            // check second endpoint of current edge
            if (better_than_old_partner(
                    g, adjusted_e_it_n2, candidate_of_node[e_it->n2].weight,
                    g.local_vertex_id_to_global_id(
                        candidate_of_node[e_it->n2].partner),
                    e_it->weight, adjusted_e_it_n1)) {
                candidate_of_node[e_it->n2].set(e_it->weight, e_it->n1);
            }
        }

        // iterate over each active local ghost edge, to set local max edges
        for (edge_iterator e_it = g.begin_active_cross_edges(),
                           it_end = g.end_active_cross_edges();
             e_it < it_end;
             e_it++) { // we have to set the candidates of both vertices -> we
                       // need this information to send the correct information

            const NodeIDType adjusted_e_it_n1 =
                g.local_vertex_id_to_global_id(e_it->n1);
            const NodeIDType adjusted_e_it_n2 =
                g.local_vertex_id_to_global_id(e_it->n2);

            // check first endpoint of current edge
            if (better_than_old_partner(
                    g, adjusted_e_it_n1, candidate_of_node[e_it->n1].weight,
                    g.local_vertex_id_to_global_id(
                        candidate_of_node[e_it->n1].partner),
                    e_it->weight, adjusted_e_it_n2)) {
                candidate_of_node[e_it->n1].set(e_it->weight, e_it->n2);
            }

            // check second endpoint of current edge
            if (better_than_old_partner(
                    g, adjusted_e_it_n2, candidate_of_node[e_it->n2].weight,
                    g.local_vertex_id_to_global_id(
                        candidate_of_node[e_it->n2].partner),
                    e_it->weight, adjusted_e_it_n1)) {
                candidate_of_node[e_it->n2].set(e_it->weight, e_it->n1);
            }
        }
    }

    void set_matched_local_nodes_and_add_edge_to_matching(
        Graph &g, std::vector<Candidate> &candidate_of_node,
        std::vector<bool> &matched, std::vector<Edge> &matching) {
        // add local matching edges to matchings list and set them inactive,
        // also set the incident nodes inactive
        //
        // We have to check for local dominant edges before we handle dominant
        // cross-edges. This makes it easy to tell partner processes about
        // locally matched nodes in the same round. Thus after each round every
        // process knows the correct amount of active cross edges.
        for (edge_iterator e_it = g.begin_active_local_edges(),
                           it_end = g.end_active_local_edges();
             e_it < it_end; e_it++) {
            // check first endpoint of current edge
            if (e_it->n2 == candidate_of_node[e_it->n1].partner &&
                e_it->n1 == candidate_of_node[e_it->n2].partner &&
                e_it->weight == candidate_of_node[e_it->n1].weight &&
                !matched[e_it->n1] &&
                !matched[e_it->n2]) { // have to check for matched, because of
                                      // multi-graphs

#ifdef LOGGING
                printf("Matched (local) (v%d (%d), v%d (%d))\n",
                       g.local_vertex_id_to_global_id(e_it->n1), e_it->n1,
                       g.local_vertex_id_to_global_id(e_it->n2), e_it->n2);
#endif

                matching.push_back(Edge(
                    g.local_vertex_id_to_global_id(e_it->n1),
                    g.local_vertex_id_to_global_id(e_it->n2), e_it->weight));

                matched[e_it->n1] = true;
                matched[e_it->n2] = true;
                e_it->color = 1;
            }
        }
    }

    void set_matched_ghost_nodes_and_add_edge_to_matching(
        Graph &g, std::vector<Candidate> &candidate_of_node,
        std::vector<std::pair<std::vector<Message>, MPI_Request>>
            &candidate_messages_of_proc,
        std::vector<std::pair<std::vector<NodeIDType>, MPI_Request>>
            &matched_messages_of_proc,
        unsigned int round, std::vector<bool> &matched,
        std::vector<Edge> &matching) {

        int procs, proc_id;
        MPI_Comm_size(MPI_COMM_WORLD, &procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

        get_ghost_candidates(g, candidate_of_node, matched,
                             candidate_messages_of_proc);

        // count the messages sent ou, because the same number has to be
        // received
        unsigned int count_msg_send_receive;

        send_messages(g, round, count_msg_send_receive,
                      candidate_messages_of_proc);

        receive_ghost_candidates(g, round, count_msg_send_receive,
                                 candidate_of_node, matched);

        clear_messages(candidate_messages_of_proc);

        add_matched_cross_edges_to_matching_and_prepare_messages_for_ghost_nodes_adjacent_to_matched_local_nodes(
            g, candidate_of_node, matched, matching, matched_messages_of_proc);

        send_messages(g, round, count_msg_send_receive,
                      matched_messages_of_proc);

        receive_information_about_matched_ghost_nodes(g, count_msg_send_receive,
                                                      round, matched);

        clear_messages(matched_messages_of_proc);
    }

    void get_ghost_candidates(
        const Graph &g, const std::vector<Candidate> &candidate_of_node,
        const std::vector<bool> &matched,
        std::vector<std::pair<std::vector<Message>, MPI_Request>>
            &messages_of_proc) {
        int procs, proc_id;
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
        MPI_Comm_size(MPI_COMM_WORLD, &procs);

        // iterate ghost edges and send appropriate messages
        for (const_edge_iterator e_it = g.begin_active_cross_edges(),
                                 it_end = g.end_active_cross_edges();
             e_it < it_end; e_it++) {
            // n1 is always local -> n2 is ghost
            if (!matched[e_it->n1] &&
                e_it->n2 == candidate_of_node[e_it->n1].partner &&
                e_it->n1 == candidate_of_node[e_it->n2].partner &&
                e_it->weight == candidate_of_node[e_it->n1]
                                    .weight) /// \todo weight comparison only
                                             /// necessary for multi-graphs
            { // we might send unnecessary messages if there are multi-edges
                // with the same weight
                messages_of_proc[g.get_proc_of_ghost_vertex(e_it->n2)]
                    .first.push_back(
                        Message(g.local_vertex_id_to_global_id(e_it->n1),
                                g.local_vertex_id_to_global_id(e_it->n2)));
            }
        }
    }

    void send_messages(const Graph &g, unsigned int &round,
                       unsigned int &send_msg_count,
                       std::vector<std::pair<std::vector<Message>, MPI_Request>>
                           &messages_of_proc) {
        int procs, proc_id;
        MPI_Comm_size(MPI_COMM_WORLD, &procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

        send_msg_count = 0;

        int tag = MPI_MSG_TAG::get_candidates_tag(round);

        // for each active partner: send messages;
        /// \todo think about using a list with active partners
        for (int p = 0; p < procs; p++) {
            if (!messages_of_proc[p]
                     .first
                     .empty()) { // check for messages_of_proc[p].size() ==
                                 // g.get_active_cross_edges_of_proc(p)??
                MPI_Isend(&messages_of_proc[p].first[0],
                          messages_of_proc[p].first.size(),
                          CANDIDATE_MESSAGE_TYPE, p, tag, MPI_COMM_WORLD,
                          &messages_of_proc[p].second);

#ifdef LOGGING
                printf("[%d] Send msg (size %ld) with tag %d to p%d\n", proc_id,
                       messages_of_proc[p].first.size(), tag, p);
#endif

                send_msg_count++;
            } else if (g.is_active_partner(p)) {
                // we have to send empty message to active partner, because they
                // expect a message from each active partner even if this
                // message is empty!!
                MPI_Isend(0, 0, CANDIDATE_MESSAGE_TYPE, p, tag, MPI_COMM_WORLD,
                          &messages_of_proc[p].second);
#ifdef LOGGING
                printf("[%d] Send empty msg with tag %d to p%d\n", proc_id, tag,
                       p);
#endif
                send_msg_count++;
            }
        }
    }

    void
    send_messages(const Graph &g, unsigned int &round,
                  unsigned int &send_msg_count,
                  std::vector<std::pair<std::vector<NodeIDType>, MPI_Request>>
                      &messages_of_proc) {
        int procs, proc_id;
        MPI_Comm_size(MPI_COMM_WORLD, &procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

        send_msg_count = 0;

        int tag = MPI_MSG_TAG::get_ghost_vertices_tag(round);

        // for each active partner: send messages;
        /// \todo think about using a list with active partners
        for (int p = 0; p < procs; p++) {
            if (!messages_of_proc[p]
                     .first
                     .empty()) { // check for messages_of_proc[p].size() ==
                                 // g.get_active_cross_edges_of_proc(p)??
                MPI_Isend(&messages_of_proc[p].first[0],
                          messages_of_proc[p].first.size(),
                          MATCHED_MESSAGE_TYPE, p, tag, MPI_COMM_WORLD,
                          &messages_of_proc[p].second);
                send_msg_count++;
#ifdef LOGGING
                printf("[%d] Send msg (size %ld) with tag %d to p%d\n", proc_id,
                       messages_of_proc[p].first.size(), tag, p);
#endif
            } else if (g.is_active_partner(p)) {
                // we have to send empty message to active partner, because they
                // expect a message from each active partner even if this
                // message is empty!!
                MPI_Isend(0, 0, MATCHED_MESSAGE_TYPE, p, tag, MPI_COMM_WORLD,
                          &messages_of_proc[p].second);
                send_msg_count++;
#ifdef LOGGING
                printf("[%d] Send empty msg with tag %d to p%d\n", proc_id, tag,
                       p);
#endif
            }
        }
    }

    void
    receive_ghost_candidates(const Graph &g, unsigned int &round,
                             unsigned int &msgs_to_receive,
                             const std::vector<Candidate> &candidate_of_node,
                             std::vector<bool> &matched) {
        int procs, proc_id;
        MPI_Comm_size(MPI_COMM_WORLD, &procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

        MPI_Status msg_status;

        int tag = MPI_MSG_TAG::get_candidates_tag(round);

        // receive messages
        //		for(int p=0; p<procs; p++)
        //		{
        //			if(!g.is_active_partner(p))
        //			{
        //				continue;
        //			}
        while (msgs_to_receive > 0) {
            MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &msg_status);

            int msg_count;
            MPI_Get_count(&msg_status, CANDIDATE_MESSAGE_TYPE, &msg_count);

            int src = msg_status.MPI_SOURCE;

            Message *messages = new Message[msg_count];

            MPI_Recv(messages, msg_count, CANDIDATE_MESSAGE_TYPE, src, tag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#ifdef LOGGING
            printf("[%d] Receive msg (size %d) with tag %d from p%d\n", proc_id,
                   msg_count, tag, src);
#endif
            msgs_to_receive--;

            for (int i = 0; i < msg_count; i++) {
                Message msg = messages[i];
                NodeIDType local_ghost_id =
                    g.global_vertex_id_to_local_id_of_ghost_vertex(
                        msg.src); /// \todo not nice, get_local_ghost_id is
                                  /// expensive
                NodeIDType local_id =
                    g.global_vertex_id_to_local_id(msg.candidate);

                if (candidate_of_node[local_id].partner == local_ghost_id) {
                    /// \todo I don't think that we have to adjust
                    /// candidate_of_node, because the candidates have been set
                    /// in an earlier step, otherwise we wouldn't have received
                    /// this message or the ghostnode wouldn't be the partner
                    matched[local_id] = true;
                    matched[local_ghost_id] = true;

#ifdef LOGGING
                    printf("Matched (cross) (v%d (%d), v%d (%d))\n", msg.src,
                           local_id, msg.candidate, local_ghost_id);
#endif
                }
            }

            delete[] messages;
        }
    }

    // forces some kind of synchronization between communicating processes
    /*!
     *
     * @param messages_of_proc
     */
    template <typename MTYPE>
    void clear_messages(std::vector<std::pair<std::vector<MTYPE>, MPI_Request>>
                            &messages_of_proc) {
        int procs;
        MPI_Comm_size(MPI_COMM_WORLD, &procs);

        // clear all message vectors
        // we might have send a few more messages but all of them were empty,
        // thus we don't have to clear messages_of_proc[p].first, because it's
        // already empty.
        for (int p = 0; p < procs; p++) {
            if (!messages_of_proc[p].first.empty()) {
                MPI_Wait(&messages_of_proc[p].second, MPI_STATUS_IGNORE);
                messages_of_proc[p]
                    .first
                    .clear(); // probably doesn't free any allocated memory, but
                              // in this case that's not too bad. I'm not
                              // totally sure about the definition of clear
            }
        }
    }

    void deactivate_local_edges_incident_to_matched_nodes(
        Graph &g, const std::vector<bool> &matched,
        std::vector<Candidate> &candidate_of_node) {
        Candidate base_candidate = Candidate();
        // deactivate local edges that are incident to matched nodes
        edge_iterator e_it = g.begin_active_local_edges();
        while (e_it < g.end_active_local_edges()) {
            // check first endpoint of current edge
            if (matched[e_it->n1] || matched[e_it->n2]) {
                g.deactivate_local_edge(e_it);
                continue; // we just set an edge to inactive, thus the current
                          // iterator position references a new active edge
            }

            // reset the candidate for the nodes that aren't matched
            // we might miss a few nodes, but those are unmatched nodes,
            // that are no longer incident to an active edge
            candidate_of_node[e_it->n1] = base_candidate;
            candidate_of_node[e_it->n2] = base_candidate;

            e_it++;
        }
    }

    void
    add_matched_cross_edges_to_matching_and_prepare_messages_for_ghost_nodes_adjacent_to_matched_local_nodes(
        Graph &g, const std::vector<Candidate> &candidate_of_node,
        const std::vector<bool> &matched, std::vector<Edge> &matching,
        std::vector<std::pair<std::vector<NodeIDType>, MPI_Request>>
            &messages_of_proc) {
        for (edge_iterator e_it = g.begin_active_cross_edges(),
                           it_end = g.end_active_cross_edges();
             e_it != it_end; e_it++) {
            // check first endpoint of current edge
            if (matched[e_it->n1] || matched[e_it->n2]) {
                /// \todo because of this condition no multi-graphs are
                /// supported which have multi-edges with the same weight
                if (matched[e_it->n1] && matched[e_it->n2] &&
                    candidate_of_node[e_it->n1].partner == e_it->n2 &&
                    candidate_of_node[e_it->n2].partner ==
                        e_it->n1) { /// \todo only add matched cross edges on
                                    /// one process, maybe on the one with the
                                    /// smaller proc-ID (or use the node ID to
                                    /// decide)

                    if ((g.is_local(e_it->n1) &&
                         g.local_vertex_id_to_global_id(e_it->n1) <
                             g.local_vertex_id_to_global_id(e_it->n2)) ||
                        (g.is_local(e_it->n2) &&
                         g.local_vertex_id_to_global_id(e_it->n2) <
                             g.local_vertex_id_to_global_id(
                                 e_it->n1))) { // only add the edge if the
                                               // global id of the local
                                               // endpoint is smaller than the
                                               // global id of the ghost
                                               // endpoint
                        matching.push_back(
                            Edge(g.local_vertex_id_to_global_id(e_it->n1),
                                 g.local_vertex_id_to_global_id(e_it->n2),
                                 e_it->weight));
                    }

                    e_it->color = 1;

                } else if (matched[e_it->n1] && g.is_local(e_it->n1)) {
                    messages_of_proc[g.get_proc_of_ghost_vertex(e_it->n2)]
                        .first.push_back(
                            g.local_vertex_id_to_global_id(e_it->n1));
                } else if (matched[e_it->n2] && g.is_local(e_it->n2)) {
                    messages_of_proc[g.get_proc_of_ghost_vertex(e_it->n1)]
                        .first.push_back(
                            g.local_vertex_id_to_global_id(e_it->n2));
                }
            }
        }
    }

    void receive_information_about_matched_ghost_nodes(
        const Graph &g, unsigned int &msgs_to_receive, unsigned int &round,
        std::vector<bool> &matched) {
        int procs, proc_id;
        MPI_Comm_size(MPI_COMM_WORLD, &procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

        MPI_Status msg_status;

        int tag = MPI_MSG_TAG::get_ghost_vertices_tag(round);

        while (msgs_to_receive > 0) {
            MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &msg_status);

            int msg_count;
            MPI_Get_count(&msg_status, MATCHED_MESSAGE_TYPE, &msg_count);

            NodeIDType *messages = new NodeIDType[msg_count];

            int src = msg_status.MPI_SOURCE;

            MPI_Recv(messages, msg_count, MATCHED_MESSAGE_TYPE, src, tag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#ifdef LOGGING
            printf("[%d] Receive msg (size %d) with tag %d from p%d\n", proc_id,
                   msg_count, tag, src);
#endif
            msgs_to_receive--;

            for (int i = 0; i < msg_count; i++) {
                const NodeIDType global_ghost_id = messages[i];
                const NodeIDType local_ghost_id =
                    g.global_vertex_id_to_local_id_of_ghost_vertex(
                        global_ghost_id); /// \todo not nice, get_local_ghost_id
                                          /// is expensive

                matched[local_ghost_id] = true;
            }

            delete[] messages;
        }
    }

    void deactivate_cross_edges_incident_to_matched_nodes(
        Graph &g, const std::vector<bool> &matched,
        std::vector<Candidate> &candidate_of_node) {
        Candidate base_candidate = Candidate();

        // deactivate all remaining cross edges that are incident to a matched
        // ghost-vertex
        edge_iterator e_it = g.begin_active_cross_edges();
        while (e_it != g.end_active_cross_edges()) {
            // check first endpoint of current edge
            if (matched[e_it->n1] || matched[e_it->n2]) {
                g.deactivate_cross_edge(e_it);

                continue; // we just set an edge to inactive, thus the current
                          // iterator position references a new active edge
            }

            // reset the candidate for the nodes that aren't matched
            // we might miss a few nodes, but those are unmatched nodes,
            // that are no longer incident to an active edge
            candidate_of_node[e_it->n1] = base_candidate;
            candidate_of_node[e_it->n2] = base_candidate;

            e_it++;
        }
    }

    inline void deactivate_edges_incident_to_matched_nodes(
        Graph &g, const std::vector<bool> &matched,
        std::vector<Candidate> &candidate_of_node) {
        deactivate_local_edges_incident_to_matched_nodes(g, matched,
                                                         candidate_of_node);
        deactivate_cross_edges_incident_to_matched_nodes(g, matched,
                                                         candidate_of_node);
    }

  private:
    MPI_Datatype CANDIDATE_MESSAGE_TYPE;
    MPI_Datatype MATCHED_MESSAGE_TYPE;
};

} // namespace parkec

#endif
