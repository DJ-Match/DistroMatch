/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef GREEDY_K_MATCHING_H
#define GREEDY_K_MATCHING_H

#include <iostream>
#include <limits>
#include <list>
#include <numeric>
#include <utility>
#include <vector>

#include <mpi.h>

#include <util/color_set.h>
#include <util/hash_functions.h>
#include <util/math_funcs.h>
#include <util/mpi_tags.h>

#include "datastructure/dis_adj_graph.h"

#include "k_edge_coloring.h"

// #define LOGGING
// #define CHECK

#ifdef LOGGING
#include <stdio.h>

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_RESET "\x1b[0m"

#endif
namespace parkec {

template <class GraphBase>
class GreedyKMatching : public KEdgeColoring<GraphBase> {
    typedef typename parkec::DisAdjGraph<> Graph;

    typedef typename Graph::WeightType WeightType;
    typedef typename Graph::NodeIDType NodeIDType;
    typedef typename Graph::EdgeIDType EdgeIDType;
    typedef typename Graph::Edge Edge;
    typedef typename Graph::EdgeE EdgeE;
    typedef typename Graph::edge_iterator edge_iterator;
    typedef typename Graph::const_edge_iterator const_edge_iterator;

    struct Message {
        NodeIDType src;
        NodeIDType candidate;
        color_set::bit_type colors;

        Message() : src(-1), candidate(-1) {}

        Message(const NodeIDType src, const NodeIDType candidate,
                color_set::bit_type colors)
            : src(src), candidate(candidate), colors(colors) {}
    };

  public:
    GreedyKMatching() {

        // init MSG sizes
        MPI_Type_contiguous(sizeof(Message), MPI_BYTE, &CANDIDATE_MESSAGE_TYPE);
        MPI_Type_commit(&CANDIDATE_MESSAGE_TYPE);

        MPI_Comm_size(MPI_COMM_WORLD, &m_procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &m_proc_id);
    }
    void set_color() { m_set_color = true; }

    void reset() {

        for (NodeIDType i = 0; i < m_color_available.size(); ++i) {
            m_color_available[i].set();
            m_candidate_of_node[i] = EdgeE();
            m_curr_candidate_id[i] = 0;
        }

        std::fill(m_ghost_is_active.begin(), m_ghost_is_active.end(), 1);

        m_local_matching_size = 0;
        m_local_matching_weight = 0.0;
    }

    void compute_k_edge_coloring(std::vector<std::vector<Edge>> &matchings,
                                 int &total_global_matching_size,
                                 int &total_number_of_rounds,
                                 double &total_global_matching_weight,
                                 GraphBase &g_ref, unsigned int k) {

        matchings.resize(k, std::vector<Edge>(0));
        compute_k_matchings(matchings, total_global_matching_size,
                            total_number_of_rounds,
                            total_global_matching_weight, g_ref, k);
    }

    void compute_k_matchings(std::vector<std::vector<Edge>> &matchings,
                             int &total_global_matching_size,
                             int &total_number_of_rounds,
                             double &total_global_matching_weight,
                             GraphBase &g_ref, unsigned int k) {

        Graph &g = static_cast<Graph &>(g_ref);

        m_color_available.resize(g.num_local_vertices(), color_set(k));
        m_active_deg.resize(g.num_local_vertices(), 0);
        m_ghost_is_active.resize(g.num_ghost_vertices(), 1);
        m_candidate_of_node.resize(g.num_local_vertices(), EdgeE());

        m_active_nodes.resize(g.num_local_vertices());
        std::iota(std::begin(m_active_nodes), std::end(m_active_nodes), 0);

        // buffers for messages that are used to send candidate requests of
        // ghostvertices

        m_candidate_messages_of_proc.resize(
            m_procs, make_pair(std::vector<Message>(), MPI_Request()));
        m_buffer_message.resize(
            m_procs, make_pair(std::vector<Message>(), MPI_Request()));

        g.sort_adj();
        m_curr_candidate_id.resize(g.num_local_vertices(), 0);

        m_num_active_vertices = g.num_local_vertices();

        for (NodeIDType i = 0; i < g.num_local_vertices(); ++i) {
            m_active_deg[i] = g.getDegree(i);
        }
#ifdef LOGGING
        for (NodeIDType i : m_active_nodes) {

            printf("[%d] v%d: active deg %lld, colors: ", m_proc_id, i,
                   m_active_deg[i]);
            std::cout << m_color_available[i] << std::endl;
        }

#endif
        m_round = 0;

        while (m_num_active_vertices) {
#ifdef LOGGING
            printf("[%d] Start round %d\n", m_proc_id, m_round);
#endif
            setCandidates(g);

            clear_messages(m_candidate_messages_of_proc);
            getMsgFromBuffer(g);

            for (NodeIDType i : m_active_nodes) {
                if (v_is_done(g, i)) {
                    continue;
                }
                auto &c = m_candidate_of_node[i];

                if (g.is_local(c.partner)) {
                    if (i < c.partner &&
                        m_candidate_of_node[c.partner].partner == i) {
                        // found match

                        auto cc = color_set::common_colors(
                                      m_color_available[i],
                                      m_color_available[c.partner])
                                      .find_first();

                        if (cc != color_set::npos) {
                            // Not colorable edge

                            m_color_available[i].setOff(cc);
                            m_color_available[c.partner].setOff(cc);

#ifdef LOGGING
                            printf(ANSI_COLOR_MAGENTA
                                   "[%d] Matched local (v%d, "
                                   "v%d)\n" ANSI_COLOR_RESET,
                                   m_proc_id, g.local_vertex_id_to_global_id(i),
                                   g.local_vertex_id_to_global_id(c.partner));
#endif

                            add_to_matching(g, i, c.partner, c.weight,
                                            matchings, cc);
                            if (m_set_color) {
                                g.colorEdge(i, m_curr_candidate_id[i]);
                                g.colorEdge(c.partner,
                                            m_curr_candidate_id[c.partner]);
                            }
                        }

                        // Update candidates
                        ++m_curr_candidate_id[i];
                        --m_active_deg[i];

                        ++m_curr_candidate_id[c.partner];
                        --m_active_deg[c.partner];
                    }

                } else {
                    auto p = g.get_proc_of_ghost_vertex(c.partner);
                    if (!g.is_active_partner(p)) {
                        throw std::logic_error("NO");
                    }
                    // prepare message
                    m_candidate_messages_of_proc[g.get_proc_of_ghost_vertex(
                                                     c.partner)]
                        .first.push_back(
                            Message(g.local_vertex_id_to_global_id(i),
                                    g.local_vertex_id_to_global_id(c.partner),
                                    m_color_available[i].get_bits()));
                }
            }

            unsigned int count_msg_send_receive = 0;
            send_messages<Message>(g, count_msg_send_receive,
                                   m_candidate_messages_of_proc,
                                   MPI_MSG_TAG::get_candidates_tag(m_round),
                                   CANDIDATE_MESSAGE_TYPE);

            receive_ghost_candidates(g, count_msg_send_receive, matchings);

            ++m_round;

            // if (m_round == 15) {
            //     throw std::logic_error("DEBUG");
            // }
        }

        finished(g, matchings);

        // printf("[%d]: round %d, size %d, weight %.2f\n", m_proc_id, m_round,
        //        m_local_matching_size, m_local_matching_weight);

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Reduce(&m_local_matching_size, &total_global_matching_size, 1,
                   MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

        MPI_Reduce(&m_local_matching_weight, &total_global_matching_weight, 1,
                   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        MPI_Reduce(&m_round, &total_number_of_rounds, 1, MPI_INT, MPI_MAX, 0,
                   MPI_COMM_WORLD);
    }

  private:
    void receive_ghost_candidates(Graph &g, unsigned int &msgs_to_receive,
                                  std::vector<std::vector<Edge>> &matchings) {

        MPI_Status msg_status;

        int tag = MPI_MSG_TAG::get_candidates_tag(m_round);

        while (msgs_to_receive > 0) {
            MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &msg_status);

            int msg_count;
            MPI_Get_count(&msg_status, CANDIDATE_MESSAGE_TYPE, &msg_count);

            int src = msg_status.MPI_SOURCE;

            Message *messages = new Message[msg_count];

            MPI_Recv(messages, msg_count, CANDIDATE_MESSAGE_TYPE, src, tag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#ifdef LOGGING
            printf("[%d] Receive msg (size %d) with tag %d from p%d\n",
                   m_proc_id, msg_count, tag, src);
#endif
            msgs_to_receive--;

            for (int i = 0; i < msg_count; ++i) {
                Message msg = messages[i];

                if (msg.src == (NodeIDType)-1) {
                    // Proc is done!
                    g.proc_is_finished(src);
#ifdef LOGGING
                    printf(ANSI_COLOR_YELLOW
                           "[%d] Received p%d is done\n" ANSI_COLOR_RESET,
                           m_proc_id, src);

#endif
                    continue;
                }
                NodeIDType local_ghost_id =
                    g.global_vertex_id_to_local_id_of_ghost_vertex(msg.src);

                if (msg.candidate == (NodeIDType)-1) {
                    // Ghost is done
                    m_ghost_is_active[local_ghost_id - g.num_local_vertices()] =
                        0;
#ifdef LOGGING
                    printf(ANSI_COLOR_YELLOW "[%d] Received from p%d that v%d "
                                             "is done\n" ANSI_COLOR_RESET,
                           m_proc_id, src, msg.src);

#endif
                    continue;
                }

                NodeIDType local_id =
                    g.global_vertex_id_to_local_id(msg.candidate);

                if (m_candidate_of_node[local_id].partner == local_ghost_id) {
                    auto cc = color_set::common_colors(
                                  m_color_available[local_id], msg.colors)
                                  .find_first();

                    if (cc != color_set::npos) {
                        m_color_available[local_id].setOff(cc);

#ifdef LOGGING
                        printf(
                            ANSI_COLOR_MAGENTA
                            "[%d] Matched cross (v%d, v%d)\n" ANSI_COLOR_RESET,
                            m_proc_id, msg.candidate, msg.src);
#endif

                        if (m_set_color) {
                            g.colorEdge(local_id,
                                        m_curr_candidate_id[local_id]);
                        }

                        if (msg.candidate < msg.src)
                            add_to_matching(
                                g, local_id, local_ghost_id,
                                m_candidate_of_node[local_id].weight, matchings,
                                cc);
                    }

                    // Update candidate
                    ++m_curr_candidate_id[local_id];
                    --m_active_deg[local_id];

#ifdef LOGGING
                    // std::cout << "Local colors: " <<
                    // m_color_available[local_id]
                    //           << ", Received: ";
                    // for (int i = 7; i >= 0; i--) {
                    //     // Extract the i-th bit using bitwise AND and shift
                    //     long int bit = (msg.colors >> i) & 1;
                    //     printf("%ld", bit);
                    // }

                    // std::cout << ", Together: "
                    //           << color_set::common_colors(
                    //                  m_color_available[local_id], msg.colors)
                    //           << std::endl;

                    printf("[%d] Matched (cross) (v%d (%d), v%d (%d), c%d)\n",
                           m_proc_id, msg.src, local_id, msg.candidate,
                           local_ghost_id, cc);
#endif
                }

#ifdef LOGGING
                else {
                    printf(
                        ANSI_COLOR_BLUE
                        "[%d] Received non candidate (v%d,v%d), local candiate "
                        "is local %d\n" ANSI_COLOR_RESET,
                        m_proc_id, msg.candidate, msg.src,

                        m_candidate_of_node[local_id].partner);
                }

#endif
            }

            delete[] messages;
        }
    }

    void finished(Graph &g, std::vector<std::vector<Edge>> &matchings) {
        clear_messages(m_candidate_messages_of_proc);
        getMsgFromBuffer(g);

        // Send finish message!
        for (int p = 0; p < m_procs; ++p) {
            if (g.is_active_partner(p)) {
                m_candidate_messages_of_proc[p].first.push_back(
                    Message(-1, -1, 0));
            }
        }

        unsigned int count_msg_send_receive = 0;
        send_messages<Message>(
            g, count_msg_send_receive, m_candidate_messages_of_proc,
            MPI_MSG_TAG::get_candidates_tag(m_round), CANDIDATE_MESSAGE_TYPE);
        receive_ghost_candidates(g, count_msg_send_receive, matchings);
        clear_messages(m_candidate_messages_of_proc);
    }

    template <typename Type>
    void send_messages(Graph &g, unsigned int &send_msg_count,
                       std::vector<std::pair<std::vector<Type>, MPI_Request>>
                           &messages_of_proc,
                       int tag, MPI_Datatype msg_type) {
        send_msg_count = 0;

        for (int p = 0; p < m_procs; ++p) {
            if (p == m_proc_id)
                continue;

            if (!messages_of_proc[p].first.empty()) {
#ifdef LOGGING
                printf("[%d] Send msg (size %ld) with tag %d to p%d\n",
                       m_proc_id, messages_of_proc[p].first.size(), tag, p);
#endif
                MPI_Isend(&messages_of_proc[p].first[0],
                          messages_of_proc[p].first.size(), msg_type, p, tag,
                          MPI_COMM_WORLD, &messages_of_proc[p].second);
                ++send_msg_count;
            } else if (g.is_active_partner(p)) {
// we have to send empty message to active partner, because they
// expect a message from each active partner even if this
// message is empty!!
#ifdef LOGGING
                printf("[%d] Send empty msg with tag %d to p%d\n", m_proc_id,
                       tag, p);
#endif
                MPI_Isend(0, 0, msg_type, p, tag, MPI_COMM_WORLD,
                          &messages_of_proc[p].second);
                ++send_msg_count;
            }
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
        for (int p = 0; p < m_procs; ++p) {
            if (!messages_of_proc[p].first.empty()) {
                MPI_Wait(&messages_of_proc[p].second, MPI_STATUS_IGNORE);
                messages_of_proc[p].first.clear();
            }
        }
    }

  private:
    inline char v_is_done(Graph &g, NodeIDType local_id) {
        if (g.is_local(local_id)) {
            return !m_active_deg[local_id] ||
                   m_color_available[local_id].none();
        } else {
            return !m_ghost_is_active[local_id - g.num_local_vertices()];
        }
    }

    inline void add_to_matching(Graph &g, NodeIDType n1, NodeIDType n2,
                                WeightType weight,
                                std::vector<std::vector<Edge>> &matchings,
                                uint_fast8_t color) {

        auto g_n1 = g.local_vertex_id_to_global_id(n1);
        auto g_n2 = g.local_vertex_id_to_global_id(n2);

        // push back edge the since this fuction in only called by local
        // edges or one pf the processes of an cross edge
        matchings[color].push_back(Edge(g_n1, g_n2, weight));

        ++m_local_matching_size;
        m_local_matching_weight += weight;
        // printf("Add edge: %d, %d, %.2f\n", n1, n2, weight);
    }

    inline void setCandidates(Graph &g) {

        auto it = m_active_nodes.begin();
        while (it < m_active_nodes.end()) {

            if (!updateCandidate(g, it)) {
                // node was deactiveded, pointer was updated
                continue;
            }

            ++it;
        }
    }

    inline bool updateCandidate(Graph &g,
                                std::vector<NodeIDType>::iterator it) {
        NodeIDType &node = *it;

        EdgeE &c = m_candidate_of_node[node];
        auto doUpdate = true;
        while (doUpdate) {
            c = g.getNextNeightbor(node, m_curr_candidate_id[node]);

#ifdef LOGGING
            printf(ANSI_COLOR_CYAN
                   "[%d] Choose new neighbor for v%d: v%d\n" ANSI_COLOR_RESET,
                   m_proc_id, g.local_vertex_id_to_global_id(node),
                   c.partner == (NodeIDType)-1
                       ? c.partner
                       : g.local_vertex_id_to_global_id(c.partner));
#endif
            doUpdate = c.partner != (NodeIDType)-1 && v_is_done(g, c.partner);
            if (doUpdate) {
                ++m_curr_candidate_id[node];
                --m_active_deg[node];
            }
        };

        if (v_is_done(g, node)) {
#ifdef LOGGING
            printf(ANSI_COLOR_RED "[%d] Vertex %d is done\n" ANSI_COLOR_RESET,
                   m_proc_id, node);
#endif
            --m_num_active_vertices;
            m_active_deg[node] = 0;
            // Reset candidate
            if (c.partner != (NodeIDType)-1)
                c = EdgeE();

            // Need to inform neighbors
            auto ps = g.getAdjProcs(node);
            for (auto &p : ps) {

                m_buffer_message[p].first.push_back(
                    Message(g.local_vertex_id_to_global_id(node), -1, 0));

#ifdef LOGGING
                printf(ANSI_COLOR_YELLOW
                       "[%d] Inform p%d that v%d is done\n" ANSI_COLOR_RESET,
                       m_proc_id, p, g.local_vertex_id_to_global_id(node));

#endif
            }

            // delete from active nodes: swap to end and delete end to keep
            // consectuve memory
            m_active_nodes[it - m_active_nodes.begin()] =
                m_active_nodes[m_active_nodes.size() - 1];
            m_active_nodes.erase(m_active_nodes.end() - 1);

            return false;
        }

        return true;
    }

    inline void getMsgFromBuffer(Graph &g) {
        // @todo improve to no copy
        // m_candidate_messages_of_proc = m_buffer_message;
        for (int p = 0; p < m_procs; ++p) {
            if (g.is_active_partner(p)) {
                m_candidate_messages_of_proc[p].first.insert(
                    m_candidate_messages_of_proc[p].first.end(),
                    m_buffer_message[p].first.begin(),
                    m_buffer_message[p].first.end());
            }
            m_buffer_message[p].first.clear();
        }
    }

    MPI_Datatype CANDIDATE_MESSAGE_TYPE;
    MPI_Datatype DONE_MESSAGE_TYPE;

    unsigned int m_local_matching_size = 0;
    double m_local_matching_weight = 0.0;

    NodeIDType m_num_active_vertices;

    int m_procs, m_proc_id;

    unsigned int m_round;

    std::vector<std::pair<std::vector<Message>, MPI_Request>>
        m_candidate_messages_of_proc;
    std::vector<std::pair<std::vector<Message>, MPI_Request>> m_buffer_message;

    std::vector<EdgeE> m_candidate_of_node;
    std::vector<EdgeIDType> m_curr_candidate_id;
    std::vector<color_set> m_color_available;
    std::vector<EdgeIDType> m_active_deg;
    std::vector<char> m_ghost_is_active;

    std::vector<NodeIDType> m_active_nodes;

    bool m_set_color = false;
};

} // namespace parkec

#endif
