/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 * This file is based on LocalMaxMatching (MIT Licence, Marcel Birn)
 *   https://github.com/LocalMaxMatching/LocalMaxMatching
 *
 *****************************************************************************/

#ifndef _PARALLEL_EDGE_GRAPH_
#define _PARALLEL_EDGE_GRAPH_

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
#include "graph.h"
#include <util/color_set.h>
#include <util/math_funcs.h>
#include <util/mpi_tags.h>

namespace parkec {

template <typename T>
void swap(const unsigned long &pos1, const unsigned long &pos2,
          std::vector<T> &vec) {
    T tmp = vec[pos1];
    vec[pos1] = vec[pos2];
    vec[pos2] = tmp;
}

/// _NodeIDType should be at least unsigned long if there are more than 2^32 - 1
/// nodes _EdgeIDType should be at least unsigned long if there are more than
/// 2^32 - 1 edges
template <typename _WeightType = double, typename _NodeIDType = unsigned int,
          typename _EdgeIDType = long long>
class ParallelEdgeGraph : public IGraph<> {

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
    ParallelEdgeGraph()
        : IGraph(), m_dummy_vertex(std::numeric_limits<NodeIDType>::max()) {
        m_dummy_edge.n1 = std::numeric_limits<NodeIDType>::max();
        m_dummy_edge.n2 = std::numeric_limits<NodeIDType>::max();
        m_dummy_edge.weight = min_val<WeightType>();
    }

    virtual ~ParallelEdgeGraph() {}

    inline EdgeIDType num_all_edges() const {
        return m_num_local_edges + m_num_cross_edges;
    }

    /*!
     * \brief Returns the edge corresponding to the given edge_id-iterator.
     *
     * Runtime: O(1)
     *
     * @param edge_id
     * @return			Returns a const reference to the edge
     * with the given ID.
     */
    inline const Edge &get_edge(const edge_iterator &edge_it) const {
        return *edge_it;
    }

    /*!
     * \brief Returns the edge corresponding to the given edge_id-iterator.
     *
     * Runtime: O(1)
     *
     * @param edge_id
     * @return			Returns a reference to the edge with the
     * given ID.
     */
    inline Edge &get_edge(const edge_iterator &edge_it) { return *edge_it; }

    /*!
     * \brief Returns the edge corresponding to the given edge_id.
     *
     * Make sure that the edge_id is still valid. It's valid as long
     * as you don't deactivate any edges. I.e. deactivating edges
     * makes currently use edge_ids invalid (it least some of them).
     *
     * Runtime: O(1)
     *
     * @param edge_id
     * @return			Returns a const reference to the edge
     * with the given ID.
     */
    inline const Edge &get_cross_edge(const EdgeIDType &edge_id) const {
        return m_cross_edges[edge_id];
    }
    inline Edge &get_cross_edge(const EdgeIDType &edge_id) {
        return m_cross_edges[edge_id];
    }
    inline const Edge &get_local_edge(const EdgeIDType &edge_id) const {
        return m_local_edges[edge_id];
    }
    inline Edge &get_local_edge(const EdgeIDType &edge_id) {
        return m_local_edges[edge_id];
    }

    /*!
     * \brief Returns corresponding id to the given iterator.
     *
     * Note that the id is NOT unique, and refers to cross/local edges
     * separately !! Use get_local_edge/get_cross_edge to receive the correct
     * edge ref
     *
     * Runtime: O(1)
     *
     * @param edge_it
     * @return	Returns the id
     */
    inline virtual EdgeIDType get_edge_id(const edge_iterator &edge_it) const {
        return is_local(edge_it->n2) ? edge_it - m_local_edges.begin()
                                     : edge_it - m_cross_edges.begin();
    }

    inline virtual EdgeIDType
    get_edge_id(const const_edge_iterator &edge_it) const {
        return is_local(edge_it->n2) ? edge_it - m_local_edges.begin()
                                     : edge_it - m_cross_edges.begin();
    }

    /*!
     * \brief Returns the weight of the edge referenced by the given
     * iterator.
     *
     * Don't use edge references!!!
     *
     * @param edge_it
     * @return			Returns the weight of the edge
     * referenced by the given iterator.
     */
    inline WeightType get_edge_weight(const edge_iterator &edge_it) const {
        return edge_it->weight;
    }

    inline WeightType get_edge_weight(const edge_iterator &edge_it) {
        return edge_it->weight;
    }

    inline virtual edge_iterator begin_active_local_edges() {
        return m_local_edges.begin();
    }
    inline const_edge_iterator begin_active_local_edges() const {
        return m_local_edges.begin();
    }

    inline virtual edge_iterator end_active_local_edges() {
        return m_local_edges.begin() + m_end_active_local_edges;
    }
    inline const_edge_iterator end_active_local_edges() const {
        return m_local_edges.begin() + m_end_active_local_edges;
    }

    inline virtual edge_iterator begin_active_edges() {
        throw std::logic_error("NOT implemented in this class");
        return m_local_edges.begin();
    }
    inline virtual edge_iterator end_active_edges() {
        throw std::logic_error("NOT implemented in this class");
        return m_local_edges.begin();
    }

    inline virtual edge_iterator begin_active_cross_edges() {
        return m_cross_edges.begin();
    }
    inline const_edge_iterator begin_active_cross_edges() const {
        return m_cross_edges.begin();
    }

    inline virtual edge_iterator end_active_cross_edges() {
        return m_cross_edges.begin() + m_end_active_cross_edges;
    }
    inline const_edge_iterator end_active_cross_edges() const {
        return m_cross_edges.begin() + m_end_active_cross_edges;
    }

    inline edge_iterator end_active_owned_cross_edges() {
        return m_cross_edges.begin() + m_end_active_owned_cross_edges;
    }
    inline const_edge_iterator end_active_owned_cross_edges() const {
        return m_cross_edges.begin() + m_end_active_owned_cross_edges;
    }

    inline edge_iterator get_cross_edge_it(const EdgeIDType edge_id) {
        return m_cross_edges.begin() + edge_id;
    }

    inline bool is_owed(edge_iterator e_it) const {
        return get_edge_id(e_it) < m_end_active_owned_cross_edges;
    }

    inline virtual edge_iterator deactivate_edge(const edge_iterator edge_it) {
        throw std::logic_error("NOT implemented in this class");
        return edge_it;
    }
    inline virtual edge_iterator
    edge_is_candidate(const edge_iterator edge_it) {
        throw std::logic_error("NOT implemented in this class");
        return edge_it;
    }
    inline virtual edge_iterator end_active_candidates() {
        throw std::logic_error("NOT implemented in this class");
        return m_local_edges.begin();
    }
    inline virtual void reset_candidates() {
        throw std::logic_error("NOT implemented in this class");
    }
    inline virtual edge_iterator
    deactivate_candidate(const edge_iterator edge_it) {
        throw std::logic_error("NOT implemented in this class");
    }

    /*!
     * \brief Deactivates the local edge, that is referenced by the given
     * iterator.
     *
     * The given edge iterator will reference a new active edge after this
     * operation, unless all edges are inactive.
     *
     * The new referenced edge will be an edge that hasn't been iterated
     * yet!!! (unless all edges are inactive).
     *
     * @param edge_it	An iterator to the local edge to be deactivated.
     *
     * @return			Returns an iterator-reference to the
     * deactivated local edge.
     */
    inline edge_iterator deactivate_local_edge(const edge_iterator edge_it) {
        EdgeIDType edge_ref = get_edge_id(edge_it);

        if (edge_ref >= m_end_active_local_edges) {
            return edge_it;
        }

        m_end_active_local_edges--;

        swap(m_end_active_local_edges, edge_ref, m_local_edges);

        return m_local_edges.begin() + m_end_active_local_edges;
    }

    /*!
     * \brief Deactivates the cross edge, that is referenced by the given
     * iterator.
     *
     * The given edge iterator will reference a new active edge after this
     * operation, unless all edges are inactive.
     *
     * The new referenced edge will be an edge that hasn't been iterated
     * yet!!! (unless all edges are inactive).
     *
     * @param edge_it	An iterator to the cross edge to be deactivated.
     *
     * @return			Returns an iterator-reference to the
     * deactivated cross edge.
     */
    inline edge_iterator
    deactivate_cross_edge(const edge_iterator edge_it,
                          bool deactivate_only_locally = false) {
        EdgeIDType edge_ref = get_edge_id(edge_it);

        if (edge_ref >= m_end_active_cross_edges) {
            return edge_it;
        }

        if (!deactivate_only_locally) {
            auto p = get_proc_of_ghost_vertex(m_cross_edges[edge_ref].n2);

            // .n2 is ghost - we ensure this in the constructor
            m_active_cross_edges_of_proc[p]--;

            m_count_active_partners -= (m_active_cross_edges_of_proc[p] == 0);
        }

        --m_end_active_cross_edges;

        if (edge_ref < m_end_active_owned_cross_edges) {

            --m_end_active_owned_cross_edges;
            swap(m_end_active_owned_cross_edges, edge_ref, m_cross_edges);

            if (m_end_active_owned_cross_edges != m_end_active_cross_edges) {
                swap(m_end_active_owned_cross_edges, m_end_active_cross_edges,
                     m_cross_edges);
            }
        } else {

            swap(m_end_active_cross_edges, edge_ref, m_cross_edges);
        }
        return m_cross_edges.begin() + m_end_active_cross_edges;
    }

    inline edge_iterator unown_edge(const edge_iterator edge_it) {
        EdgeIDType edge_ref = get_edge_id(edge_it);

        if (edge_ref >= m_end_active_owned_cross_edges) {
            return edge_it;
        }
        --m_end_active_owned_cross_edges;
        swap(m_end_active_owned_cross_edges, edge_ref, m_cross_edges);

        return m_cross_edges.begin() + m_end_active_owned_cross_edges;
    }

    // inline edge_iterator
    // deactivate_unsorted_edge_locally(const edge_iterator edge_it) {
    //     EdgeIDType edge_ref = get_edge_id(edge_it);

    //     if (edge_ref >= m_end_active_cross_edges) {
    //         return edge_it;
    //     }

    //     // this is the overall end of edges now
    //     m_end_active_cross_edges--;

    //     swap(m_end_active_cross_edges, edge_ref, m_local_edges);

    //     return m_cross_edges.begin() + m_end_active_cross_edges;
    // }

    /*!
     * \brief Returns true if there are active edges, otherwise false.
     *
     * @return Returns true if there are active edges, otherwise false.
     */
    inline virtual bool has_active_edges() const {
        return m_end_active_local_edges | m_end_active_cross_edges;
    }

    /*!
     * \brief Returns the number of active edges.
     *
     * @return Returns the number of active edges.
     */
    inline EdgeIDType num_active_edges() const {
        return m_end_active_local_edges + m_end_active_cross_edges;
    }

    WeightType get_weight() const;

    inline bool increase_active_cross_edges_of_proc(const int proc_id) {
        if (m_active_cross_edges_of_proc[proc_id] == 0) {
            m_count_active_partners++;
        }
        return ++m_active_cross_edges_of_proc[proc_id];
    }

    inline EdgeIDType
    get_active_cross_edges_count_of_proc(const int proc_id) const {
        return m_active_cross_edges_of_proc[proc_id];
    }

    inline NodeIDType
    global_vertex_id_to_local_id_check(const NodeIDType n_id) const {
        if (global_id_is_local(n_id)) {
            return global_vertex_id_to_local_id_of_local_vertex(n_id);
        } else {
            if (m_ghost_global_to_local_hash.contains(n_id))
                return global_vertex_id_to_local_id_of_ghost_vertex(n_id);
            else {
                printf("There is no ghost node glob v%d!\n", n_id);
                throw std::invalid_argument("Invalid request.");
                return -1;
            }
        }
    }

    inline virtual NodeIDType get_ghost_id(const NodeIDType n_id) const {
        throw std::logic_error("NOT implemented in this class");
        return -1;
    }

    inline virtual void print() const {
        for (EdgeIDType pos = 0; pos < m_num_local_edges; pos++) {
            const auto &e = m_local_edges[pos];
            printf("(%d<->%d, %f),", local_vertex_id_to_global_id(e.n1),
                   local_vertex_id_to_global_id(e.n2), e.weight);
        }
        for (EdgeIDType pos = 0; pos < m_num_cross_edges; pos++) {
            const auto &e = m_cross_edges[pos];
            printf("(%d<->%d, %f),", local_vertex_id_to_global_id(e.n1),
                   local_vertex_id_to_global_id(e.n2), e.weight);
        }
        printf("\n");
    }

    inline virtual void add_local_edge(Edge &e) override {

        if (m_initialized) {
            printf("Graph is already final!");
            exit(1);
        }
        if (!(global_id_is_local(e.n1) && global_id_is_local(e.n2))) {
            // This is not my edge!
            printf("[%d] Edge (%d, %d) does not belong the assigned vertex "
                   "range.",
                   m_proc_id, e.n1, e.n2);
            throw std::logic_error("Initialization error");
        }

        // Adjust local ids
        e.n1 -= m_first_global_vertex;
        e.n2 -= m_first_global_vertex;

        m_local_edges.push_back(e);

        if (e.weight > m_max_weight) {
            m_max_weight = e.weight;
        }

        // printf("Added local edge (v%d(%d),v%d(%d))\n",
        //        e.n1 + m_first_global_vertex, e.n1, e.n2 +
        //        m_first_global_vertex, e.n2);
    }

    inline virtual void add_cross_edge(Edge &e, int proc_of_ghost) {
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

        if (m_active_cross_edges_of_proc[proc_of_ghost] == 0) {
            ++m_count_active_partners;
        }

        ++m_active_cross_edges_of_proc[proc_of_ghost];

        // printf("Added cross edge (v%d(%d),v%d(%d))\n", e.n1,
        //        e.n1 - m_first_global_vertex, e.n2, local_ghost_id);

        // Adjust local ids
        e.n1 -= m_first_global_vertex;
        e.n2 = local_ghost_id;

        m_cross_edges.push_back(e);

        if (e.weight > m_max_weight) {
            m_max_weight = e.weight;
        }
    }

    inline virtual void finalize_initialization() {
        IGraph::finalize_initialization();

        m_num_local_edges = m_end_active_local_edges = m_local_edges.size();
        m_num_cross_edges = m_end_active_cross_edges =
            m_end_active_owned_cross_edges = m_cross_edges.size();

        // add one dummy node
        m_local_ghost_to_global.push_back(m_dummy_vertex);

        m_local_edges.push_back(m_dummy_edge);
        m_cross_edges.push_back(m_dummy_edge);
    }

    inline void sort_edges_by_weight() {
        std::sort(m_local_edges.begin(), m_local_edges.end() - 1,
                  [](Edge &e1, Edge &e2) { return e1.weight < e2.weight; });

        std::sort(m_cross_edges.begin(), m_cross_edges.end() - 1,
                  [](Edge &e1, Edge &e2) { return e1.weight < e2.weight; });
    }

    void virtual reset() {
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
    }

    bool valid_k_matchings(std::vector<std::vector<Edge>> &matchings,
                           bool check_maximal = true) override {

        // reset graph
        reset();

        int local_result = true;

        std::vector<color_set> colors_avilable(num_all_vertices(),
                                               color_set(matchings.size()));

        struct EdgeToSend {
            NodeIDType n1;
            NodeIDType n2;
            int color;
        };

        struct NodeToSend {
            NodeIDType node;
            color_set colors;
        };

        MPI_Datatype EDGE_TYPE;
        MPI_Type_contiguous(sizeof(EdgeToSend), MPI_BYTE, &EDGE_TYPE);
        MPI_Type_commit(&EDGE_TYPE);

        MPI_Datatype NODE_TYPE;
        MPI_Type_contiguous(sizeof(NodeToSend), MPI_BYTE, &NODE_TYPE);
        MPI_Type_commit(&NODE_TYPE);

        std::vector<std::pair<std::vector<EdgeToSend>, MPI_Request>>
            msg_edges_for_proc(
                m_procs, make_pair(std::vector<EdgeToSend>(), MPI_Request()));
        std::vector<std::pair<std::vector<NodeToSend>, MPI_Request>>
            msg_nodes_for_proc(
                m_procs, make_pair(std::vector<NodeToSend>(), MPI_Request()));

        absl::flat_hash_set<std::pair<NodeIDType, NodeIDType>> colored_edges;

        std::vector<absl::flat_hash_set<NodeIDType>> nodes_for_proc;
        nodes_for_proc.resize(m_procs);

        for (unsigned int c = 0; c < matchings.size(); ++c)
            for (Edge &e : matchings[c]) {
                e.color = c;
            }

        // Helper function to log and update matched nodes
        auto check_and_update_matched = [&](NodeIDType local_id,
                                            NodeIDType global_id, int color) {
            if (!colors_avilable[local_id][color]) {
                printf("proc %d: Vertex v%d (%d) is incident to more "
                       "than one edge with color %d!\n",
                       m_proc_id, global_id, local_id, color);
                local_result = false;
            }
            colors_avilable[local_id].setOff(color);
        };

        auto check_if_colored = [&](NodeIDType n1, NodeIDType n2) {
            if (n1 > n2) {
                auto temp = n2;
                n2 = n1;
                n1 = temp;
            }
            if (colored_edges.contains(std::make_pair(n1, n2))) {
                printf("proc %d: Edge (v%d,v%d) has more than one color!\n",
                       m_proc_id, n1, n2);
                local_result = false;
            }
            colored_edges.insert(std::make_pair(n1, n2));
        };

        // Process matching and prepare messages
        for (auto &matching : matchings)
            for (const Edge &e : matching) {

                // printf("p%d CHECK: (v%d,v%d, c%d)\n", m_proc_id, e.n1, e.n2,
                //        e.color);

                check_if_colored(e.n1, e.n2);

                NodeIDType local_n1 = global_vertex_id_to_local_id(e.n1);
                NodeIDType local_n2 = global_vertex_id_to_local_id(e.n2);

                check_and_update_matched(local_n1, e.n1, e.color);
                check_and_update_matched(local_n2, e.n2, e.color);

                if (is_ghost(local_n2)) {
                    auto p = get_proc_of_ghost_vertex(local_n2);

                    EdgeToSend edge;
                    edge.n1 = e.n1;
                    edge.n2 = e.n2;
                    edge.color = e.color;
                    msg_edges_for_proc[p].first.push_back(edge);
                }
            }

        std::function<void(EdgeToSend &)> process_req_edges =
            [this, &check_if_colored,
             &check_and_update_matched](EdgeToSend &edge) {
                // printf("p%d CHECK r: (v%d,v%d, c%d)\n", m_proc_id,
                //    edge.n1, edge.n2, edge.color);

                check_if_colored(edge.n1, edge.n2);

                NodeIDType local_node = global_vertex_id_to_local_id(edge.n2);
                NodeIDType local_ghost = global_vertex_id_to_local_id(edge.n1);

                check_and_update_matched(local_node, edge.n2, edge.color);
                check_and_update_matched(local_ghost, edge.n1, edge.color);
            };

        send_receive(msg_edges_for_proc, EDGE_TYPE, MPI_MSG_TAG::test_matching1,
                     process_req_edges);

        if (!check_maximal) {

            // Reduce local results across all processes
            int global_result;
            MPI_Allreduce(&local_result, &global_result, 1, MPI_INT, MPI_LAND,
                          MPI_COMM_WORLD);

            return global_result;
        }

        auto e_it = m_cross_edges.begin();
        while (e_it != m_cross_edges.end() - 1) {
            auto p = get_proc_of_ghost_vertex(e_it->n2);
            nodes_for_proc[p].insert(e_it->n1);
            ++e_it;
        }

        if (m_cross_edges.size() == 1 && m_num_ghost_vertices) {
            e_it = m_local_edges.begin();
            while (e_it != m_local_edges.end() - 1) {
                if (!is_local(e_it->n2)) {
                    auto p = get_proc_of_ghost_vertex(e_it->n2);
                    nodes_for_proc[p].insert(e_it->n1);
                }
                ++e_it;
            }
        }

        for (int p = 0; p < m_procs; ++p) {
            for (auto node_local : nodes_for_proc[p]) {
                NodeToSend node;
                node.node = local_vertex_id_to_global_id(node_local);
                node.colors = colors_avilable[node_local];
                msg_nodes_for_proc[p].first.push_back(node);
            }
        }

        std::function<void(NodeToSend &)> process_req_nodes =
            [this, &colors_avilable](NodeToSend &node) {
                NodeIDType local_node = global_vertex_id_to_local_id(node.node);

                // printf("proc %d: v%d had free colors: (", m_proc_id,
                // node.node); std::cout << colors_avilable[local_node] << ")
                // and received ("; std::cout << colors_avilable[local_node] <<
                // ") = now (";

                colors_avilable[local_node].make_intersect(node.colors);

                // std::cout << colors_avilable[local_node] << ")\n";
            };

        send_receive(msg_nodes_for_proc, NODE_TYPE, MPI_MSG_TAG::test_matching2,
                     process_req_nodes);

        auto check_if_colorable = [&](NodeIDType n1, NodeIDType n2,
                                      WeightType weight) {
            auto cc_set = color_set::common_colors(colors_avilable[n1],
                                                   colors_avilable[n2]);
            if (!(cc_set.find_first() == color_set::npos)) {

                printf("proc %d: Edge (v%d, v%d, %f) has free colors: (",
                       m_proc_id, local_vertex_id_to_global_id(n1),
                       local_vertex_id_to_global_id(n2), weight);

                std::cout << cc_set << ")\n";
                local_result = false;
            }
        };

        for (edge_iterator e_it = begin_active_local_edges();
             e_it != end_active_local_edges(); e_it++) {
            check_if_colorable(e_it->n1, e_it->n2, e_it->weight);
        }

        for (edge_iterator e_it = begin_active_cross_edges();
             e_it != end_active_cross_edges(); e_it++) {
            check_if_colorable(e_it->n1, e_it->n2, e_it->weight);
        }

        // Reduce local results across all processes
        int global_result;
        MPI_Allreduce(&local_result, &global_result, 1, MPI_INT, MPI_LAND,
                      MPI_COMM_WORLD);

        return global_result;
    }

    inline virtual void activate_edges() {
        m_end_active_local_edges = m_num_local_edges;
        m_end_active_cross_edges = m_num_cross_edges;
        m_end_active_owned_cross_edges = m_num_cross_edges;

        // count the active cross edges for each proc
        m_active_cross_edges_of_proc.assign(m_procs, 0);

        for (EdgeIDType e_id = 0; e_id < m_num_cross_edges; ++e_id) {
            // .n2 is ghost node
            m_active_cross_edges_of_proc[get_proc_of_ghost_vertex(
                m_cross_edges[e_id].n2)]++;
        }

        m_count_active_partners = 0;
        for (int p = 0; p < m_procs; p++) {
            if (m_active_cross_edges_of_proc[p] > 0) {
                m_count_active_partners++;
            }
        }
    }

    inline EdgeIDType getNumLocalEdges() const override {
        return m_initialized ? m_num_local_edges : m_local_edges.size();
    }
    inline EdgeIDType getNumCrossEdges() const override {
        return m_initialized ? m_num_cross_edges : m_cross_edges.size();
    }
    inline EdgeIDType getNumAllEdges() const override {
        return m_initialized ? m_num_local_edges + m_num_cross_edges
                             : m_local_edges.size() + m_cross_edges.size();
    }

    inline WeightType getMaxWeight() const { return m_max_weight; }

  protected:
    EdgeIDType m_end_active_local_edges;
    EdgeIDType m_end_active_cross_edges;
    EdgeIDType m_end_active_owned_cross_edges;
    EdgeIDType m_end_active_candidates;

    // first part local edges, second part cross edges
    std::vector<Edge> m_local_edges;
    std::vector<Edge> m_cross_edges;

    WeightType m_max_weight;

    Edge m_dummy_edge;
    NodeIDType m_dummy_vertex;
};

} // namespace parkec

#endif
