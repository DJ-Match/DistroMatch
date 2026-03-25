/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef _ADJ_GRAPH_
#define _ADJ_GRAPH_

#include <algorithm>
#include <concepts>
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
#include <util/hash_functions.h>
#include <util/math_funcs.h>
#include <util/mpi_tags.h>

namespace parkec {

struct Info {
    long long partner;
    float weight;
};

inline bool heapComp(Info left, Info right) {
    return (left.weight < right.weight ||
            (left.weight == right.weight && left.partner < right.partner));
}

class alignas(64) Node {
  public:
    long maxSize;
    long curSize;
    Info *heap;
    Info minEntry;

    void print();
    inline long long min_id() {
        return (curSize > 0 && curSize == maxSize) ? heap[0].partner
                                                   : (long long)-1;
    }

    inline float min_weight() {
        return (curSize > 0 && curSize == maxSize) ? heap[0].weight : 0.0;
    }

    inline long find_id(long long idx) {
        for (long i = 0; i < curSize; i++)
            if (heap[i].partner == idx)
                return 1;
        return 0;
    }

    void AddHeap(float wt, long long idx) {
        if (curSize == maxSize) {
            if (maxSize > 2) {
                // heap[0].weight=wt;
                // heap[0].partner=idx;

                /// Only heapify one branch of the heap tree
                long long small, ri, li, pi = 0;
                long long done = 0;

                if (heap[2].weight > heap[1].weight ||
                    (heap[2].weight == heap[1].weight &&
                     heap[2].partner > heap[1].partner))
                    small = 1;
                else
                    small = 2;

                if (wt > heap[small].weight ||
                    (wt == heap[small].weight && idx > heap[small].partner)) {
                    heap[0].weight = heap[small].weight;
                    heap[0].partner = heap[small].partner;
                    heap[small].weight = wt;
                    heap[small].partner = idx;
                } else {
                    heap[0].weight = wt;
                    heap[0].partner = idx;
                }

                pi = small;
                while (!done) {
                    li = 2 * pi + 1;
                    ri = 2 * pi + 2;
                    small = pi;

                    if (li < maxSize &&
                        (heap[li].weight < heap[small].weight ||
                         (heap[li].weight == heap[small].weight &&
                          heap[li].partner < heap[small].partner)))
                        small = li;
                    if (ri < maxSize &&
                        (heap[ri].weight < heap[small].weight ||
                         (heap[ri].weight == heap[small].weight &&
                          heap[ri].partner < heap[small].partner)))
                        small = ri;

                    if (small != pi) {
                        wt = heap[pi].weight;
                        idx = heap[pi].partner;

                        heap[pi].weight = heap[small].weight;
                        heap[pi].partner = heap[small].partner;

                        heap[small].weight = wt;
                        heap[small].partner = idx;
                    } else
                        done = 1;

                    pi = small;
                }
            } else {
                if (maxSize != 1 &&
                    (wt > heap[1].weight ||
                     (wt == heap[1].weight && idx > heap[1].partner))) {
                    heap[0].weight = heap[1].weight;
                    heap[0].partner = heap[1].partner;
                    heap[1].weight = wt;
                    heap[1].partner = idx;
                } else {
                    heap[0].weight = wt;
                    heap[0].partner = idx;
                }
            }

            minEntry.partner = heap[0].partner;
            minEntry.weight = heap[0].weight;
        } else {
            heap[curSize].weight = wt;
            heap[curSize].partner = idx;
            curSize++;
            if (curSize == maxSize) {
                std::sort(heap, heap + curSize, heapComp);
                minEntry.weight = heap[0].weight;
                minEntry.partner = heap[0].partner;
            }
        }
    }
};

class BSuitor;

/// _NodeIDType should be at least unsigned long if there are more than 2^32 - 1
/// nodes _EdgeIDType should be at least unsigned long if there are more than
/// 2^32 - 1 edges
template <typename _WeightType = double, typename _NodeIDType = unsigned int,
          typename _EdgeIDType = long long>
class DisAdjGraph : public IGraph<> {
    friend class BSuitor;
    // typedef _WeightType WeightType;
    // typedef _NodeIDType NodeIDType;
    // typedef _EdgeIDType EdgeIDType;

  public:
    // using WeightType = _WeightType;
    // using NodeIDType = _NodeIDType;
    // using EdgeIDType = _EdgeIDType;
    // using EdgeE = EdgeType;
    typedef parkec::EdgeE<WeightType, NodeIDType> EdgeE;

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

        // printf("Added local edge (v%d(%d),v%d(%d),w%.2f)\n",
        //        e.n1 + m_first_global_vertex, e.n1, e.n2 +
        //        m_first_global_vertex, e.n2, e.weight);

        // printf("Added local edge (v%d(%d),v%d(%d),w%.2f)\n",
        //        e.n2 + m_first_global_vertex, e.n2, e.n1 +
        //        m_first_global_vertex, e.n1, e.weight);

        m_init_edges[e.n1].push_back(EdgeE(e.n2, e.weight));
        m_init_edges[e.n2].push_back(EdgeE(e.n1, e.weight));

        ++m_num_local_edges;
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

        // printf("Added cross edge (v%d(%d),v%d(%d),w%.2f)\n", e.n1,
        //        e.n1 - m_first_global_vertex, e.n2, local_ghost_id, e.weight);

        // Adjust local ids
        e.n1 -= m_first_global_vertex;
        e.n2 = local_ghost_id;

        m_init_edges[e.n1].push_back(EdgeE(e.n2, e.weight));
        ++m_num_cross_edges;
    }

    inline virtual void finalize_initialization() {
        IGraph::finalize_initialization();

        m_edges_adj.reserve(m_num_local_edges * 2 + m_num_cross_edges);
        for (NodeIDType i = 0; i < m_num_local_vertices; ++i) {
            m_start_edges_of_v[i + 1] =
                m_start_edges_of_v[i] + m_init_edges[i].size();
            m_end_edges_of_v[i] = m_start_edges_of_v[i + 1];

            m_edges_adj.insert(m_edges_adj.end(), m_init_edges[i].begin(),
                               m_init_edges[i].end());
        }

        // m_init_edges.clear();
    }

    inline virtual void
    start_initialization(std::vector<NodeIDType> &first_global_vertex_of_proc) {
        IGraph::start_initialization(first_global_vertex_of_proc);

        m_init_edges.resize(m_num_local_vertices);
        m_start_edges_of_v.resize(m_num_local_vertices + 1, 0);
        m_end_edges_of_v.resize(m_num_local_vertices, 0);
        m_first_global_vertex_of_proc = first_global_vertex_of_proc;
    }

    inline int get_proc_of_global_vertex(const NodeIDType global_id) const {
        // guess proc_id based on equal distribution
        unsigned int proc_id =
            std::floor(global_id / (m_num_global_vertices / m_procs));
        proc_id = std::min(proc_id, (unsigned int)m_procs - 1);

        // correct guessed proc_id
        while (global_id < m_first_global_vertex_of_proc[proc_id])
            --proc_id;
        while (global_id >= m_first_global_vertex_of_proc[proc_id + 1])
            ++proc_id;

        return proc_id;
    }

    inline void print() const {
        for (NodeIDType i = 0; i < m_num_local_vertices; ++i) {
            printf("(%d: ", local_vertex_id_to_global_id(i));
            for (EdgeIDType e_id = m_start_edges_of_v[i];
                 e_id < m_end_edges_of_v[i]; ++e_id) {
                auto &e = m_edges_adj[e_id];
                printf("%d, w%.2f; ", local_vertex_id_to_global_id(e.partner),
                       e.weight);
            }
            printf("\n");
        }
    }

    inline void sort_adj() {
        for (NodeIDType i = 0; i < m_num_local_vertices; ++i) {
            std::sort(m_edges_adj.begin() + m_start_edges_of_v[i],
                      m_edges_adj.begin() + m_end_edges_of_v[i],
                      [this, i](const EdgeE &a, const EdgeE &b) {
                          if (a.weight != b.weight)
                              return a.weight > b.weight;

                          // Break ties!
                          unsigned int gid_i = local_vertex_id_to_global_id(i);
                          unsigned int gid_a =
                              local_vertex_id_to_global_id(a.partner);
                          unsigned int gid_b =
                              local_vertex_id_to_global_id(b.partner);

                          unsigned int old_hash = hash(gid_a ^ gid_i);
                          unsigned int new_hash = hash(gid_b ^ gid_i);

                          if (new_hash != old_hash)
                              return new_hash < old_hash;

                          return (gid_a ^ gid_i) < (gid_b ^ gid_i);
                      });
        }
    }
    inline void sort_adj_by_id() {
        for (NodeIDType i = 0; i < m_num_local_vertices; ++i) {
            std::sort(m_edges_adj.begin() + m_start_edges_of_v[i],
                      m_edges_adj.begin() + m_end_edges_of_v[i],
                      [this](const EdgeE &a, const EdgeE &b) {
                          return a.weight != b.weight
                                     ? a.weight > b.weight
                                     : local_vertex_id_to_global_id(a.partner) >
                                           local_vertex_id_to_global_id(
                                               b.partner);
                      });
        }
    }

    inline EdgeE &getNextNeightbor(NodeIDType localNode, NodeIDType next) {

        if (next >=
            m_end_edges_of_v[localNode] - m_start_edges_of_v[localNode]) {
            return dummy_e;
        }
        return m_edges_adj[m_start_edges_of_v[localNode] + next];
    }

    inline const EdgeIDType getDegree(NodeIDType localNode) {
        return m_end_edges_of_v[localNode] - m_start_edges_of_v[localNode];
    }

    inline absl::flat_hash_set<int> getAdjProcs(NodeIDType localNode) {

        absl::flat_hash_set<int> procs_adj_to_node;
        for (EdgeIDType e_id = m_start_edges_of_v[localNode];
             e_id < m_end_edges_of_v[localNode]; ++e_id) {
            auto &v = m_edges_adj[e_id].partner;
            if (!is_local(v)) {
                auto p = get_proc_of_ghost_vertex(v);
                if (is_active_partner(p))
                    procs_adj_to_node.insert(p);
            }
        }
        return procs_adj_to_node;
    }

    void proc_is_finished(int proc) { m_active_cross_edges_of_proc[proc] = 0; }
    void activate_all_procs() {
        for (int p = 0; p < m_procs; ++p)
            if (p != m_proc_id)
                m_active_cross_edges_of_proc[p] = 1;
    }

    void reduceGraph(Node *S_Queue) {
        // Make sure adj is sorted
        sort_adj_by_id();
        unsigned int start_i;
        for (NodeIDType i = 0; i < m_num_local_vertices; ++i) {
            start_i = m_start_edges_of_v[i];
            auto gid = i + m_first_global_vertex;
            auto n_matched = S_Queue[gid].curSize;

            auto heap = S_Queue[gid].heap;
            // Just make sure it is sorted!
            std::sort(heap, heap + n_matched, heapComp);
            auto idx = 0;

            for (NodeIDType j = n_matched - 1; j < n_matched; --j) {
                auto lid_partner =
                    global_vertex_id_to_local_id(heap[j].partner);

                while (m_edges_adj[start_i + idx].partner != lid_partner) {
                    ++idx;
                }

                if (idx != (n_matched - 1 - j)) {
                    // Need to swap to correct position; coorect postition is
                    // reverse to heap
                    std::swap(m_edges_adj[start_i + n_matched - 1 - j],
                              m_edges_adj[start_i + idx]);
                }
                ++idx;
            }
            m_end_edges_of_v[i] = start_i + n_matched;
        }
    }
    // template <typename std::enable_if<
    //               std::is_same<decltype(EdgeType : color),
    //               uint_fast8_t>::value, int>::type = 0>
    void colorEdge(NodeIDType localNode, NodeIDType neighborId) {

        if (neighborId >=
            m_end_edges_of_v[localNode] - m_start_edges_of_v[localNode]) {
            return;
        }
        m_edges_adj[m_start_edges_of_v[localNode] + neighborId].color = 1;
    }

    // template <typename std::enable_if<
    //               std::is_same<decltype(EdgeType : color),
    //               uint_fast8_t>::value, int>::type = 0>
    void reduceGraph() {
        for (NodeIDType localNode = 0; localNode < m_num_local_vertices;
             ++localNode) {
            EdgeIDType e_id = m_start_edges_of_v[localNode];
            while (e_id < m_end_edges_of_v[localNode]) {

                if (m_edges_adj[e_id].color != UINT_FAST8_MAX) {
                    // Deactive colored edge!
                    --m_end_edges_of_v[localNode];
                    std::swap(m_edges_adj[e_id],
                              m_edges_adj[m_end_edges_of_v[localNode]]);
                    continue;
                }

                ++e_id;
            }
        }
    }

    bool valid_k_matchings(std::vector<std::vector<Edge>> &matchings,
                           bool check_maximal = true) override {

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

        // Reset ends
        for (NodeIDType v = 0; v < m_num_local_vertices - 1; ++v) {
            m_end_edges_of_v[v] = m_start_edges_of_v[v + 1];
        }
        m_end_edges_of_v[m_num_local_vertices - 1] = m_edges_adj.size();

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

                NodeIDType local_n1 = global_vertex_id_to_local_id(e.n1);
                NodeIDType local_n2 = global_vertex_id_to_local_id(e.n2);

                check_if_colored(local_n1, local_n2);

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

                NodeIDType local_node = global_vertex_id_to_local_id(edge.n2);
                NodeIDType local_ghost = global_vertex_id_to_local_id(edge.n1);

                check_if_colored(local_node, local_ghost);

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

        for (NodeIDType localNode = 0; localNode < m_num_local_vertices;
             ++localNode) {
            for (EdgeIDType e_id = m_start_edges_of_v[localNode];
                 e_id < m_end_edges_of_v[localNode]; ++e_id) {
                auto &e = m_edges_adj[e_id];
                if (!is_local(e.partner)) {
                    nodes_for_proc[get_proc_of_ghost_vertex(e.partner)].insert(
                        localNode);
                }
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
            if (n1 > n2) {
                auto temp = n2;
                n2 = n1;
                n1 = temp;
            }
            if (colored_edges.contains(std::make_pair(n1, n2))) {
                // edge is colored!
                return;
            }
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

        for (NodeIDType localNode = 0; localNode < m_num_local_vertices;
             ++localNode) {
            for (EdgeIDType e_id = m_start_edges_of_v[localNode];
                 e_id < m_end_edges_of_v[localNode]; ++e_id) {
                auto &e = m_edges_adj[e_id];
                check_if_colorable(localNode, e.partner, e.weight);
            }
        }

        // Reduce local results across all processes
        int global_result;
        MPI_Allreduce(&local_result, &global_result, 1, MPI_INT, MPI_LAND,
                      MPI_COMM_WORLD);

        return global_result;
    }

  private:
    EdgeE dummy_e = EdgeE();
    EdgeIDType m_num_edges = 0;
    std::vector<EdgeIDType> m_start_edges_of_v;
    std::vector<EdgeIDType> m_end_edges_of_v;
    std::vector<EdgeE> m_edges_adj;
    std::vector<std::vector<EdgeE>> m_init_edges;
    std::vector<NodeIDType> m_first_global_vertex_of_proc;
};

} // namespace parkec

#endif