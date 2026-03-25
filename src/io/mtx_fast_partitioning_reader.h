/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef MTX_FAST_PART_READER
#define MTX_FAST_PART_READER

#include <algorithm>
#include <iterator>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

#include "datastructure/parallel_edge_graph.h"
#include "io/mtx_reader.h"
#include "io/reader.h"
#include "util/random.h"

#include "absl/container/flat_hash_set.h"

namespace graph_io {

template <class Graph, typename Graph::NodeIDType (*adjust_id)(
                           const typename Graph::NodeIDType) =
                           parkec::identity<typename Graph::NodeIDType>>
class MtxPartitionGraphReader : public MtxGraphReader<Graph> {
    typedef typename Graph::WeightType WeightType;
    typedef typename Graph::NodeIDType NodeIDType;
    typedef typename Graph::EdgeIDType EdgeIDType;
    typedef typename Graph::Edge Edge;

  public:
    struct EndPoint {
        NodeIDType node;
        WeightType weight;

        EndPoint(NodeIDType node, WeightType weight)
            : node(node), weight(weight) {}
    };

    void read_header(std::ifstream &in, NodeIDType &n) override {
        NodeIDType m;
        this->read_header_with_m(in, n, m);
        fennel_alpha = m / std::pow(n, fennel_gamma);
        m_sum_degrees = m * 2;
    }

    void read_graph(std::ifstream &in, Graph &g, bool use_partitioning,
                    NodeIDType num_vertices, unsigned int proc_id,
                    unsigned int procs,
                    std::vector<NodeIDType> &first_global_vertex_of_proc,
                    std::vector<NodeIDType> &vertex_mapping) override {

        NodeIDType src_vertex, tar_vertex;

        std::string line;
        unsigned int self_loops = 0;

        m_adj.clear();

        // note the insertion of dummy nodes to fit the all_gather
        unsigned int nodes_per_proc =
            num_vertices / procs + (num_vertices % procs > 0);

        // distribute nodes over processes
        NodeIDType first_in = nodes_per_proc * proc_id;
        NodeIDType last_ex = nodes_per_proc * (proc_id + 1);

        m_adj.resize(nodes_per_proc);
        m_node_weights.resize(nodes_per_proc * procs, 0);
        m_vertex_partition.resize(nodes_per_proc * procs, -1);

        while (std::getline(in, line)) {
            if (line[0] == '%') {
                // ignore comment line
                continue;
            }

            std::stringstream strs(line);

            strs >> src_vertex >> tar_vertex;

            // Ids should start with zero
            --src_vertex;
            --tar_vertex;

            double weight = 1;
            if (this->m_weighted) {
                strs >> weight;
            }

            auto src_is_local = src_vertex >= first_in && src_vertex < last_ex;
            auto tar_is_local = tar_vertex >= first_in && tar_vertex < last_ex;

            if (!(src_is_local || tar_is_local)) {
                continue;
            }

            if (src_vertex == tar_vertex) {
                // ignore self edges
                ++self_loops;
                continue;
            }

            if (tar_is_local) {
                m_adj[tar_vertex - first_in].push_back(
                    EndPoint(src_vertex, weight));
                ++m_node_weights[tar_vertex];
            }
            if (src_is_local) {
                m_adj[src_vertex - first_in].push_back(
                    EndPoint(tar_vertex, weight));
                ++m_node_weights[src_vertex];
            }
        }

        if (self_loops)
            printf("[%d] Ignored %d self loops\n", proc_id, self_loops);

        fennel_alpha *= std::pow(procs, fennel_gamma - 1);
        // * 1/std::pow(fennel_t_alpha, m_fennel_rounds);
        if (proc_id == 0)
            printf("Partition using fennel (fennel_alpha: %.8f, fennel_gamma: "
                   "%.3f, initial partitioning seed: %d)\n",
                   fennel_alpha, fennel_gamma, m_seed);
        double start_time = MPI_Wtime();

        MPI_Allgather(&m_node_weights[first_in], nodes_per_proc, MPI_INT,
                      &m_node_weights[0], nodes_per_proc, MPI_INT,
                      MPI_COMM_WORLD);
        m_sum_weights = 0;
        for (unsigned int i = 0; i < num_vertices; ++i) {
            m_sum_weights += m_node_weights[i];
        }

        // Start with initial partition
        std::vector<unsigned int> partition_current_buf;

        parkec::generate_initial_partitioning(
            partition_current_buf, nodes_per_proc, procs, proc_id, m_seed);

        MPI_Allgather(&partition_current_buf[0], nodes_per_proc, MPI_INT,
                      &m_vertex_partition[0], nodes_per_proc, MPI_INT,
                      MPI_COMM_WORLD);

        m_buffer.resize(procs);

        // Partitioning
        m_edges_in_proc.resize(procs, 0);
        m_partition_sizes.resize(procs, 0);
        m_p_scores.resize(procs, 0);

        for (int round = 0; round < m_fennel_rounds; ++round) {
            calcPartitionSizes();

            if (proc_id == 0) {
                std::cout << "--- Round: " << round << " ---\n";
                std::cout << "Partition sizes at start: (";
                for (unsigned int i = 0; i < procs; ++i) {
                    std::cout << m_partition_sizes[i] << ", ";
                }
                std::cout << "), Balance: " << calcBalance(procs)
                          << ", fennel_alpha: " << fennel_alpha << "\n";
            }

            for (unsigned int local_node = 0; local_node < m_adj.size();
                 ++local_node) {

                auto node = local_node + first_in;

                auto old_p = m_vertex_partition[node];

                auto new_p = old_p;

                if (m_adj[local_node].size() != 0) {
                    new_p = decide_proc_for_neighborhood(local_node, node,
                                                         proc_id, procs);
                }

                MPI_Allgather(&new_p, 1, MPI_INT, &m_buffer[0], 1, MPI_INT,
                              MPI_COMM_WORLD);

                for (unsigned int proc = 0; proc < procs; ++proc) {
                    auto updated_node = local_node + (proc * nodes_per_proc);

                    auto p_before = m_vertex_partition[updated_node];
                    auto p_now = m_buffer[proc];

                    m_partition_sizes[p_now] += m_node_weights[updated_node];
                    m_partition_sizes[p_before] -= m_node_weights[updated_node];

                    m_vertex_partition[updated_node] = p_now;
                }
            }

            fennel_alpha *= fennel_t_alpha;
        }

        m_edges_in_proc.clear();
        m_p_scores.clear();

        MPI_Barrier(MPI_COMM_WORLD);
        double end_time = MPI_Wtime();
        double elapsed_time_is_s(end_time - start_time);
        if (proc_id == 0)
            printf("Partitioning took %.5f s\n", elapsed_time_is_s);

        // Reset partition size to node count
        auto final_balance = calcBalance(procs);
        std::fill(m_partition_sizes.begin(), m_partition_sizes.end(), 0);
        for (unsigned int i = 0; i < num_vertices; ++i) {
            auto proc_of_vertex = m_vertex_partition[i];
            vertex_mapping.push_back(m_partition_sizes[proc_of_vertex] +
                                     proc_of_vertex * num_vertices);
            ++m_partition_sizes[proc_of_vertex];
        }

        if (proc_id == 0) {
            std::cout << "--- END " << " ---\n";
            std::cout << "Partition sizes: (";
            for (unsigned int i = 0; i < procs; ++i) {
                std::cout << m_partition_sizes[i] << ", ";
            }
            std::cout << "), Balance: " << final_balance << "\n";
        }

        // First: count nodes per proc and set ranges/vertex mapping
        for (unsigned int i = 0; i < procs - 1; ++i) {

            first_global_vertex_of_proc.push_back(
                first_global_vertex_of_proc[i] + m_partition_sizes[i]);
        }
        first_global_vertex_of_proc.push_back(num_vertices);

        for (unsigned int i = 0; i < num_vertices; ++i) {
            vertex_mapping[i] = vertex_mapping[i] % num_vertices +
                                first_global_vertex_of_proc[std::floor(
                                    vertex_mapping[i] / num_vertices)];
        }

        // Initialize graph
        g.start_initialization(first_global_vertex_of_proc);

        // Second: read edges into data structure or send edges to correct proc
        std::vector<std::vector<Edge>> msgs_for_proc(procs);

        for (unsigned int i = 0; i < m_adj.size(); ++i) {
            if (m_adj[i].size() == 0) {
                continue;
            }

            auto u = i + first_in;
            auto node = vertex_mapping[u];

            auto p = m_vertex_partition[u];

            if (p == proc_id) {
                for (auto v : m_adj[i]) {
                    auto e = Edge(node, vertex_mapping[v.node], v.weight);
                    auto neighbor_p = m_vertex_partition[v.node];
                    if (neighbor_p == proc_id) {
                        if (node < vertex_mapping[v.node])
                            g.add_local_edge(e);
                    } else {
                        g.add_cross_edge(e, neighbor_p);
                    }
                }
            } else {
                // send edges to correct proc !!
                auto dest_p_u = m_vertex_partition[u];

                for (auto v : m_adj[i]) {
                    msgs_for_proc[dest_p_u].push_back(
                        Edge(u, v.node, v.weight));
                }
            }
        }
        send_edges_to_correct_proc(msgs_for_proc, vertex_mapping, proc_id,
                                   procs, g);
    }

  private:
    unsigned int decide_proc_for_neighborhood(unsigned int local_node,
                                              unsigned int node,
                                              unsigned int proc_id,
                                              unsigned int procs) {

        std::fill(m_p_scores.begin(), m_p_scores.end(), 0);
        std::fill(m_edges_in_proc.begin(), m_edges_in_proc.end(), 0);

        for (auto edge : m_adj[local_node]) {
            auto p = m_vertex_partition[edge.node];

            m_edges_in_proc[p] += edge.weight;
        }

        for (unsigned p = 0; p < procs; ++p) {
            m_p_scores[p] =
                m_edges_in_proc[p] / m_sum_weights -
                m_node_weights[node] * fennel_alpha * fennel_gamma *
                    std::pow(m_partition_sizes[p], fennel_gamma - 1) /
                    m_sum_degrees;

            // LDG:
            // m_p_scores[p] =
            //     m_edges_in_proc[p] *
            //     (1 - (m_partition_sizes[p] * procs) / m_sum_weights);
        }

        return std::distance(
            m_p_scores.begin(),
            std::max_element(m_p_scores.begin(), m_p_scores.end()));
    }

    void calcPartitionSizes() {
        std::fill(m_partition_sizes.begin(), m_partition_sizes.end(), 0);

        for (unsigned int i = 0; i < m_vertex_partition.size(); ++i) {
            m_partition_sizes[m_vertex_partition[i]] += m_node_weights[i];
        }
    }

    double calcBalance(int procs) {
        double balance = 1;

        for (unsigned int i = 0; i < m_partition_sizes.size(); ++i) {
            auto balance_i = m_partition_sizes[i] / (m_sum_weights / procs);
            if (balance < balance_i) {
                balance = balance_i;
            }
        }
        return balance;
    }

    void
    send_edges_to_correct_proc(std::vector<std::vector<Edge>> &msgs_for_proc,
                               std::vector<NodeIDType> &vertex_mapping,
                               unsigned int proc_id, unsigned int procs,
                               Graph &g) {
        MPI_Datatype MPI_EDGE_TYPE;
        MPI_Type_contiguous(sizeof(Edge), MPI_BYTE, &MPI_EDGE_TYPE);
        MPI_Type_commit(&MPI_EDGE_TYPE);

        int tag = MPI_MSG_TAG::send_nodes_to_correct_proc;

        std::vector<MPI_Request> reqs(procs, MPI_REQUEST_NULL);

        // send cross edge_container to partners
        for (unsigned int p = 0; p < procs; p++) {
            if (p != proc_id) {
                if (!msgs_for_proc[p]
                         .empty()) { // send cross edge_container to p
                    MPI_Isend(&msgs_for_proc[p][0], msgs_for_proc[p].size(),
                              MPI_EDGE_TYPE, p, tag, MPI_COMM_WORLD, &reqs[p]);
                } else { // send empty message
                    MPI_Isend(0, 0, MPI_EDGE_TYPE, p, tag, MPI_COMM_WORLD,
                              &reqs[p]);
                }
            }
        }

        // receive cross edges
        for (unsigned int p = 0; p < procs; p++) {
            if (p != proc_id) {
                MPI_Status status;
                MPI_Probe(p, tag, MPI_COMM_WORLD, &status);

                int msg_count;
                MPI_Get_count(&status, MPI_EDGE_TYPE, &msg_count);

                std::vector<Edge> incoming_edges(msg_count);

                MPI_Recv(&incoming_edges[0], msg_count, MPI_EDGE_TYPE, p, tag,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (int m = 0; m < msg_count; m++) {
                    auto &e = incoming_edges[m];
                    auto p1 = m_vertex_partition[e.n1];
                    auto p2 = m_vertex_partition[e.n2];

                    e.n1 = vertex_mapping[e.n1];
                    e.n2 = vertex_mapping[e.n2];
                    if (p1 != proc_id) {
                        printf("[%d] Edge (%d,%d), with p1:%d and p2:%d\n",
                               proc_id, e.n1, e.n2, p1, p2);
                        throw std::invalid_argument("Received wrong edge.");
                    }
                    if (p2 == proc_id) {
                        if (e.n1 < e.n2)
                            g.add_local_edge(e);
                    } else {
                        g.add_cross_edge(incoming_edges[m], p2);
                    }
                }
            }
        }

        MPI_Waitall(procs, &reqs[0], MPI_STATUS_IGNORE);
    }

    void write_vertex_partitioning(unsigned int num_nodes) {

        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        if (rank == 0) {
            std::ofstream partitioning_file;
            partitioning_file.open("vertex_partitioning.txt");

            for (NodeIDType i = 0; i < num_nodes; ++i) {
                partitioning_file << m_vertex_partition[i] << "\n";
            }
            partitioning_file.close();
        }
    }
    void set_seed(int seed) override { m_seed = seed; }

    std::vector<std::vector<EndPoint>> m_adj;
    std::vector<unsigned int> m_vertex_partition;

    std::vector<unsigned int> m_buffer;

    std::vector<unsigned int> m_node_weights;
    std::vector<unsigned int> m_partition_sizes;

    double m_sum_weights;
    int m_sum_degrees;

    float fennel_alpha = 0;
    float fennel_gamma = 2.0f;
    int m_fennel_rounds = 5;
    float fennel_t_alpha = 2;

    int m_seed = 0;

    std::vector<int> m_edges_in_proc;
    std::vector<float> m_p_scores;
};

} // namespace graph_io

#endif
