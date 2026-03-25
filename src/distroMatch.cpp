/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#include <algorithm/k_edge_coloring.h>
#include <config.h>
#include <datastructure/graph.h>
#include <io/parallel_range_reader.h>
#include <io/parse_parameter.h>
#include <select_algorithm.h>
#include <util/memory.h>
#include <util/mpi_system_support.h>

// #define LOGGING
// #define WRITETOFILE
// #define LOGGING_V

#ifdef LOGGING_V
unsigned int watched_id = 386830;
#endif

int main(int argc, char **argv) {
    // {
    //     int i = 0;
    //     char hostname[256];
    //     gethostname(hostname, sizeof(hostname));
    //     printf("PID %d on %s ready for attach\n", getpid(), hostname);
    //     fflush(stdout);
    //     while (0 == i)
    //         sleep(5);
    // }

    MPI_Init(&argc, &argv);

    auto baseline = getPeakRSS();
    int proc_id, procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

#ifdef WRITETOFILE
    std::string out = "log/output" + std::to_string(proc_id) + ".ansi";
    freopen(out.c_str(), "w", stdout);
#endif
    // warm up the network - necessary for the IC1
    parkec::warm_up_network(MPI_COMM_WORLD);

    bool is_root = proc_id == 0;

    // parse input parameters
    Config config;

    std::string infile;
    std::string outfile;

    parse_parameters(argc, argv, infile, outfile, config);
    if (is_root) {
        printf("Reading graph from file %s\n", infile.c_str());
    }

    typedef parkec::IGraph<> Graph;

    typedef typename Graph::Edge Edge;
    typedef typename Graph::EdgeIDType EdgeIDType;
    typedef typename Graph::NodeIDType NodeIDType;
    typedef typename Graph::WeightType WeightType;

    unsigned int k = config.k;
    NodeIDType num_vertices = 0;

    double start_read_file = MPI_Wtime();

    std::unique_ptr<Graph> g = parkec::get_graph(config);

    graph_io::read_parallel_range<WeightType, NodeIDType, EdgeIDType, Edge,
                                  Graph>(infile, *g, num_vertices, config);

    // g->print();

    MPI_Barrier(MPI_COMM_WORLD);
    double duration_read_file = MPI_Wtime() - start_read_file;

    if (is_root) {
        printf("using %d processes\nread time: %.6f s\ncomputing %d disjoint "
               "matchings\n",
               procs, duration_read_file, k);
    }

    unsigned long global_local_edge_count, global_cross_edge_count;
    unsigned long local_edge_count = g->getNumLocalEdges();
    unsigned long local_cross_edge_count = g->getNumCrossEdges();

    MPI_Reduce(&local_edge_count, &global_local_edge_count, 1,
               MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_cross_edge_count, &global_cross_edge_count, 1,
               MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (is_root) {
        printf("vertices: %u  local_edges: %lu cross_edges: %lu\n",
               num_vertices, global_local_edge_count,
               global_cross_edge_count / 2);
    }

    // synchronize before computing the matching, other process might still
    // read the input graphs
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    int total_global_matching_size = 0;
    double total_global_matching_weight = 0.0;
    int total_number_of_rounds = 0;

    std::vector<std::vector<Edge>> matchings;
    std::shared_ptr<parkec::KEdgeColoring<Graph>> algorithm =
        parkec::get_algorithm<Graph>(config);

    algorithm->compute_k_edge_coloring(matchings, total_global_matching_size,
                                       total_number_of_rounds,
                                       total_global_matching_weight, *g, k);

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    double elapsed_time_is_s(end_time - start_time);

    std::vector<int> matching_sizes(k, 0);
    std::vector<double> matching_weight(k, 0);

    for (unsigned int c = 0; c < matchings.size(); ++c) {
        for (Edge &e : matchings[c]) {
            matching_weight[c] += e.weight;
        }

        matching_sizes[c] = matchings[c].size();
    }
    if (proc_id == 0) {
        MPI_Reduce(MPI_IN_PLACE, matching_sizes.data(), k, MPI_INT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, matching_weight.data(), k, MPI_DOUBLE, MPI_SUM,
                   0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(matching_sizes.data(), nullptr, k, MPI_INT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
        MPI_Reduce(matching_weight.data(), nullptr, k, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);
    }

    if (proc_id == 0) {
        for (unsigned int i = 0; i < k; ++i) {
            printf("%d: %d, %.4f\n", i, matching_sizes[i], matching_weight[i]);
        }
    }

    auto memory = getPeakRSS() - baseline;
    auto total_memory = 0;

    MPI_Reduce(&memory, &total_memory, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (is_root) {
        printf("running_time: %.6f s\n", elapsed_time_is_s);
        printf("total_size: %d\n", total_global_matching_size);
        printf("total_weight: %.2f\n", total_global_matching_weight);
        printf("num_rounds: %d\n", total_number_of_rounds);
        printf("mem: %d\n", total_memory);
    }

    if (!g->valid_k_matchings(matchings, true)) {
        if (is_root) {
            printf("Invalid or not maximal matching!\n");
        }
    }
#ifdef LOGGING
#ifdef LOGGING_V
    printf("(watched: v%d) ", watched_id);
#endif
    printf("[%d] ", proc_id);
    for (unsigned int c = 0; c < matchings.size(); ++c) {
        printf("c%d: ", c);

        for (Edge &e : matchings[c]) {
#ifdef LOGGING_V
            if (e.n1 == watched_id || e.n2 == watched_id)
#endif
                printf("(v%d,v%d) ", e.n1, e.n2);
        }
        printf("|");
    }
    printf("\n");
#endif

    // printf("[%d] ", proc_id);
    // for (unsigned int c = 0; c < matchings.size(); ++c) {
    //     printf("c%d: ", c);
    //     for (Edge &e : matchings[c]) {
    //         printf("(v%d,v%d,w%.2f) ", e.n1, e.n2, e.weight);
    //     }
    //     printf("|");
    // }
    // printf("\n");

    MPI_Finalize();

    return 0;
}