/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef PARALLEL_RANGE_READER_H
#define PARALLEL_RANGE_READER_H

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <mpi.h>

#include <config.h>
#include <datastructure/parallel_edge_graph.h>
#include <io/bin_reader.h>
#include <io/metis_reader.h>
#include <io/mtx_fast_partitioning_reader.h>
#include <io/mtx_reader.h>
#include <io/reader.h>
#include <util/mpi_tags.h>
#include <util/random.h>

namespace graph_io {

template <class Graph>
std::unique_ptr<GraphReader<Graph>> get_reader(const char *format) {

    std::unique_ptr<GraphReader<Graph>> reader;
    if (strcmp(".graph", format) == 0) {
        reader = std::make_unique<MetisGraphReader<Graph>>();
    } else if (strcmp(".mtx", format) == 0) {
        reader = std::make_unique<MtxGraphReader<Graph>>();
    } else if (strcmp(".bin", format) == 0) {
        reader = std::make_unique<BinGraphReader<Graph>>();
    } else {
        throw std::invalid_argument("Graph input format not known. Supported "
                                    "are METIS (\".graph\") and "
                                    "Matrix Market (\".mtx\").");
    }

    return reader;
}

template <typename WeightType, typename NodeIDType, typename EdgeIDType,
          typename Edge, class Graph>
void read_parallel_range(const std::string &filename, Graph &g, NodeIDType &n,
                         const Config &config) {

    std::filesystem::path fsPath(filename);
    std::string directory = fsPath.parent_path().string();
    std::string graph_name = fsPath.stem().string();
    std::string graph_ext = fsPath.extension().string();

    std::ifstream in;

    if (strcmp(".bin", graph_ext.c_str()) == 0) {
        in.open(filename.c_str(), std::ios::binary);
    } else {
        in.open(filename.c_str()); // Text mode
    }

    int procs = 1;
    int proc_id = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    std::unique_ptr<GraphReader<Graph>> reader;

    if (config.partitioning == DIST_STREAM) {
        if (strcmp(".mtx", graph_ext.c_str()) != 0) {
            throw std::invalid_argument(
                "With distributed partitioning only Matrix Market (\".mtx\") "
                "graph format supported.");
        }
        reader = std::make_unique<MtxPartitionGraphReader<Graph>>();

    } else {
        reader = get_reader<Graph>(graph_ext.c_str());
    }

    reader->read_header(in, n);

    std::string partition_filename;

    switch (config.partitioning) {
    case METIS:
        partition_filename = directory + "/part/" + graph_name +
                             ".graph.part." + std::to_string(procs) +
                             ".metis.s" + std::to_string(config.seed);
        break;
    case METIS_UNWEIGHTED:
        partition_filename = directory + "/part/" + graph_name +
                             "_u.graph.part." + std::to_string(procs) +
                             ".metis.s" + std::to_string(config.seed);
        break;
    case METIS_NODE_WEIGHTS_UNWEIGHTED:
        partition_filename = directory + "/part/" + graph_name +
                             "_uw.graph.part." + std::to_string(procs) +
                             ".metis.s" + std::to_string(config.seed);
        break;
    case KAHIP_UNWEIGHTED:
        partition_filename = directory + "/part/" + graph_name +
                             "_u.graph.part." + std::to_string(procs) +
                             ".kahip.s" + std::to_string(config.seed);
        break;
    default:
        break;
    }

    std::vector<NodeIDType> vertex_mapping;
    std::vector<NodeIDType> num_vertices_of_proc;
    std::vector<NodeIDType> first_global_vertex_of_proc;
    NodeIDType base_size = n / procs;

    first_global_vertex_of_proc.reserve(procs + 1);
    first_global_vertex_of_proc.push_back(0);

    bool partition_file_exists = std::filesystem::exists(partition_filename);

    bool read_partition =
        config.partitioning == METIS ||
        config.partitioning == METIS_UNWEIGHTED ||
        config.partitioning == METIS_NODE_WEIGHTS_UNWEIGHTED ||
        config.partitioning == KAHIP_UNWEIGHTED;

    bool use_partitioning = config.partitioning != NO_PARTITIONING;

    if (read_partition && !partition_file_exists) {
        if (proc_id == 0)
            printf("Could not find partition file %s, continuing without ...\n",
                   partition_filename.c_str());
        use_partitioning = false;
    }

    if (read_partition && partition_file_exists) {
        std::ifstream in_part(partition_filename.c_str());
        if (!in_part) {
            if (proc_id == 0)
                std::cerr << "Error opening " << filename << std::endl;
            exit(1);
        }
        if (proc_id == 0)
            printf("Read partition from %s\n", partition_filename.c_str());
        vertex_mapping.reserve(n);
        num_vertices_of_proc.resize(procs, 0);

        std::string part_line;
        while (std::getline(in_part, part_line)) {
            std::stringstream strs(part_line);
            int proc_of_vertex;
            strs >> proc_of_vertex;
            vertex_mapping.push_back(num_vertices_of_proc[proc_of_vertex] +
                                     proc_of_vertex * n);
            ++num_vertices_of_proc[proc_of_vertex];
        }

        for (int i = 0; i < procs - 1; ++i) {
            first_global_vertex_of_proc.push_back(
                first_global_vertex_of_proc[i] + num_vertices_of_proc[i]);
        }

        for (unsigned int i = 0; i < n; ++i) {
            vertex_mapping[i] =
                vertex_mapping[i] % n +
                first_global_vertex_of_proc[std::floor(vertex_mapping[i] / n)];
        }
    } else {
        if (config.partitioning == RANDOM_PARTITIONING) {
            if (proc_id == 0)
                printf("Generate random partitioning with seed %d\n",
                       config.seed);
            double start_time = MPI_Wtime();

            parkec::generate_vertex_mapping<NodeIDType>(vertex_mapping, n,
                                                        config.seed);

            MPI_Barrier(MPI_COMM_WORLD);
            double end_time = MPI_Wtime();
            double elapsed_time_is_s(end_time - start_time);
            if (proc_id == 0)
                printf("Random generation took %.5f s\n", elapsed_time_is_s);
        }

        if (config.partitioning != DIST_STREAM) {

            for (int p = 0; p < procs - 1; p++) {
                NodeIDType range_size = base_size + (p < (int)(n % procs));
                first_global_vertex_of_proc.push_back(
                    first_global_vertex_of_proc[p] + range_size);
            }
        }
    }

    if (config.partitioning != DIST_STREAM) {
        first_global_vertex_of_proc.push_back(n);

        g.start_initialization(first_global_vertex_of_proc);
    }

    if (config.partitioning == DIST_STREAM) {
        reader->set_seed(config.seed);
    }

    reader->read_graph(in, g, use_partitioning, n, proc_id, procs,
                       first_global_vertex_of_proc, vertex_mapping);

    in.close();

    g.finalize_initialization();
}

} // namespace graph_io

#endif
