/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef VERTEX_MAPPING_H
#define VERTEX_MAPPING_H

#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

namespace parkec {

template <typename NodeIDType>
void generate_vertex_mapping(std::vector<NodeIDType> &vertex_mapping,
                             int num_vertices, unsigned int seed) {

    vertex_mapping.resize(num_vertices);

    std::iota(std::begin(vertex_mapping), std::end(vertex_mapping), 0);

    // Mersenne Twister RNG with the provided seed
    std::mt19937 rng(seed);

    std::shuffle(vertex_mapping.begin(), vertex_mapping.end(), rng);
}

void generate_initial_partitioning(std::vector<unsigned int> &partitioning,
                                   unsigned int element_cout,
                                   unsigned int max_num, unsigned int start,
                                   unsigned int seed) {

    partitioning.reserve(element_cout);

    for (unsigned int i = start; i < element_cout + start; ++i) {
        partitioning.push_back(i % max_num);
    }

    // Mersenne Twister RNG with the provided seed
    std::mt19937 rng(seed);
    std::shuffle(partitioning.begin(), partitioning.end(), rng);
}
} // namespace parkec

#endif // VERTEX_MAPPING_H
