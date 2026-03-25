/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef SELECT_ALGORITHM_H
#define SELECT_ALGORITHM_H

#include <memory>

#include "config.h"
#include <algorithm/ccm_conflicts.h>
#include <algorithm/k_edge_coloring.h>
#include <algorithm/k_match.h>
#include <algorithm/match_repeat.h>
#include <algorithm/repeat_suitor.h>
#include <datastructure/dis_adj_graph.h>
#include <datastructure/graph.h>
#include <datastructure/min_parallel_edge_graph.h>
#include <datastructure/parallel_edge_graph.h>

namespace parkec {

template <class Graph>
std::unique_ptr<KEdgeColoring<Graph>> get_algorithm(Config &config) {

    std::unique_ptr<KEdgeColoring<Graph>> algo;
    switch (config.algo) {
    case MREP_LM:
        algo = std::make_unique<GreedyMatchingRepeat<Graph>>();
        break;
    case K_CS:
        algo = std::make_unique<GreedyKMatching<Graph>>();
        break;
    case MREP_S:
        algo = std::make_unique<RepeatSuitor<Graph>>();
        break;
    case K_MM:
        algo = std::make_unique<DirectCCM<Graph>>(config.ccm_epsilon);
        break;
    default:
        exit(0);
        break;
    }

#ifdef TESTING
    algo->set_test_param_1(config.test_param_1);
    algo->set_test_param_2(config.test_param_2);
#endif

    return algo;
}

std::unique_ptr<IGraph<>> get_graph(Config &config) {

    switch (config.algo) {
    case MREP_LM:
        return std::make_unique<parkec::ParallelEdgeGraph<>>();
    case K_MM:
        return std::make_unique<parkec::MinParallelEdgeGraph<>>();
        break;
    case K_CS:
    case MREP_S:
        return std::make_unique<parkec::DisAdjGraph<>>();
    default:
        exit(0);
        break;
    }
}

}; // namespace parkec

#endif /* end of include guard: FACTORY_7NX3X */