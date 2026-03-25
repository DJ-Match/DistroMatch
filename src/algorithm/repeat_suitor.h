/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 * The code for a single matching is based on the description in
 *
 * F. Manne and M. Halappanavar “New Effective Multithreaded Matching
 * Algorithms,” in 2014 IEEE 28th International Parallel and Distributed
 * Processing Symposium, Phoenix, AZ, USA, May 19-23, 2014, IEEE Computer
 * Society, 2014, pp. 519–528. doi: 10.1109/IPDPS.2014.61.
 *
 *****************************************************************************/

#ifndef REPEAT_SUITOR_H
#define REPEAT_SUITOR_H

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
#include <algorithm/k_match.h>

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
class RepeatSuitor : public KEdgeColoring<GraphBase> {
    typedef typename parkec::DisAdjGraph<> Graph;

    typedef typename Graph::WeightType WeightType;
    typedef typename Graph::NodeIDType NodeIDType;
    typedef typename Graph::EdgeIDType EdgeIDType;
    typedef typename Graph::EdgeE EdgeE;
    typedef typename Graph::Edge Edge;

  public:
    RepeatSuitor() {}

    void compute_k_edge_coloring(std::vector<std::vector<Edge>> &matchings,
                                 int &total_global_matching_size,
                                 int &total_number_of_rounds,
                                 double &total_global_matching_weight,
                                 GraphBase &g_ref, unsigned int k) {

        Graph &g = static_cast<Graph &>(g_ref);
        auto k_dm = std::make_unique<GreedyKMatching<Graph>>();
        k_dm->set_color();

        matchings.resize(k, std::vector<Edge>(0));

        // Run 1-DM
        // std::vector<std::vector<Edge>> matching_i;

        total_global_matching_size = 0;
        total_number_of_rounds = 0;
        total_global_matching_weight = 0;

        int one_global_matching_size;
        int one_number_of_rounds;
        double one_global_matching_weight;

        for (unsigned int i = 0; i < k; ++i) {
            Wrapper wrapped(matchings[i]);

            one_global_matching_size = 0;
            one_number_of_rounds = 0;
            one_global_matching_weight = 0;

            // matching_i.clear();
            k_dm->compute_k_matchings(
                reinterpret_cast<std::vector<std::vector<Edge>> &>(wrapped),
                one_global_matching_size, one_number_of_rounds,
                one_global_matching_weight, g, 1);

            // Remove colored edges
            g.reduceGraph();
            g.activate_all_procs();

            k_dm->reset();

            // matchings[i] = matching_i[0];

            total_global_matching_size += one_global_matching_size;
            total_number_of_rounds += one_number_of_rounds;
            total_global_matching_weight += one_global_matching_weight;
        }
        // throw std::logic_error("DEBUG");
    }

    // Helper to avid copying
    class Wrapper {
      public:
        explicit Wrapper(std::vector<Edge> &vec) : wrapped_vec(vec) {}

        std::vector<Edge> &operator[](size_t index) {
            if (index != 0) {
                throw std::out_of_range("Wrapper: Index out of range");
            }
            return wrapped_vec;
        }

        const std::vector<Edge> &operator[](size_t index) const {
            if (index != 0) {
                throw std::out_of_range("Wrapper: Index out of range");
            }
            return wrapped_vec;
        }

        // Pretend to have a size of 1.
        size_t size() const { return 1; }

      private:
        std::vector<Edge> &wrapped_vec;
    };
};

} // namespace parkec

#endif
