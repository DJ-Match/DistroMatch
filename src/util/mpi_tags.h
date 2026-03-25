/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef MPI_MSG_TAG_H
#define MPI_MSG_TAG_H

struct MPI_MSG_TAG {
  public:
    // test_for_matching_tag
    static const unsigned int test_matching1 = 1000;
    static const unsigned int test_matching2 = 2000;

    // send_cross_edges_to_neighbors
    static const unsigned int send_cross_edges_to_neighbors = 7;
    static const unsigned int send_nodes_to_correct_proc = 8;

    // GreedyMatchingRepeat/ ColoredCandidateMatching Tags
    static const unsigned int candidates = 10000000;
    static const unsigned int candidates_rep = 110000;
    static const unsigned int colored_edges = 120000;
    static const unsigned int cross_degrees = 130000;
    static const unsigned int used_colors = 140000;
    static const unsigned int ghost_vertices = 200000;
    static const unsigned int vertex_weights = 600000;
    static const unsigned int sanity = 300000;

    static const unsigned int procs_of_vertex = 400000;
    static const unsigned int colors_of_vertex = 500000;

    // K-Split Tags
    static const unsigned int req_color_cross_edge = 20;
    static const unsigned int rep_color_cross_edge = 21;
    static const unsigned int inform_color_cross_edge = 22;
    static const unsigned int inform_adj_ghost = 23;

    static int get_procs_of_vertex_tag(unsigned int phase) {
        return procs_of_vertex + phase;
    }
    static int get_colors_of_vertex_tag(unsigned int phase) {
        return colors_of_vertex + phase;
    }

    static int get_vertex_weights_tag(unsigned int phase, unsigned int round) {
        return vertex_weights + phase * 10000 + round;
    }

    static int get_candidates_tag(unsigned int phase, unsigned int round,
                                  unsigned int cr) {
        return candidates * (phase + 1) + round * 100000 + cr;
    }

    static int get_candidates_tag(unsigned int round, unsigned int cr) {
        return candidates + round * 100000 + cr;
    }

    static int get_candidates_tag(unsigned int round) {
        return candidates + round;
    }
    static int get_candidates_rep_tag(unsigned int round) {
        return candidates_rep + round;
    }
    static int get_colored_edges(unsigned int round) {
        return colored_edges + round;
    }

    static int get_used_colors_tag(unsigned int round) {
        return used_colors + round;
    }

    static int get_ghost_vertices_tag(unsigned int round) {
        return ghost_vertices + round;
    }
    static int get_ghost_vertices_tag(unsigned int phase, unsigned int round) {
        return ghost_vertices + phase * 10000 + round;
    }
    static int sanity_check(unsigned int round) { return sanity + round; }
};

#endif
