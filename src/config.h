/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 *****************************************************************************/

#ifndef CONFIG_H
#define CONFIG_H

#include <string>

// #define TESTING

typedef enum {
    MREP_LM,
    MREP_S,
    K_MM,
    K_CS,
} KEC_ALGO;

typedef enum {
    NO_PARTITIONING,
    METIS,
    METIS_UNWEIGHTED,
    METIS_NODE_WEIGHTS_UNWEIGHTED,
    KAHIP_UNWEIGHTED,
    RANDOM_PARTITIONING,
    DIST_STREAM
} PARTITIONING;

struct Config {
  public:
    KEC_ALGO algo;
    unsigned int k;
    PARTITIONING partitioning;
    int seed;
    double ccm_epsilon = 0.1;

#ifdef TESTING
    bool test_param_1;
    bool test_param_2;
#endif
};

#endif