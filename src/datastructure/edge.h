/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 * This file is partially based on LocalMaxMatching (MIT Licence, Marcel Birn)
 *   https://github.com/LocalMaxMatching/LocalMaxMatching
 *
 *****************************************************************************/

#ifndef _PARKEC_EDGE_H
#define _PARKEC_EDGE_H

#include "util/color_set.h"
#include <iostream>

namespace parkec {

template <typename WeightType, typename NodeIDType> struct Edge {
    static constexpr uint_fast8_t no_color = UINT_FAST8_MAX;
    // static constexpr uint_fast8_t not_colorable = UINT_FAST8_MAX - 1;

    NodeIDType n1;
    NodeIDType n2;

    WeightType weight;

    uint_fast8_t color;

    color_set palette;

    int tries;

    Edge() : n1(0), n2(0), weight(0), color(no_color) {}

    Edge(const NodeIDType n1, const NodeIDType n2, const WeightType weight)
        : n1(n1), n2(n2), weight(weight), color(no_color),
          palette(color_set(0)) {}

    Edge(const NodeIDType n1, const NodeIDType n2, const WeightType weight,
         int k)
        : n1(n1), n2(n2), weight(weight), color(no_color),
          palette(color_set(k)) {}
};

template <typename WeightType, typename NodeIDType> struct EdgeE {
    NodeIDType partner;
    WeightType weight;
    uint_fast8_t color = UINT_FAST8_MAX;

    EdgeE() : partner(-1), weight(0) {}

    EdgeE(const NodeIDType partner, const WeightType weight)
        : partner(partner), weight(weight) {}
};

// template <typename WeightType, typename NodeIDType> struct EdgeX {
//     NodeIDType partner;
//     WeightType weight;
//     uint_fast8_t color = UINT_FAST8_MAX;

//     EdgeX() : partner(-1), weight(0) {}

//     EdgeX(const NodeIDType partner, const WeightType weight)
//         : partner(partner), weight(weight) {}
// };
} // namespace parkec
#endif
