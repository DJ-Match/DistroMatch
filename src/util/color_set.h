/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 * This file is based on DyDJ Match (Copyright (C) 2022-2023 : Kathrin Hanauer,
 * Lara Ost)
 *   https://github.com/DJ-Match/DyDJ-Match/blob/main/src/tools/color_set.h
 *
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <assert.h>
#include <limits>
#include <ostream>
#include <random>

// Holds the available colors of a vertex
struct color_set {
    using size_type = uint_fast8_t;
    static constexpr size_type npos = -1;
    static constexpr size_type not_colorable = -2;
    using bit_type = unsigned long;
    // using bit_type = unsigned __int128;

  private:
    static constexpr bit_type one = 1u;

    bit_type bits;
    bit_type all_bits;
    size_type bit_size;

    inline color_set(bit_type bits, bit_type all, size_type size)
        : bits(bits), all_bits(all), bit_size(size) {}

    friend color_set operator&(const color_set &a, const color_set &b);

  public:
    color_set() : bits(0), all_bits(0), bit_size(0) {}

    // Create a new `color_set`. By default, all colors are included.
    color_set(size_type size)
        : all_bits(
              size == sizeof(bit_type) * 8
                  ? std::numeric_limits<bit_type>::max()
                  : (one << size) -
                        1), // 2^bit_size - 1 -> bit_size many bits are set to 1
          bit_size(size) {

        if (size > sizeof(bit_type) * 8) {

            throw std::invalid_argument("Number of colors is too large.");

            // std::format("Cannot handle color sets with > {} colors.",
            // sizeof(bit_type) * 8) -> gcc 13+
        }

        bits = all_bits;
    }

    inline size_type find_first() const {
        if (bits == 0) {
            return npos;
        } else {
            return lowestBit(bits);
        }
    }

    inline size_type find_next(size_type pos) const {
        auto shifted = bits >> pos;
        if (shifted == 0) {
            return npos;
        } else {
            return lowestBit(shifted) + pos;
        }
    }

    inline size_type find_random(std::mt19937 &rng) const {
        if (bits == 0)
            return npos; // No bits set

        // Count the number of active bits
        size_type active_count = __builtin_popcountl(bits);

        // Generate a random active bit index
        std::uniform_int_distribution<size_type> dist(0, active_count - 1);
        size_type random_index = dist(rng);

        // Locate the `random_index`-th active bit efficiently
        bit_type temp = bits;
        while (random_index > 0) {
            // Remove the lowest set bit and decrement the index
            temp &= temp - 1;
            --random_index;
        }

        // The remaining lowest set bit is our target
        return __builtin_ctzl(temp);
    }

    inline color_set find_random(std::mt19937 &rng, size_type amount) const {
        if (bits == 0 || amount == 0)
            return color_set(0); // No bits set or no numbers requested

        // Count the number of active bits
        size_type active_count = __builtin_popcountl(bits);

        // Limit the number of random bits to at most `active_count`
        amount = std::min(amount, active_count);

        // Create a new color_set to hold the random selection
        color_set random_set(bit_size);
        random_set.bits = 0; // Clear all bits initially

        // Extract active bit positions into a vector
        std::vector<size_type> active_positions;
        active_positions.reserve(active_count);

        bit_type temp = bits;
        while (temp != 0) {
            size_type pos = __builtin_ctzl(temp);
            active_positions.push_back(pos);
            temp &= temp - 1; // Remove the lowest set bit
        }

        // Shuffle the active positions and pick the first `amount` positions
        std::shuffle(active_positions.begin(), active_positions.end(), rng);

        for (size_type i = 0; i < amount; ++i) {
            random_set.bits |= (bit_type(1) << active_positions[i]);
        }

        return random_set;
    }

    inline void flip() {
        // Flip all bits, then set unused bits to 0
        bits = (~bits) & all_bits;
    }

    inline void set() { bits = all_bits; }

    inline void setOn(size_type i) { bits = bits | (one << i); }

    inline void setOff(size_type i) { bits = bits & ~(one << i); }

    inline void unset_all() { bits = 0; }

    inline bool none() const { return bits == 0; }

    inline bool any() const { return !none(); }

    inline bool all() const { return bits == all_bits; }

    inline void make_union(const color_set &a) { bits = bits | a.bits; }
    inline void make_intersect(const color_set &a) { bits = bits & a.bits; }

    inline void make_intersect(const bit_type &a_bits) { bits = bits & a_bits; }

    inline size_type count() const {
        size_type count = 0;
        for (size_type i = 0; i < bit_size; ++i) {
            count += ((bits & (one << i)) != 0);
        }
        return count;
    }

    inline bit_type get_bits() const { return bits; }

    inline size_type size() const { return bit_size; }

    static inline color_set common_colors(const color_set &a,
                                          const color_set &b) {
        return (a & b);
    }

    static inline color_set common_colors(const color_set &a,
                                          const bit_type &b_colors) {
        return {a.bits & b_colors, a.all_bits, a.bit_size};
    }

    inline bool operator[](size_type i) const { return bits & (one << i); }

    friend std::ostream &operator<<(std::ostream &stream, const color_set &c);

  private:
    static size_type lowestBit(bit_type b) { return __builtin_ctzl(b); }
};

inline color_set operator&(const color_set &a, const color_set &b) {
    assert(a.size() == b.size());
    return {a.bits & b.bits, a.all_bits, a.bit_size};
}

inline std::ostream &operator<<(std::ostream &stream, const color_set &c) {
    for (color_set::size_type i = c.bit_size; i > 0; --i) {
        stream << (int)(c[i - 1]);
    }
    return stream;
}
