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

#include <array>
#include <assert.h>
#include <limits>
#include <ostream>
#include <stdexcept>

struct color_set {
    using size_type = uint_fast8_t;
    static constexpr size_type npos = -1;
    static constexpr size_type not_colorable = -2;
    static constexpr size_type num_segments = 2;
    static constexpr uint64_t one = 1u;
    using bit_type =
        std::array<uint64_t, num_segments>; // Fixed-size array for 2 segments

  private:
    bit_type bits{};
    bit_type all_bits{};
    size_type bit_size;

    inline color_set(bit_type bits, bit_type all, size_type size)
        : bits(bits), all_bits(all), bit_size(size) {}

    friend color_set operator&(const color_set &a, const color_set &b);

  public:
    color_set() : bits{}, all_bits{}, bit_size(0) {}

    color_set(size_type size) : bits{}, all_bits{}, bit_size(size) {
        if (size > num_segments * 64) {
            throw std::invalid_argument("Number of colors is too large.");
        }
        // Calculate full segments
        size_type full_segments = size / 64;
        size_type remaining_bits = size % 64;

        // Set each full segment to all 1s
        for (size_type i = 0; i < full_segments; ++i) {
            all_bits[i] = UINT64_MAX;
        }

        if (remaining_bits > 0) {
            all_bits[full_segments] = (one << remaining_bits) - 1;
        }

        // Initialize bits with all_bits to include all colors by default
        bits = all_bits;
    }

    inline size_type find_first() const {
        if (bits[0] != 0) {
            return lowestBit(bits[0]);
        } else if (bits[1] != 0) {
            return lowestBit(bits[1]) + 64;
        }
        return npos;
    }

    inline size_type find_next(size_type pos) const {
        size_type segment = pos / 64;
        size_type offset = pos % 64;
        uint64_t shifted = bits[segment] >> offset;
        if (shifted != 0) {
            return lowestBit(shifted) + pos;
        }
        if (segment + 1 < num_segments && bits[segment + 1] != 0) {
            return lowestBit(bits[segment + 1]) + (segment + 1) * 64;
        }
        return npos;
    }

    inline void flip() {
        bits[0] = ~bits[0] & all_bits[0];
        bits[1] = ~bits[1] & all_bits[1];
    }

    inline void set() { bits = all_bits; }

    inline void setOn(size_type i) { bits[i / 64] |= (one << (i % 64)); }

    inline void setOff(size_type i) { bits[i / 64] &= ~(one << (i % 64)); }

    inline bool none() const { return bits[0] == 0 && bits[1] == 0; }

    inline bool any() const { return bits[0] != 0 || bits[1] != 0; }

    inline bool all() const {
        return bits[0] == all_bits[0] && bits[1] == all_bits[1];
    }

    inline size_type count() const {
        return __builtin_popcountll(bits[0]) + __builtin_popcountll(bits[1]);
    }

    inline size_type size() const { return bit_size; }

    static inline color_set common_colors(const color_set &a,
                                          const color_set &b) {
        return a & b;
    }

    inline bool operator[](size_type i) const {
        return bits[i / 64] & (one << (i % 64));
    }

    friend std::ostream &operator<<(std::ostream &stream, const color_set &c);

  private:
    static size_type lowestBit(uint64_t b) { return __builtin_ctzl(b); }
};

inline color_set operator&(const color_set &a, const color_set &b) {
    assert(a.size() == b.size());
    color_set::bit_type new_bits = {a.bits[0] & b.bits[0],
                                    a.bits[1] & b.bits[1]};
    return {std::move(new_bits), a.all_bits, a.bit_size};
}

inline std::ostream &operator<<(std::ostream &stream, const color_set &c) {
    for (color_set::size_type i = c.bit_size - 1; i >= 0; --i) {
        stream << (int)(c[i / 64][i % 64]);
    }
    return stream;
}
