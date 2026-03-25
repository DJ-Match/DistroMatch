/******************************************************************************
 * Source of DistroMatch - Distributed Disjoint Weighted Matchings.
 *
 * Contact information:
 *   https://github.com/DJ-Match/DistroMatch
 *
 * This file is based on LocalMaxMatching (MIT Licence, Marcel Birn)
 *   https://github.com/LocalMaxMatching/LocalMaxMatching
 *****************************************************************************/

#ifndef HASH_FUNCTIONS_H_
#define HASH_FUNCTIONS_H_

namespace parkec {

/// http://www.concentric.net/~ttwang/tech/inthash.htm 32 bit Mix Functions
// is a reversible function (computed all possible values and counted them
// -> no duplicates)
inline unsigned int hash(int val) {
    int key = val;
    key = ~key + (key << 15); // key = (key << 15) - key - 1;
    key = key ^ ((unsigned int)key >> 12);
    key = key + (key << 2);
    key = key ^ ((unsigned int)key >> 4);
    key = key * 2057; // key = (key + (key << 3)) + (key << 11);
    key = key ^ ((unsigned int)key >> 16);
    return key;
}

inline unsigned int hash(unsigned int val) { return hash((int)val); }

// added this wrapper just for compatible reasons - was lazy
inline unsigned int hash2(int val) { return hash(val); }

/// http://www.concentric.net/~ttwang/tech/inthash.htm 64 bit Mix Functions
inline long hash(long key) {
    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ ((unsigned long)key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ ((unsigned long)key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ ((unsigned long)key >> 28);
    key = key + (key << 31);
    return key;
}

inline long hash(unsigned long key) { return hash((long)key); }

unsigned int r;

unsigned int simple_random_hash(unsigned int val) {
    unsigned long h = val * r;
    return (unsigned int)h;
}

} // namespace parkec

#endif
