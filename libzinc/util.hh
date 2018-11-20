#pragma once

#include <cstdint>
#include <cassert>

#include <libzinc/encoding.hh>
#include <immintrin.h>

namespace aether {

namespace morton {

static uint64_t fast_log2(const uint64_t x) {
    assert(x != 0);
    return sizeof(x) * 8 - 1 - __builtin_clzll(x);
}

template<uint32_t Dimension, uint32_t BitsPerDimension>
// This gives the maximum alignment level possible for a given number, 
// e.g. given 2 it returns 0, given 8 it returns 1, given 0 it returns 32.
static uint64_t get_max_align_level(const uint64_t code){
    return code == 0 ? BitsPerDimension : __builtin_ctzll(code)/ Dimension;
}

template<uint32_t Dimension>
// Returns the level size required for two points to be in the same cell
static uint64_t get_unifying_level(const uint64_t min, const uint64_t max){
    assert(max >= min);
    // the gets the maximum allowed for this range, max+1
    return max == min ? 0 : (fast_log2(max ^ min)/Dimension) +1;
}

template<uint32_t Dimension>
// Returns the  morton code for a given level, starting from 0;
static uint64_t get_morton_code(const uint64_t level) {
    return (1LLU << level * Dimension) -1;
}

template<uint32_t Dimension>
// given a morton value it will return the morton aligned value of size level that contains it.
// e.g. given 14,1 it will give 12. given 14,2 it will give 0.
static uint64_t get_parent_morton_aligned(const uint64_t code, uint32_t level) {
    return (code >> (Dimension * level)) << (Dimension * level);
}

template<uint32_t Dimension, uint32_t BitsPerDimension>
// given a range this will return the next morton aligned value
static uint64_t get_align_max(const uint64_t min, const uint64_t max) {
    assert(max >= min);
    if (max == min) return min;
    // this gets the maximum allowed for this value
    uint64_t align_max = (min != 0) ? 
          min + get_morton_code<Dimension>(get_max_align_level<Dimension,BitsPerDimension>(min))
        : std::numeric_limits<uint64_t>::max();
    uint64_t max_align = min + get_morton_code<Dimension>((fast_log2(max+1 - min)/Dimension));
    return std::min(align_max,max_align);
}


static const uint64_t __morton_2_x_mask = 0x5555555555555555;
static const uint64_t __morton_2_y_mask = 0xaaaaaaaaaaaaaaaa;

static const uint64_t __morton_3_x_mask = 0x1249249249249249;
static const uint64_t __morton_3_y_mask = 0x2492492492492492;
static const uint64_t __morton_3_z_mask = 0x4924924924924924;

#ifdef __BMI2__
template<typename Integer>
static inline constexpr Integer pdep(Integer v, Integer mask) {
    static_assert(sizeof(Integer) == 4 || sizeof(Integer) == 8);
    if constexpr (sizeof(Integer) == 4)
        return _pdep_u32(v, mask);
    else if constexpr (sizeof(Integer) == 8)
        return _pdep_u64(v, mask);
}

template<typename Integer>
static inline constexpr Integer pext(Integer v, Integer mask) {
    static_assert(sizeof(Integer) == 4 || sizeof(Integer) == 8);
    if constexpr (sizeof(Integer) == 4)
        return _pext_u32(v, mask);
    else if constexpr (sizeof(Integer) == 8)
        return _pext_u64(v, mask);
}

template<typename Integer>
static inline constexpr Integer expand_bits_2(Integer v) {
    return pdep(v, __morton_2_x_mask);
}

template<typename Integer>
static inline constexpr Integer compact_bits_2(Integer v) {
    return pext(v, __morton_2_x_mask);
}

template<typename Integer>
static inline constexpr Integer expand_bits_3(Integer v) {
    return pdep(v, __morton_3_x_mask);
}

template<typename Integer>
static inline constexpr Integer compact_bits_3(Integer v) {
    return pext(v, __morton_3_x_mask);
}

#else
template<typename Integer>
static_assert(sizeof(Integer) == 8, "non BMI2 expand/compact_bits_2 are not yet implemented");
static inline constexpr Integer expand_bits_2(Integer x) {
    x = (x ^ (x << 16)) & 0x0000ffff0000ffff;
    x = (x ^ (x << 8))  & 0x00ff00ff00ff00ff;
    x = (x ^ (x << 4))  & 0x0f0f0f0f0f0f0f0f;
    x = (x ^ (x << 2))  & 0x3333333333333333;
    x = (x ^ (x << 1))  & 0x5555555555555555;
    return x;
}

template<typename Integer>
static inline constexpr Integer compact_bits_2(Integer x) {
    x &= 0x5555555555555555;
    x = (x ^ (x >>  1))  & 0x3333333333333333;
    x = (x ^ (x >>  2))  & 0x0f0f0f0f0f0f0f0f;
    x = (x ^ (x >>  4))  & 0x00ff00ff00ff00ff;
    x = (x ^ (x >>  8))  & 0x0000ffff0000ffff;
    x = (x ^ (x >>  16)) & 0x00000000ffffffff;
    return x;
}
#endif

} //::morton

} //::aether
