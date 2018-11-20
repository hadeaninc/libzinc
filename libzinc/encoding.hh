#pragma once

#include <cstdint>
#include <cassert>
#include <immintrin.h>
#include <array>
#include <tuple>

#include <libzinc/util.hh>

template<uint32_t Dimension, uint32_t BitsPerDimension>
struct morton_code {
    static_assert((Dimension == 2 && BitsPerDimension == 32) || (Dimension == 3 && BitsPerDimension == 21),
        "only 2D and 32 bits, or 3D and 21 bits are currently supported.");
};

using aether::morton::__morton_2_x_mask;
using aether::morton::__morton_2_y_mask;
using aether::morton::__morton_3_x_mask;
using aether::morton::__morton_3_y_mask;
using aether::morton::__morton_3_z_mask;
using aether::morton::expand_bits_2;
using aether::morton::compact_bits_2;
using aether::morton::expand_bits_3;
using aether::morton::compact_bits_3;

template<>
struct morton_code<2, 32> {
    uint64_t data;
    static constexpr uint32_t dimension = 2;
    static constexpr uint32_t max_level = 32;
    operator uint64_t() const {
        return data;
    }
    morton_code<2, 32>(uint64_t _data): data(_data) {};
    static morton_code<2, 32> encode(std::array<uint32_t, 2> p) {
        return {
            (expand_bits_2<uint64_t>(std::get<0>(p)) << 0) |
            (expand_bits_2<uint64_t>(std::get<1>(p)) << 1)
        };
    }
    static std::array<uint32_t, 2> decode(const struct morton_code<2, 32> code) {
        return {
            static_cast<uint32_t>(compact_bits_2<uint64_t>(code.data >> 0)),
            static_cast<uint32_t>(compact_bits_2<uint64_t>(code.data >> 1)),
        };
    }
    friend void operator-=(morton_code<2, 32>& lhs, const morton_code<2, 32>& rhs) {
        uint64_t x = (lhs.data & __morton_2_x_mask) - (rhs.data & __morton_2_x_mask);
        uint64_t y = (lhs.data & __morton_2_y_mask) - (rhs.data & __morton_2_y_mask);
        lhs.data = (x & __morton_2_x_mask) | (y & __morton_2_y_mask);
    }

    friend void operator+=(morton_code<2, 32>& lhs, const morton_code<2, 32>& rhs) {
        uint64_t x = (lhs.data | ~__morton_2_x_mask) + (rhs.data & __morton_2_x_mask);
        uint64_t y = (lhs.data | ~__morton_2_y_mask) + (rhs.data & __morton_2_y_mask);
        lhs.data = (x & __morton_2_x_mask) | (y & __morton_2_y_mask);
    }
};

template<>
struct morton_code<3, 21> {
    uint64_t data;
    static constexpr uint32_t dimension = 3;
    static constexpr uint32_t max_level = 21;

    static morton_code<3, 21> encode(std::array<uint32_t, 3> p) {
        return {
            (expand_bits_2<uint64_t>(std::get<0>(p)) << 0) |
            (expand_bits_2<uint64_t>(std::get<1>(p)) << 1) |
            (expand_bits_2<uint64_t>(std::get<2>(p)) << 2)
        };
    }
    static std::array<uint32_t, 3> decode(const struct morton_code<3, 21> code) {
        return {
            static_cast<uint32_t>(compact_bits_2<uint64_t>(code.data >> 0)),
            static_cast<uint32_t>(compact_bits_2<uint64_t>(code.data >> 1)),
            static_cast<uint32_t>(compact_bits_2<uint64_t>(code.data >> 2)),
        };
    }
    friend void operator+=(morton_code<3, 21>& lhs, const morton_code<3, 21>& rhs) {
        uint64_t x = (lhs.data | ~__morton_3_x_mask) + (rhs.data & __morton_3_x_mask);
        uint64_t y = (lhs.data | ~__morton_3_y_mask) + (rhs.data & __morton_3_y_mask);
        uint64_t z = (lhs.data | ~__morton_3_z_mask) + (rhs.data & __morton_3_z_mask);
        lhs.data = (x & __morton_3_x_mask) | (y & __morton_3_y_mask) | (z & __morton_3_z_mask);
    }

    friend void operator-=(morton_code<3, 21>& lhs, const morton_code<3, 21>& rhs) {
        uint64_t x = (lhs.data & __morton_3_x_mask) - (rhs.data & __morton_3_x_mask);
        uint64_t y = (lhs.data & __morton_3_y_mask) - (rhs.data & __morton_3_y_mask);
        uint64_t z = (lhs.data & __morton_3_z_mask) - (rhs.data & __morton_3_z_mask);
        lhs.data = (x & __morton_3_x_mask) | (y & __morton_3_y_mask) | (z & __morton_3_z_mask);
    }
};
