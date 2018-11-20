#pragma once

#include <cstdint>

#include <libzinc/AABB.hh>

//The coordinate of an octree cell.
struct tree_cell {
    uint64_t code;
    uint64_t level;

    void fix_code() {
        uint64_t dimension = 2;
        code = (code >> (dimension * level)) << (dimension * level);
    };

    bool check_overlap(uint64_t dim, tree_cell y) const {
        uint64_t l = (level > y.level) ? level : y.level;
        return (code >> (dim * l)) == (y.code >> (dim * l));
    }

    template<uint32_t Dimension, uint32_t BitsPerDimension>
    bool contains(const morton_code<Dimension, BitsPerDimension> c) const {
        return (c.data >> (level * c.dimension)) == (c.data >> (level * c.dimension));
    }

    aether::morton::region<2,32> region() const {
    	return {{{code, code + (1 << (level * 2)) - 1}}};
    }

    template<typename T>
    aether::morton::region<2,32, T> region(T data) const {
        return {{{code, code + (1 << (level * 2)) - 1}, data}};
    }
};
