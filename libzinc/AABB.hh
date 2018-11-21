#pragma once

#include <cstdint>
#include <cassert>

#include <algorithm>
#include <vector>
#include <tuple>
#include <variant>
#include <type_traits>

#include "encoding.hh"
#include "region.hh"
#include <immintrin.h>

namespace zinc {

namespace morton {

template<uint32_t Dimension, uint32_t BitsPerDimension>
struct AABB {
    // AABBs can have min = max represent a single interval.
    // That is min and max are inclusive values
    morton_code<Dimension, BitsPerDimension> min, max;

    AABB(morton_code<Dimension, BitsPerDimension> _min, morton_code<Dimension, BitsPerDimension> _max): min(_min), max(_max) {};

    class iterator_intervals {
        public:
            typedef std::input_iterator_tag iterator_category;
            typedef morton::detail::interval<Dimension, BitsPerDimension> value_type;
            typedef const value_type *pointer;
            typedef const value_type &reference;

        private:

            size_t iterator_index;
            value_type value;
            value_type curr;
            const AABB &parent_aabb;

            std::vector<AABB> inputs;

            bool is_end() const {
                return is_finished;
            }

            void progress() {
                while (!inputs.empty()) {
                    AABB aabb = inputs.back();
                    inputs.pop_back();
                    if (aabb.is_morton_aligned()) {
                        // This is checking that curr is uninitialised
                        // in the case that parent_aabb = {0,0} this still works as
                        // curr is always assigned to the value at the end.
                        if (curr == morton::detail::interval<Dimension, BitsPerDimension>{0,0}){
                            curr = aabb.to_cell();
                            continue;
                        } else if (curr.end + 1 == aabb.min) {
                            curr.end = aabb.max;
                            continue;
                        } else {
                            value = curr;
                            curr = aabb.to_cell();
                            iterator_index++;
                            return;
                        }
                    }
                    auto [litmax, bigmin] = aabb.morton_get_next_address();
                    AABB first = {aabb.min, litmax};
                    AABB second = {bigmin, aabb.max};
                    assert(first.max >= first.min);
                    assert(second.max >= second.min);
                    inputs.push_back(second);
                    inputs.push_back(first);
                }
                value = curr;
                is_finished = true;
                return;
            }


        public:
            bool is_finished;

            iterator_intervals(const AABB &_parent): value({0, 0}), curr({0, 0}), parent_aabb(_parent) {
                iterator_index = 0;
                inputs = {_parent};
                is_finished = false;
                progress();
            }

            void set_end() {
                is_finished = true;
            }

            iterator_intervals &operator++() {
                progress();
                return *this;
            }

            bool operator==(const iterator_intervals &i) const {
                return parent_aabb == i.parent_aabb && is_finished == i.is_finished && (is_finished == true || iterator_index == i.iterator_index);
            }

            bool operator!=(const iterator_intervals &i) const {
                return !(*this == i);
            }

            reference operator*() const {
                return value;
            }

            pointer operator->() const {
                return &value;
            }
    };

    typedef iterator_intervals iterator;

    iterator begin() const {
        return iterator(*this);
    }

    iterator end() const {
        auto i = iterator(*this);
        i.set_end();
        return i;
    }

    friend bool operator==(const AABB& lhs, const AABB& rhs){
        return lhs.min == rhs.min && lhs.max == rhs.max;
    }

    bool is_morton_aligned() const;

    morton::detail::interval<Dimension, BitsPerDimension> to_cell() const;

    // This generates a list of all morton aligned intervals (morton cells)
    // that are within the AABB
    // it generates them in a sorted order, from lowest interval to highest.
    region<Dimension, BitsPerDimension> to_cells() const;

    // This generates a list of all contiguous morton intervals (these are not necessarily aligned)
    // that are within the AABB
    // it generates them in a sorted order, from lowest interval to highest.
    region<Dimension, BitsPerDimension> to_intervals() const;

    uint64_t get_next_morton_outside(uint64_t m) const;

    uint64_t get_next_morton_inside(uint64_t m) const;

    //morton_get_next_address is complex, read these for more detail
    //https://en.wikipedia.org/wiki/Z-order_curve#Use_with_one-dimensional_data_structures_for_range_searching
    //https://stackoverflow.com/questions/30170783/how-to-use-morton-orderz-order-curve-in-range-search/34956693#34956693
    //https://raima.com/wp-content/uploads/COTS_embedded_database_solving_dynamic_pois_2012.pdf
    //http://cppedinburgh.uk/slides/201603-zcurves.pdf
    std::pair<morton_code<Dimension, BitsPerDimension>, morton_code<Dimension, BitsPerDimension>> morton_get_next_address();
};

template<uint32_t Dimension, uint32_t BitsPerDimension>
bool AABB<Dimension, BitsPerDimension>::is_morton_aligned() const {
    assert(max >= min);
    uint64_t align_max = min != 0 ? __builtin_ctzll(min) : std::numeric_limits<uint64_t>::max();
    uint64_t diff = max - min + 1;
    uint64_t align = __builtin_ctzll(diff);
    return
        align / Dimension <= align_max / Dimension &&
        __builtin_popcountll(diff) == 1 &&
        align % Dimension == 0;
}

template<uint32_t Dimension, uint32_t BitsPerDimension>
morton::detail::interval<Dimension, BitsPerDimension> AABB<Dimension, BitsPerDimension>::to_cell() const {
    assert(this->is_morton_aligned());
    return morton::detail::interval<Dimension, BitsPerDimension>{min, max};
}

// This generates a list of all morton aligned intervals (morton cells)
// that are within the AABB
// it generates them in a sorted order, from lowest interval to highest.
template<uint32_t Dimension, uint32_t BitsPerDimension>
region<Dimension, BitsPerDimension> AABB<Dimension, BitsPerDimension>::to_cells() const {
    assert(max >= min);
    std::vector<AABB> inputs {*this};
    std::vector<morton::detail::interval<Dimension, BitsPerDimension>> outputs;
    while (!inputs.empty()) {
        AABB aabb = inputs.back();
        inputs.pop_back();
        if (aabb.is_morton_aligned()) {
            outputs.push_back(aabb.to_cell());
            continue;
        }
        auto [litmax, bigmin] = aabb.morton_get_next_address();
        AABB first = {aabb.min, litmax};
        AABB second = {bigmin, aabb.max};
        assert(first.max >= first.min);
        assert(second.max >= second.min);
        inputs.push_back(second);
        inputs.push_back(first);
    }
    return {outputs};
}

// This generates a list of all contiguous morton intervals (these are not necessarily aligned)
// that are within the AABB
// it generates them in a sorted order, from lowest interval to highest.
template<uint32_t Dimension, uint32_t BitsPerDimension>
region<Dimension, BitsPerDimension> AABB<Dimension, BitsPerDimension>::to_intervals() const {
    assert(max >= min);
    std::vector<AABB> inputs {*this};
    std::vector<morton::detail::interval<Dimension, BitsPerDimension>> outputs;
    while (!inputs.empty()) {
        AABB aabb = inputs.back();
        inputs.pop_back();
        if (aabb.is_morton_aligned()) {
            //if the cell generated connects to the previous cell, merge.
            if (!outputs.empty() && outputs.back().end + 1 == aabb.min) {
                outputs.back().end = aabb.max;
            } else {
                outputs.push_back(aabb.to_cell());
            }
            continue;
        }
        auto [litmax, bigmin] = aabb.morton_get_next_address();
        AABB first = {aabb.min, litmax};
        AABB second = {bigmin, aabb.max};
        assert(first.max >= first.min);
        assert(second.max >= second.min);
        inputs.push_back(second);
        inputs.push_back(first);
    }
    return {outputs};
}

template<uint32_t Dimension, uint32_t BitsPerDimension>
uint64_t AABB<Dimension, BitsPerDimension>::get_next_morton_outside(uint64_t m) const {
    uint64_t min_x, min_y, max_x, max_y, x, y;
    auto pack = morton_code<Dimension, BitsPerDimension>::decode({min});
    min_x = pack[0];
    min_y = pack[1];
    pack = morton_code<Dimension, BitsPerDimension>::decode({max});
    max_x = pack[0];
    max_y = pack[1];
    pack = morton_code<Dimension, BitsPerDimension>::decode({m});
    x = pack[0];
    y = pack[1];
    assert(x == min_x || y == min_y);
    max_x += 1;
    max_y += 1;
    min_x = std::min(min_x, (uint64_t)(1LLU << 63));
    min_y = std::min(min_y, (uint64_t)(1LLU << 63));
    uint64_t l = std::min({__builtin_ctzll(min_x), __builtin_ctzll(min_y), __builtin_ctzll(max_x), __builtin_ctzll(max_x)}) * Dimension;
    m = (m >> l) << l;
    m = m + (1LLU << l);
    return m;
}

template<uint32_t Dimension, uint32_t BitsPerDimension>
uint64_t AABB<Dimension, BitsPerDimension>::get_next_morton_inside(uint64_t m) const {
    return 0;
}

//morton_get_next_address is complex, read these for more detail
//https://en.wikipedia.org/wiki/Z-order_curve#Use_with_one-dimensional_data_structures_for_range_searching
//https://stackoverflow.com/questions/30170783/how-to-use-morton-orderz-order-curve-in-range-search/34956693#34956693
//https://raima.com/wp-content/uploads/COTS_embedded_database_solving_dynamic_pois_2012.pdf
//http://cppedinburgh.uk/slides/201603-zcurves.pdf
template<uint32_t Dimension, uint32_t BitsPerDimension>
std::pair<morton_code<Dimension, BitsPerDimension>, morton_code<Dimension, BitsPerDimension>> AABB<Dimension, BitsPerDimension>::morton_get_next_address() {
    static_assert(Dimension == 2, "this needs more care to convert to 3D");
    uint64_t litmax = max;
    uint64_t bigmin = min;
    uint64_t index = 65 - __builtin_clzll(min ^ max);
    uint64_t mask = ~((1LLU << (index / 2)) - 1);
    uint64_t inc = 1LLU << ((index / 2) - 1);
    index %= 2;

    uint32_t part = (compact_bits_2(min >> index) & mask) + inc;
    bigmin &= ~(__morton_2_x_mask << index);
    bigmin |= expand_bits_2((uint64_t)part) << index;

    part -= 1;
    litmax &= ~(__morton_2_x_mask << index);
    litmax |= expand_bits_2((uint64_t)part) << index;

    return std::pair<morton_code<Dimension, BitsPerDimension>, morton_code<Dimension, BitsPerDimension>>({litmax}, {bigmin});
}

} //::morton

} //::zinc
