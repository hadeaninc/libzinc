#pragma once

#include <cstdint>
#include <cassert>

#include <vector>
#include <tuple>
#include <optional>

#include <libzinc/encoding.hh>
#include <libzinc/util.hh>
#include <immintrin.h>

namespace zinc {

namespace morton {

namespace detail {

template<uint32_t Dimension, uint32_t BitsPerDimension, typename T = std::monostate>
struct interval {
    morton_code<2, 32> start, end;
    T data {};
    interval(morton_code<2, 32> _start, morton_code<2, 32> _end): start(_start), end(_end) {};
    interval(morton_code<2, 32> _start, morton_code<2, 32> _end, T _t): start(_start), end(_end), data(_t) {};

    template<typename M>
    friend bool operator==(const interval& lhs, const interval<Dimension, BitsPerDimension, M>& rhs) {
        if constexpr (std::is_same<T, std::monostate>::value || std::is_same<M, std::monostate>::value) {
            return std::tie(lhs.start, lhs.end) == std::tie(rhs.start, rhs.end);
        } else if constexpr (std::is_same<T, M>::value) {
            return std::tie(lhs.start, lhs.end, lhs.data) == std::tie(rhs.start, rhs.end, rhs.data);
        } else {
            static_assert(std::is_same<T, std::monostate>::value || std::is_same<M, std::monostate>::value || std::is_same<T, M>::value, "the data types on intervals should be either the same or std::monostate");
        }
    }

    friend bool operator!=(const interval& lhs, const interval& rhs) {
        return !(lhs == rhs);
    }

    template<typename M>
    friend bool operator<(const interval& lhs, const interval<Dimension, BitsPerDimension, M>& rhs) {
        if constexpr (std::is_same<T, std::monostate>::value || std::is_same<M, std::monostate>::value) {
            return std::tie(lhs.start, lhs.end) < std::tie(rhs.start, rhs.end);
        } else if constexpr (std::is_same<T, M>::value) {
            return std::tie(lhs.start, lhs.end, lhs.data) < std::tie(rhs.start, rhs.end, rhs.data);
        } else {
            static_assert(std::is_same<T, std::monostate>::value || std::is_same<M, std::monostate>::value || std::is_same<T, M>::value, "the data types on intervals should be either the same or std::monostate");
        }
    }

    friend bool operator>(const interval& lhs, const interval& rhs) {
        return rhs < lhs;
    }

    friend bool operator<=(const interval& lhs, const interval& rhs) {
        return !(lhs > rhs);
    }

    friend bool operator>=(const interval& lhs, const interval& rhs) {
        return !(lhs < rhs);
    }

    bool data_equals(const interval& rhs) const;

    std::optional<interval> intersect(const interval& rhs) const;

    bool contains(const morton_code<Dimension, BitsPerDimension> c) const {
        return c >= start && c <= end;
    }

    uint64_t area() const;

    uint64_t start_alignment() const;

    uint64_t end_alignment() const;

    std::vector<interval> to_cells() const;

    std::vector<interval> to_cells(size_t max_level) const;
    // this returns an sorted map of the cells and size.
    // e.g. 3 cells of size 1, 2 cells of size 2, and one cell of size 3.
    // where size 1 = 1 on each side, size 2 = 2, size 3  = 4
    std::vector<std::pair<uint64_t,uint64_t>> count_cells() const;
};

template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
bool interval<Dimension, BitsPerDimension, T>::data_equals(const detail::interval<Dimension,BitsPerDimension,T>& rhs) const {
    if constexpr (std::is_same<T, std::monostate>::value) {
        return true;
    } else {
        return data == rhs.data;
    }
}

template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
std::optional< interval<Dimension,BitsPerDimension,T> > interval<Dimension, BitsPerDimension, T>::intersect(const interval<Dimension,BitsPerDimension,T>& rhs) const {
    uint64_t i_start = std::max(start, rhs.start);
    uint64_t i_end = std::min(end, rhs.end);
    if (i_start > i_end) {
      return std::nullopt;
    }
    return std::optional{interval<Dimension, BitsPerDimension>{i_start,i_end,data}};
}

template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
uint64_t interval<Dimension, BitsPerDimension, T>::area() const {
    assert(start <= end);
    return end + 1 - start;
}

template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
uint64_t interval<Dimension, BitsPerDimension, T>::start_alignment() const {
    return start != 0 ? __builtin_ctzll(start) / Dimension : std::numeric_limits<uint64_t>::max();
}

template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
uint64_t interval<Dimension, BitsPerDimension, T>::end_alignment() const {
    return end != 0 ? __builtin_ctzll(end) / Dimension : std::numeric_limits<uint64_t>::max();
}

template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
std::vector<interval<Dimension,BitsPerDimension,T>> interval<Dimension, BitsPerDimension, T>::to_cells() const {
    uint64_t s = start;
    std::vector<interval<Dimension,BitsPerDimension,T>> v = {};
    while(s <= end){
        auto amax = get_align_max<Dimension, BitsPerDimension>(s,end);
        v.push_back({s,amax,data});
        s = amax+1;
    }
    return v;
}

template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
std::vector<interval<Dimension,BitsPerDimension,T>> interval<Dimension, BitsPerDimension, T>::to_cells(size_t max_level) const {
    uint64_t s = start;
    std::vector<interval<Dimension,BitsPerDimension,T>> v = {};
    while(s <= end){
        auto amax = std::min(s + (1 << max_level), get_align_max<Dimension, BitsPerDimension>(s,end));
        v.push_back({s,amax,data});
        s = amax+1;
    }
    return v;
}
// this returns an sorted map of the cells and size.
// e.g. 3 cells of size 1, 2 cells of size 2, and one cell of size 3.
// where size 1 = 1 on each side, size 2 = 2, size 3  = 4
template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
std::vector<std::pair<uint64_t,uint64_t>> interval<Dimension, BitsPerDimension, T>::count_cells() const {
    uint64_t s = start;
    std::vector<std::pair<uint64_t,uint64_t>> v = {};
    while(s <= end){
        auto amax = get_align_max<Dimension, BitsPerDimension>(s,end);
        auto diff = fast_log2(1 + amax - s)/Dimension;
        size_t i = 0;
        bool found = false;
        for(; i < v.size(); i++){
            if (v[i].first == diff) {
                found = true;
                break;
            }
        }
        if (found) {
            v[i].second += 1;
        } else {
            v.push_back({diff,1});
        }
        s = amax+1;
    }
    std::sort(v.begin(), v.end());
    return v;
}

template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
static std::pair<uint64_t,uint64_t> get_parent_cell(interval<Dimension, BitsPerDimension, T> interval){
    auto level = get_unifying_level<Dimension>(interval.start, interval.end);
    auto parent = get_parent_morton_aligned<Dimension>(interval.start, level);
    return {parent, level};
}

template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
static interval<Dimension, BitsPerDimension, T> get_parent_cell(const interval<Dimension, BitsPerDimension, T> interval){
    std::pair<uint64_t,uint64_t> pair = get_parent_cell<Dimension>(interval);
    if (std::is_same<T, std::monostate>::value){
        return {pair.first, pair.first + get_morton_code<Dimension>(pair.second)};
    } else {
        return {pair.first, pair.first + get_morton_code<Dimension>(pair.second), interval.data};
    }
}


} //::detail

} //::morton

} //::zinc
