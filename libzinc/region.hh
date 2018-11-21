#pragma once

#include <cstdint>
#include <cassert>

#include <algorithm>
#include <vector>
#include <tuple>
#include <variant>
#include <type_traits>

#include <libzinc/encoding.hh>
#include <libzinc/interval.hh>
#include <immintrin.h>

namespace zinc {

namespace morton {

//https://en.wikipedia.org/wiki/Linear_octree
//https://geidav.wordpress.com/2014/08/18/advanced-octrees-2-node-representations/
//(see Linear (hashed) Octrees section
template <uint32_t Dimension, uint32_t BitsPerDimension, typename T = std::monostate>
struct region {
    using interval_type = morton::detail::interval<Dimension, BitsPerDimension, T>;
    std::vector<interval_type> intervals;

    template<typename M = std::monostate>
    friend bool operator==(const region& lhs, const region<Dimension, BitsPerDimension, M>& rhs) {
        return lhs.intervals == rhs.intervals;
    }

    friend region operator|(const region& lhs, const region& rhs) {
        region r = lhs;
        r |= rhs;
        return r;
    }

    template<typename M = std::monostate>
    friend region operator&(const region& lhs, const region<Dimension, BitsPerDimension, M>& rhs) {
        region r = lhs;
        r &= rhs;
        return r;
    }

    friend void operator|=(region& lhs, const region& rhs) {
        {
            assert(std::is_sorted(lhs.intervals.begin(), lhs.intervals.end()) && std::is_sorted(rhs.intervals.begin(), rhs.intervals.end()));
            size_t mid = lhs.intervals.size();
            lhs.intervals.insert(lhs.intervals.end(), rhs.intervals.begin(), rhs.intervals.end());
            std::inplace_merge(lhs.intervals.begin(), lhs.intervals.begin() + mid, lhs.intervals.end());
        }
        std::vector<bool> to_delete(lhs.intervals.size());
        for (size_t i = 0; i < lhs.intervals.size() - 1; i++) {
            auto end = lhs.intervals[i].end;
            size_t j;
            for (j = i + 1;
                j < lhs.intervals.size() &&
                lhs.intervals[j].start-1 <= end && lhs.intervals[j].data_equals(lhs.intervals[i])
            ; j++) {
                end = std::max(end, lhs.intervals[j].end);
            }
            lhs.intervals[i].end = end;
            for (size_t k = i+1; k < j; k++) {
                to_delete[k] = true;
            }
            i = j-1; // to skip over the merged intervals
        }
        auto it = to_delete.begin();
        lhs.intervals.erase(
            std::remove_if(lhs.intervals.begin(), lhs.intervals.end(),
                [&it](const morton::detail::interval<Dimension, BitsPerDimension, T>& d) -> bool { return *it++; }
            )
        , lhs.intervals.end());
    }

    // this assumes regions contain a sorted list of morton intervals
    template<typename M = std::monostate>
    friend void operator&=(region& lhs, const region<Dimension, BitsPerDimension, M>& rhs) {
        assert(std::is_sorted(lhs.intervals.begin(), lhs.intervals.end()) && std::is_sorted(rhs.intervals.begin(), rhs.intervals.end()));
        auto lhs_it = lhs.intervals.begin();
        auto rhs_it = rhs.intervals.begin();
        std::vector<morton::detail::interval<Dimension, BitsPerDimension, T>> out = {};
        while(lhs_it != lhs.intervals.end() && rhs_it != rhs.intervals.end()){
            if (lhs_it->end < rhs_it->start) {
                ++lhs_it;
                continue;
            } else if (rhs_it->end < lhs_it->start) {
                ++rhs_it;
                continue;
            }
            auto s = std::max(lhs_it->start, rhs_it->start);
            auto e = std::min(lhs_it->end, rhs_it->end);
            out.push_back(morton::detail::interval<Dimension, BitsPerDimension, T>{s,e, lhs_it->data});
            if (lhs_it->end < rhs_it->end){
                ++lhs_it;
            } else if(rhs_it->end < lhs_it->end) {
                ++rhs_it;
            } else {
                ++rhs_it;
                ++lhs_it;
            }
        }
        lhs.intervals = out;
    }

    template<typename M = std::monostate>
    friend void operator-=(region& lhs, const region<Dimension, BitsPerDimension, M>& rhs) {
        assert(std::is_sorted(lhs.intervals.begin(), lhs.intervals.end()) && std::is_sorted(rhs.intervals.begin(), rhs.intervals.end()));
        auto lhs_it = lhs.intervals.begin();
        auto rhs_it = rhs.intervals.begin();
        std::vector<morton::detail::interval<Dimension, BitsPerDimension, T>> out {};
        morton_code<2, 32> s{0};
        if (lhs_it != lhs.intervals.end()) {
             s = lhs_it->start;
        }
        while(lhs_it != lhs.intervals.end() && rhs_it != rhs.intervals.end()){
            if (lhs_it->end < rhs_it->start ){ //if the lhs is entirely behind the rhs push the lhs
                // lhs: |------|
                // rhs:          |--|
                out.push_back(morton::detail::interval<Dimension, BitsPerDimension, T>{s,lhs_it->end, lhs_it->data});
                ++lhs_it;
                if (lhs_it != lhs.intervals.end()) {
                    s = lhs_it->start;
                }
                continue;
            }
            if (s > rhs_it->end){ // if the rhs is entirely behind the lhs, push the rhs
                // lhs:      |---|
                // rhs:|--|
                ++rhs_it;
                continue;
            }
            // we know now that the rhs collides with the lhs
            if (s >= rhs_it->start){
                //      A       |  |     B
                // lhs:  |---|  |or|     |---|
                // rhs:|--|     |  | |-------|
                if (lhs_it->end <= rhs_it->end){ // case B - move onto the next interval
                    ++lhs_it;
                    if (lhs_it != lhs.intervals.end()) {
                        s = lhs_it->start;
                    }
                } else { // case A - split the interval and check the next rhs segment
                    s = rhs_it->end + 1;
                    ++rhs_it;
                }
                continue;
            }
            // s must be < rhs
            //      A     |  |     B
            // lhs:|---|  |or| |---|
            // rhs: |-|   |  |   |----|
            out.push_back({s,rhs_it->start-1,lhs_it->data});
            if (rhs_it->end < lhs_it->end) { // case A
                s = rhs_it->end +1;
                ++rhs_it;
            } else { // case B
                ++lhs_it;
                if (lhs_it != lhs.intervals.end()) {
                    s = lhs_it->start;
                }
            }
        }
        if (lhs_it != lhs.intervals.end()) {
            //push the last one, then copy the rest
            out.push_back({s,lhs_it->end, lhs_it->data});
            out.insert(out.end(), ++lhs_it, lhs.intervals.end());
        }
        lhs.intervals = out;
    }

    template<typename M = std::monostate>
    friend region operator-(const region& lhs, const region<Dimension, BitsPerDimension, M>& rhs) {
        region x = lhs;
        x -= rhs;
        return x;
    }

    template<typename M>
    bool intersects(const region<Dimension, BitsPerDimension, M>& rhs) const;
    bool empty() const;
    uint64_t area() const;
    bool contains(const morton_code<Dimension, BitsPerDimension> c) const {
        for (auto& i: intervals) {
            if (c.code < i.start) {
                return false;
            } else if (c.code <= i.end) {
                return true;
            }
        }
        return false;
    };
    std::vector<detail::interval<Dimension, BitsPerDimension>> to_cells() const;
    std::vector<detail::interval<Dimension, BitsPerDimension>> to_cells(size_t max_level) const;
    std::vector<std::pair<uint64_t,uint64_t>> count_cells() const;
};

    template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
    template<typename M>
    bool region<Dimension, BitsPerDimension, T>::intersects(const region<Dimension, BitsPerDimension, M>& rhs) const {
        assert(std::is_sorted(intervals.begin(), intervals.end()) && std::is_sorted(rhs.intervals.begin(), rhs.intervals.end()));
        auto lhs_it = intervals.begin();
        auto rhs_it = rhs.intervals.begin();
        while(lhs_it != intervals.end() && rhs_it != rhs.intervals.end()){
            if (lhs_it->end < rhs_it->start) {
                ++lhs_it;
                continue;
            } else if (rhs_it->end < lhs_it->start) {
                ++rhs_it;
                continue;
            }
            return true;
        }
        return false;
    }

    template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
    bool region<Dimension, BitsPerDimension, T>::empty() const {
        return intervals.empty();
    }

    template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
    uint64_t region<Dimension, BitsPerDimension, T>::area() const {
        uint64_t a = 0;
        for(auto &cell : intervals){
            a += cell.area();
        }
        return a;
    }

    template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
    std::vector<detail::interval<Dimension, BitsPerDimension>> region<Dimension, BitsPerDimension, T>::to_cells() const {
        std::vector<detail::interval<Dimension, BitsPerDimension>> v = {};
        for (auto &i : intervals){
            auto c = i.to_cells();
            v.insert(v.end(), c.begin(), c.end());
        }
        return v;
    }

    template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
    std::vector<detail::interval<Dimension, BitsPerDimension>> region<Dimension, BitsPerDimension, T>::to_cells(size_t max_level) const {
        std::vector<detail::interval<Dimension, BitsPerDimension>> v = {};
        for (auto &i : intervals){
            auto c = i.to_cells(max_level);
            v.insert(v.end(), c.begin(), c.end());
        }
        return v;
    }

    template<uint32_t Dimension, uint32_t BitsPerDimension, typename T>
    std::vector<std::pair<uint64_t,uint64_t>> region<Dimension, BitsPerDimension, T>::count_cells() const {
        std::vector<std::pair<uint64_t,uint64_t>> counts = {};
        if (intervals.empty()) { return counts;}
        counts = intervals[0].count_cells();
        for(size_t i = 1; i < intervals.size(); i++){
            auto v = intervals[i].count_cells();
            for(auto it = v.begin(); it != v.end(); ++it){
                size_t j = 0;
                for(; j < counts.size(); j++){
                    if(counts[j].first < it->first && j < counts.size() -1 && counts[j+1].first > it->first){
                        //insert if between two existing values
                        counts.emplace(counts.begin()+j,*it);
                        break;
                    }
                    if(counts[j].first == it->first){
                        counts[j].second += it->second;
                        break;
                    }
                }
                if (j == counts.size()){
                    counts.push_back(*it);
                }
            }
        }
        return counts;
    }

template<typename T>
static region<2,32,T> cell_to_region(uint64_t code, uint64_t level, T data) {
    if constexpr (std::is_same<T, std::monostate>::value) {
        return {{{code, code + (1 << (level * 2)) - 1}}};
    } else {
        return {{{code, code + (1 << (level * 2)) - 1, data}}};
    }
}

} //::morton

}//zinc
