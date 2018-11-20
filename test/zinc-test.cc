#include <cstdio>
#include <cstdint>
#include <cassert>
#include <iostream>
#include <cmath>

#include <libzinc/zinc.hh>

int main() {
    {
        assert(!(aether::morton::AABB<2, 32>{3, 12}).is_morton_aligned());
        assert(!(aether::morton::AABB<2, 32>{15, 48}).is_morton_aligned());
        assert(!(aether::morton::AABB<2, 32>{1, 2}).is_morton_aligned());
        assert(!(aether::morton::AABB<2, 32>{16, 23}).is_morton_aligned());
        assert((aether::morton::AABB<2, 32>{0, 3}).is_morton_aligned());
        assert((aether::morton::AABB<2, 32>{8, 11}).is_morton_aligned());
        assert((aether::morton::AABB<2, 32>{12, 15}).is_morton_aligned());
        assert((aether::morton::AABB<2, 32>{4, 7}).is_morton_aligned());
        assert((aether::morton::AABB<2, 32>{0, 0}).is_morton_aligned());
        assert((aether::morton::AABB<2, 32>{2, 2}).is_morton_aligned());
        assert((aether::morton::AABB<2, 32>{7, 7}).is_morton_aligned());
        assert((aether::morton::AABB<2, 32>{0, 15}).is_morton_aligned());
    }
    {
        uint64_t x;
        x = 1; assert(aether::morton::fast_log2(x) == (uint64_t)log2(x));
        x = 2; assert(aether::morton::fast_log2(x) == (uint64_t)log2(x));
        x = 3; assert(aether::morton::fast_log2(x) == (uint64_t)log2(x));
        x = 4; assert(aether::morton::fast_log2(x) == (uint64_t)log2(x));
        x = 8; assert(aether::morton::fast_log2(x) == (uint64_t)log2(x));
    }
    {
        aether::morton::AABB<2, 32> aabb = {51, 193};
        uint64_t litmax = 0, bigmin = 0;
        std::tie(litmax, bigmin) = aabb.morton_get_next_address();
        assert(litmax == 107);
        assert(bigmin == 145);
        aabb.max = 107;
        std::tie(litmax, bigmin) = aabb.morton_get_next_address();
        assert(litmax == 63);
        assert(bigmin == 98);
        aabb.min = 98;
        std::tie(litmax, bigmin) = aabb.morton_get_next_address();
        assert(litmax == 99);
        assert(bigmin == 104);
        aabb.min = 145;
        aabb.max = 193;
        std::tie(litmax, bigmin) = aabb.morton_get_next_address();
        assert(litmax == 149);
        assert(bigmin == 192);
    }

    {
        bool good =
            (aether::morton::region<2, 32>{{
                {0, 1},
                {1, 2},
                {2, 3},
                {3, 4},
            }} | aether::morton::region<2, 32>{{
            }}) == aether::morton::region<2, 32>{{
                {0, 4},
            }};
        assert(good);
    }

    {
        aether::morton::AABB<2, 32> aabb {0, 12};

        aether::morton::region<2, 32> region = aabb.to_cells();

        std::vector<aether::morton::detail::interval<2, 32>>::iterator it = region.intervals.begin();
        assert(it->start == 0);
        assert(it->end == 3);
        ++it;
        assert(it->start == 4);
        assert(it->end == 4);
        ++it;
        assert(it->start == 6);
        assert(it->end == 6);
        ++it;
        assert(it->start == 8);
        assert(it->end == 8);
        ++it;
        assert(it->start == 9);
        assert(it->end == 9);
        ++it;
        assert(it->start == 12);
        assert(it->end == 12);
        ++it;
        assert(it == region.intervals.end());
    }

    {
        aether::morton::AABB<2, 32> aabb {0, 12};

        aether::morton::region<2, 32> region = aabb.to_intervals();

        std::vector<aether::morton::detail::interval<2, 32>>::iterator it = region.intervals.begin();
        assert(it->start == 0);
        assert(it->end == 4);
        ++it;
        assert(it->start == 6);
        assert(it->end == 6);
        ++it;
        assert(it->start == 8);
        assert(it->end == 9);
        ++it;
        assert(it->start == 12);
        assert(it->end == 12);
        ++it;
        assert(it == region.intervals.end());
    }

    {
        aether::morton::AABB<2, 32> aabb {0, 12};

        auto it = aabb.begin();
        assert(it->start == 0);
        assert(it->end == 4);
        ++it;
        assert(it->start == 6);
        assert(it->end == 6);
        ++it;
        assert(it->start == 8);
        assert(it->end == 9);
        ++it;
        assert(it->start == 12);
        assert(it->end == 12);
        ++it;
        assert(it == aabb.end());
    }

    {
        // test that 0,0 works with the uninitialised check in the iterator.progress() function
        aether::morton::AABB<2, 32> aabb {0, 0};

        auto it = aabb.begin();
        assert(it->start == 0);
        assert(it->end == 0);
        ++it;
        assert(it == aabb.end());
    }

    {
        auto l = aether::morton::detail::interval<2, 32>{0,5};
        auto r = aether::morton::detail::interval<2, 32>{2,7};
        auto i = l.intersect(r);
        assert(i.has_value());
        assert(i.value().start == 2);
        assert(i.value().end == 5);
        r = aether::morton::detail::interval<2, 32>{7,23};
        i = l.intersect(r);
        assert(!i.has_value());
    }

    {
        auto a = aether::morton::region<2, 32>{
            {aether::morton::detail::interval<2, 32>{0,1},aether::morton::detail::interval<2, 32>{3,4}}
        };
        auto b = aether::morton::region<2, 32>{
            {aether::morton::detail::interval<2, 32>{0,3}}
        };
        a &= b;
        assert(a.intervals.size() == 2);
        auto x = aether::morton::detail::interval<2, 32>{0,1};
        auto y = aether::morton::detail::interval<2, 32>{3,3};
        assert(a.intervals[0] == x);
        assert(a.intervals[1] == y);
    }

    {
        auto a = aether::morton::region<2, 32>{
            {aether::morton::detail::interval<2, 32>{0,63}}
        };
        auto b = aether::morton::region<2, 32>{
            {aether::morton::detail::interval<2, 32>{0,3},aether::morton::detail::interval<2, 32>{24,27}, aether::morton::detail::interval<2, 32>{48,63}}
        };
        a &= b;
        assert(a.intervals.size() == 3);
        auto x = aether::morton::detail::interval<2, 32>{0,3};
        auto y = aether::morton::detail::interval<2, 32>{24,27};
        auto z = aether::morton::detail::interval<2, 32>{48,63};
        assert(a.intervals[0] == x);
        assert(a.intervals[1] == y);
        assert(a.intervals[2] == z);
    }

    {
        auto a = aether::morton::region<2, 32>{{}};
        auto b = aether::morton::region<2, 32>{{}};
        a &= b;
        assert(a.intervals.size() == 0);
    }

    {
        auto a = aether::morton::region<2, 32>{{aether::morton::detail::interval<2, 32>{0,1}}};
        auto b = aether::morton::region<2, 32>{{}};
        a &= b;
        assert(a.intervals.size() == 0);
    }

    {
        auto a = aether::morton::region<2, 32>{{}};
        auto b = aether::morton::region<2, 32>{{aether::morton::detail::interval<2, 32>{0,1}}};
        a &= b;
        assert(a.intervals.size() == 0);
    }

    {
        auto a = aether::morton::region<2, 32>{
            {aether::morton::detail::interval<2, 32>{0,15}
            ,aether::morton::detail::interval<2, 32>{16,17}
            ,aether::morton::detail::interval<2, 32>{20,21}
            ,aether::morton::detail::interval<2, 32>{24,25}
            ,aether::morton::detail::interval<2, 32>{28,29}}
        };
        auto b = aether::morton::region<2, 32>{
            {aether::morton::detail::interval<2, 32>{0,1}
            ,aether::morton::detail::interval<2, 32>{4,5}
            ,aether::morton::detail::interval<2, 32>{8,9}
            ,aether::morton::detail::interval<2, 32>{12,13}
            ,aether::morton::detail::interval<2, 32>{16,31}}
        };
        a &= b;
        assert(a.intervals.size() == 8);
        auto r1 = aether::morton::detail::interval<2, 32>{0,1};
        auto r2 = aether::morton::detail::interval<2, 32>{4,5};
        auto r3 = aether::morton::detail::interval<2, 32>{8,9};
        auto r4 = aether::morton::detail::interval<2, 32>{12,13};
        auto r5 = aether::morton::detail::interval<2, 32>{16,17};
        auto r6 = aether::morton::detail::interval<2, 32>{20,21};
        auto r7 = aether::morton::detail::interval<2, 32>{24,25};
        auto r8 = aether::morton::detail::interval<2, 32>{28,29};
        assert(a.intervals[0] == r1);
        assert(a.intervals[1] == r2);
        assert(a.intervals[2] == r3);
        assert(a.intervals[3] == r4);
        assert(a.intervals[4] == r5);
        assert(a.intervals[5] == r6);
        assert(a.intervals[6] == r7);
        assert(a.intervals[7] == r8);
    }

    {
        auto a = aether::morton::region<2, 32>{
            {aether::morton::detail::interval<2, 32>{0,1},aether::morton::detail::interval<2, 32>{3,4}}
        };
        auto b = aether::morton::region<2, 32>{
            {aether::morton::detail::interval<2, 32>{0,3}}
        };
        a -= b;
        assert(a.intervals.size() == 1);
        auto r1 = aether::morton::detail::interval<2, 32>{4,4};
        assert(a.intervals[0] == r1);
    }

    {
        auto a = aether::morton::region<2, 32>{
            {aether::morton::detail::interval<2, 32>{0,1},aether::morton::detail::interval<2, 32>{3,4}}
        };
        auto b = aether::morton::region<2, 32>{{}};
        a -= b;
        assert(a.intervals.size() == 2);
        auto r1 = aether::morton::detail::interval<2, 32>{0,1};
        auto r2 = aether::morton::detail::interval<2, 32>{3,4};
        assert(a.intervals[0] == r1);
        assert(a.intervals[1] == r2);
    }

    {
        auto a = aether::morton::region<2, 32>{{}};
        auto b = aether::morton::region<2, 32>{{}};
        a -= b;
        assert(a.intervals.size() == 0);
    }

    {
        auto a = aether::morton::region<2, 32>{
            {aether::morton::detail::interval<2, 32>{0,2},aether::morton::detail::interval<2, 32>{4,6}}
        };
        auto b = aether::morton::region<2, 32>{{
            {0,0},{2,4},{6,6}
        }};
        a -= b;
        assert(a.intervals.size() == 2);
        auto r1 = aether::morton::detail::interval<2, 32>{1,1};
        auto r2 = aether::morton::detail::interval<2, 32>{5,5};
        assert(a.intervals[0] == r1);
        assert(a.intervals[1] == r2);
        b -= aether::morton::region<2, 32>{{aether::morton::detail::interval<2, 32>{0,2},aether::morton::detail::interval<2, 32>{4,6}}};
        assert(b.intervals.size() == 1);
        r1 = aether::morton::detail::interval<2, 32>{3,3};
        assert(b.intervals[0] == r1);
    }

    {
        auto a = aether::morton::region<2, 32>{{
            {0,5}
        }};
        auto b = aether::morton::region<2, 32>{{
            {1,1},{4,4}
        }};
        a -= b;
        assert(a.intervals.size() == 3);
        auto r1 = aether::morton::detail::interval<2, 32>{0,0};
        auto r2 = aether::morton::detail::interval<2, 32>{2,3};
        auto r3 = aether::morton::detail::interval<2, 32>{5,5};
        assert(a.intervals[0] == r1);
        assert(a.intervals[1] == r2);
        assert(a.intervals[2] == r3);
        b -= aether::morton::region<2, 32>{{
            {0,5}
        }};
        assert(b.intervals.size() == 0);
    }

    {
        auto a = aether::morton::region<2, 32>{{
            {0,2}
        }};
        auto b = aether::morton::region<2, 32>{{
            {4,6}
        }};
        a -= b;
        assert(a.intervals.size() == 1);
        auto r1 = aether::morton::detail::interval<2, 32>{0,2};
        assert(a.intervals[0] == r1);
        b -= a;
        assert(b.intervals.size() == 1);
        r1 = aether::morton::detail::interval<2, 32>{4,6};
        assert(b.intervals[0] == r1);
    }

    {
        aether::morton::region<2,32> a = {{
            {0,3},
            {8,11},
            {13,17},
            {21,24},
            {26,28},
            {31,34},
            {36,40},
            {42,51}
        }};
        aether::morton::region<2,32> b = {{
            {4,7},
            {9,12},
            {14,15},
            {19,22},
            {25,29},
            {33,38},
            {42,42},
            {45,46},
            {48,48},
            {50,52},
        }};
        aether::morton::region<2,32> r = {{
            {0, 3},
            {8, 8},
            {13, 13},
            {16, 17},
            {23, 24},
            {31, 32},
            {39, 40},
            {43, 44},
            {47, 47},
            {49, 49}
        }};
        a -= b;
        assert(a==r);
    }

    {
        aether::morton::region<2,32> a = {{
            {0,3},
            {8,11},
            {13,17},
            {21,24},
            {26,28},
            {31,34},
            {36,40},
            {42,51}
        }};
        aether::morton::region<2,32> b = {{
            {4,7},
            {9,12},
            {14,15},
            {19,22},
            {25,29},
            {33,38},
            {42,42},
            {45,46},
            {48,48},
            {50,52},
        }};
        aether::morton::region<2,32> r = {{
            {0, 3},
            {8, 8},
            {13, 13},
            {16, 17},
            {23, 24},
            {31, 32},
            {39, 40},
            {43, 44},
            {47, 47},
            {49, 49}
        }};
        a -= b;
        assert(a==r);
    }

    {
        auto interval = aether::morton::detail::interval<2, 32>{0,0};
        assert(interval.area() == 1);
        interval = aether::morton::detail::interval<2, 32>{0,1};
        assert(interval.area() == 2);
        interval = aether::morton::detail::interval<2, 32>{1,2};
        assert(interval.area() == 2);
    }

    {
        auto region = aether::morton::region<2, 32>{{
            {0,0},{2,5}
        }};
        assert(region.area() == 5);
        region = aether::morton::region<2, 32>{{}};
        assert(region.area() == 0);
    }

    {
        assert((aether::morton::detail::interval<2, 32>{1, 0}).start_alignment() == 0);
        assert((aether::morton::detail::interval<2, 32>{2, 0}).start_alignment() == 0);
        assert((aether::morton::detail::interval<2, 32>{4, 0}).start_alignment() == 1);
        assert((aether::morton::detail::interval<2, 32>{6, 0}).start_alignment() == 0);
        assert((aether::morton::detail::interval<2, 32>{8, 0}).start_alignment() == 1);
        assert((aether::morton::detail::interval<2, 32>{16, 0}).start_alignment() == 2);
        assert((aether::morton::detail::interval<2, 32>{0, 1}).end_alignment() == 0);
        assert((aether::morton::detail::interval<2, 32>{0, 2}).end_alignment() == 0);
        assert((aether::morton::detail::interval<2, 32>{0, 4}).end_alignment() == 1);
        assert((aether::morton::detail::interval<2, 32>{0, 6}).end_alignment() == 0);
        assert((aether::morton::detail::interval<2, 32>{0, 8}).end_alignment() == 1);
        assert((aether::morton::detail::interval<2, 32>{0, 16}).end_alignment() == 2);
    }

    {
        auto i = aether::morton::get_align_max<2,32>(0,21);
        assert(i == 15);
        i = aether::morton::get_align_max<2,32>(1,21);
        assert(i == 1);
        i = aether::morton::get_align_max<2,32>(3,21);
        assert(i == 3);
        i = aether::morton::get_align_max<2,32>(4,21);
        assert(i == 7);
        i = aether::morton::get_align_max<2,32>(0,64);
        assert(i == 63);
        i = aether::morton::get_align_max<2,32>(0,0);
        assert(i == 0);
        i = aether::morton::get_align_max<2,32>(16,21);
        assert(i == 19);
        i = aether::morton::get_align_max<2,32>(0,3);
        assert(i == 3);
        i = aether::morton::get_align_max<2,32>(12,31);
        assert(i == 15);
        i = aether::morton::get_align_max<2,32>(16,63);
        assert(i == 31);
    }

    {
        auto v = aether::morton::detail::interval<2, 32>{0,21}.count_cells();
        std::vector<std::pair<uint64_t,uint64_t>> i = {{0,2},{1,1},{2,1}};
        assert(v == i);

        v = aether::morton::detail::interval<2, 32>{0,3}.count_cells();
        i = {{1,1}};
        assert(v == i);

        v = aether::morton::detail::interval<2, 32>{0,63}.count_cells();
        i = {{3,1}};
        assert(v == i);

        v = aether::morton::detail::interval<2, 32>{1,63}.count_cells();
        i = {{0,3},{1,3},{2,3}};
        assert(v == i);
    }

    {
        auto region = aether::morton::region<2, 32>{{{0,21}, {23,31}}};
        std::vector<std::pair<uint64_t,uint64_t>> i = {{0,3},{1,3},{2,1}};
        assert(region.count_cells() == i);


        region = aether::morton::region<2, 32>{{}};
        i = {};
        assert(region.count_cells() == i);

        region = aether::morton::region<2, 32>{{{0,21}, {23,31}, {43,63}}};
        i = {{0,4},{1,4},{2,2}};
        assert(region.count_cells() == i);

        auto interval = aether::morton::detail::interval<2, 32>{0,21};
        region = aether::morton::region<2, 32>{{interval}};
        assert(region.count_cells() == interval.count_cells());
    }

    {
        // this tests for UB in bit shifts
        aether::morton::AABB<2,32> aabb = {4611686018427387648, 13835058055282164480ULL};
        auto region = aabb.to_intervals();
        assert(region.intervals.size() == 51);
    }

    {
        aether::morton::detail::interval<2,32> i = {0,15};
        auto r = i.to_cells();
        std::vector<aether::morton::detail::interval<2,32>> t = {i};
        assert(r == t);
        i = {1,15};
        t = {{1,1},{2,2},{3,3},{4,7},{8,11},{12,15}};
        r = i.to_cells();
        assert(r == t);
    }

    {
        aether::morton::region<2,32> r = {{
            {1,15},
            {57,57},
            {59,63}
        }};
        auto intervals = r.to_cells();
        std::vector<aether::morton::detail::interval<2,32>> t = {{1,1},{2,2},{3,3},{4,7},{8,11},{12,15},{57,57},{59,59},{60,63}};
        assert(intervals == t);
    }

    {
        aether::morton::region<2,32,uint64_t> a = {{
            {1,4,0},
            {17,31,1},
        }};
        aether::morton::region<2,32,uint64_t> b = {{
            {3,7,0},
            {13,16,1},
        }};
        aether::morton::region<2,32,uint64_t> r = {{
            {1,7,0},
            {13,31,1},
        }};
        a |= b;
        assert(a == r);
    }

    {
        uint64_t a = 254;
        auto a1 = aether::morton::get_parent_morton_aligned<2>(a, 1);
        auto a2 = aether::morton::get_parent_morton_aligned<2>(a, 2);
        auto a3 = aether::morton::get_parent_morton_aligned<2>(a, 3);
        auto a4 = aether::morton::get_parent_morton_aligned<2>(a, 4);
        assert(a1 == 252);
        assert(a2 == 240);
        assert(a3 == 192);
        assert(a4 == 0);
    }

    {
        assert(aether::morton::get_morton_code<2>(0) == 0);
        assert(aether::morton::get_morton_code<2>(1) == 3);
        assert(aether::morton::get_morton_code<2>(2) == 15);
        assert(aether::morton::get_morton_code<2>(3) == 63);
    }

    {
        assert((aether::morton::get_max_align_level<2,32>(0))  == 32);
        assert((aether::morton::get_max_align_level<2,32>(4))  == 1);
        assert((aether::morton::get_max_align_level<2,32>(16)) == 2);
        assert((aether::morton::get_max_align_level<2,32>(64)) == 3);
        assert((aether::morton::get_max_align_level<2,32>(3))  == 0);
    }

    {
        uint64_t a1 = aether::morton::get_unifying_level<2>(152, 156);
        uint64_t a2 = aether::morton::get_unifying_level<2>(152, 153);
        uint64_t a3 = aether::morton::get_unifying_level<2>(133, 152);
        uint64_t a4 = aether::morton::get_unifying_level<2>(0, 255);
        uint64_t a5 = aether::morton::get_unifying_level<2>(127, 128);
        uint64_t a6 = aether::morton::get_unifying_level<2>(127, 127);
        assert(a1 == 2);
        assert(a2 == 1);
        assert(a3 == 3);
        assert(a4 == 4);
        assert(a5 == 4);
        assert(a6 == 0);
    }
    
    return 0;
}
