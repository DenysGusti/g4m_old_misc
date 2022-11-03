#ifndef G4M_NEW_MISC_H
#define G4M_NEW_MISC_H

#include <map>
#include <vector>

#include <cmath>
#include <algorithm>
#include <ranges>

#include <iostream>
#include <string>
#include <format>

using namespace std;
namespace rv = ranges::views;

using ld = long double;

template<class T>
concept FloatingPointVector = same_as<vector<float>, T> || same_as<vector<double>, T> || same_as<vector<ld>, T>;

template<class T>
concept KeyArgument = floating_point<T> || FloatingPointVector<T>;

namespace g4m {
    template<KeyArgument IDX, floating_point VAL>
    class vipol {
    public:
        virtual ~vipol() = default;

        virtual VAL ip(const IDX &i) const noexcept = 0;

        VAL operator()(const IDX &i) const noexcept {
            return ip(i);
        };
    };

    template<floating_point Key, floating_point Value>
    class ipol : public vipol<Key, Value> {
    private:
        // find min or max key
        Key minOrMaxKey(const bool min) const noexcept {
            if (data.empty())
                return {};

            return (min ? data.begin()->first : data.rbegin()->first);
        }

        // find min or max value
        Value minOrMaxValue(const bool min) const noexcept {
            if (data.empty())
                return {};

            const auto values = data | rv::values;

            return (min ? ranges::min(values) : ranges::max(values));
        }

        // find min or max value < or <= x
        Value
        minOrMaxValueLessOrNotGreater(const Key x, const bool min, const bool strict) const noexcept {
            auto it = (strict ? data.lower_bound(x) : data.upper_bound(x));

            if (it != data.begin())
                --it;

            if (strict ? it->first >= x : it->first > x)
                return {};

            const auto subRange = ranges::subrange(data.begin(), ++it) | rv::values;

            return (min ? ranges::min(subRange) : ranges::max(subRange));
        }

        // find min or max value > or >= x
        Value
        minOrMaxValueGreaterOrNotLess(const Key x, const bool min, const bool strict) const noexcept {
            auto it = (strict ? data.upper_bound(x) : data.lower_bound(x));

            if (it == data.end())
                return {};

            const auto subRange = ranges::subrange(it, data.end()) | rv::values;

            return (min ? ranges::min(subRange) : ranges::max(subRange));
        }

        // find min or max value in (x, y) or [x, y]
        Value
        minOrMaxValueRange(const Key x, const Key y, const bool min, const bool strict) const noexcept {
            auto itA = (strict ? data.upper_bound(x) : data.lower_bound(x));
            auto itB = (strict ? data.lower_bound(y) : data.upper_bound(y));

            if (itB != data.begin())
                --itB;

            if (itA == data.end() || (strict ? itB->first >= y : itB->first > y))
                return {};

            const auto subRange = ranges::subrange(itA, ++itB) | rv::values;

            if (!subRange)
                return {};

            return (min ? ranges::min(subRange) : ranges::max(subRange));
        }

        // interpolate intermediate value
        Value interpolate(const Key i, const Key i0, const Key i1, const Value v0, const Value v1) const noexcept {
            return v0 + (i - i0) / (i1 - i0) * (v1 - v0); // i0 != i1
        }

    public:
        // access: [myObject].data.[method]
        // other methods: https://en.cppreference.com/w/cpp/container/map
        map<Key, Value> data;

        ipol() noexcept = default;

        // assign to map
        ipol(const map<Key, Value> &data_map) noexcept: data{data_map} {}

        // assign to map
        ipol &operator=(const map<Key, Value> &data_map) noexcept {
            data = data_map;
            return *this;
        }

        // add x to all
        ipol &operator+=(const Value x) noexcept {
            for (auto &value: data | rv::values)
                value += x;
        }

        // multiply all by x
        ipol &operator*=(const Value x) noexcept {
            for (auto &value: data | rv::values)
                value *= x;
        }

        // string representation
        string str() const noexcept {
            string s = "Interpolation data:\n";
            for (const auto &[key, value]: data)
                s += format("{}: {}\n", key, value);
            return s;
        }

        // print to a stream
        friend ostream &operator<<(ostream &os, const ipol &obj) {
            os << obj.str();
            return os;
        }

        // find min key
        Key minKey() const noexcept {
            return minOrMaxKey(true);
        }

        // find max key
        Key maxKey() const noexcept {
            return minOrMaxKey(false);
        }

        // find min value
        Value minValue() const noexcept {
            return minOrMaxValue(true);
        }

        // find max value
        Value maxValue() const noexcept {
            return minOrMaxValue(false);
        }

        // find min value <= x
        Value minValueNotGreater(const Key x) const noexcept {
            return minOrMaxValueLessOrNotGreater(x, true, false);
        }

        // find min value < x
        Value minValueLess(const Key x) const noexcept {
            return minOrMaxValueLessOrNotGreater(x, true, true);
        }

        // find max value <= x
        Value maxValueNotGreater(const Key x) const noexcept {
            return minOrMaxValueLessOrNotGreater(x, false, false);
        }

        // find max value < x
        Value maxValueLess(const Key x) const noexcept {
            return minOrMaxValueLessOrNotGreater(x, false, true);
        }

        // find min value >= x
        Value minValueNotLess(const Key x) const noexcept {
            return minOrMaxValueGreaterOrNotLess(x, true, false);
        }

        // find min value > x
        Value minValueGreater(const Key x) const noexcept {
            return minOrMaxValueGreaterOrNotLess(x, true, true);
        }

        // find max value >= x
        Value maxValueNotLess(const Key x) const noexcept {
            return minOrMaxValueGreaterOrNotLess(x, false, false);
        }

        // find max value > x
        Value maxValueGreater(const Key x) const noexcept {
            return minOrMaxValueGreaterOrNotLess(x, false, true);
        }

        // find min value in [x, y]
        Value minValueRangeNotStrict(const Key x, const Key y) const noexcept {
            return minOrMaxValueRange(x, y, true, false);
        }

        // find min value in (x, y)
        Value minValueRangeStrict(const Key x, const Key y) const noexcept {
            return minOrMaxValueRange(x, y, true, true);
        }

        // find max value in [x, y]
        Value maxValueRangeNotStrict(const Key x, const Key y) const noexcept {
            return minOrMaxValueRange(x, y, false, false);
        }

        // find max value in (x, y)
        Value maxValueRangeStrict(const Key x, const Key y) const noexcept {
            return minOrMaxValueRange(x, y, false, true);
        }

        // Returns true if the map is filled in with at least one non-zero value
        bool nonZero() const noexcept {
            return ranges::any_of(data | ranges::views::values, [](const auto x) -> bool { return x != 0; });
        }

        // interpolate i (better for pointers, default ())
        Value ip(const Key &i) const noexcept override {
            if (data.empty())
                return {};

            auto itUp = data.lower_bound(i);

            if (itUp == data.end())
                return data.rbegin()->second;

            if (itUp == data.begin())
                return data.begin()->second;

            auto itLo = prev(itUp);

            const auto &[i0, v0] = *itLo;
            const auto &[i1, v1] = *itUp;
            return interpolate(i, i0, i1, v0, v1);
        }

/*        Value operator()(const vector<Key> &vec) {
            const size_t regions = 1 << vec.size(); // 2 ^ size, cheaper, faster, more precise
            vector<vector<Value> > y(regions);
            vector<Value> dist(regions, -1);

            Value c;
            for (const auto &[key, value]: data) {

            }
        }*/
    };

    //Multidimensional interpolation
    template<floating_point Key, floating_point Value>
    class ipolm : public vipol<vector<Key>, Value> {
    private:
        //returns a minimal or maximal key
        vector<Value> minOrMaxKey(const bool min) const noexcept {
            if (data.empty())
                return vector<Key>{}; // returns empty vector!

            vector<Key> ret(data.begin()->first.size()); // every vector must be the same size! cheap construction

            if (min)
                for (size_t i = 0; i < ret.size(); ++i)
                    ret[i] = min_element(data.begin(), data.end(), [&](const auto &p) {return p.first[i]; })->first[i];
            else
                for (size_t i = 0; i < ret.size(); ++i)
                    ret[i] = max_element(data.begin(), data.end(), [&](const auto &p) { return p.first[i]; })->first[i];

            return ret;
        }

    public:
        map<vector<Key>, Value> data;

        ipolm() noexcept = default;

        // assign to map
        ipolm(const map<vector<Key>, Value> &data_map) noexcept: data{data_map} {}

        // assign to map
        ipolm &operator=(const map<vector<Key>, Value> &data_map) noexcept {
            data = data_map;
            return *this;
        }

        // add x to all
        ipolm &operator+=(const Value x) noexcept {
            for (auto &value: data | rv::values)
                value += x;
        }

        // multiply all by x
        ipolm &operator*=(const Value x) noexcept {
            for (auto &value: data | rv::values)
                value *= x;
        }

        //returns a minimal key
        vector<Value> minKey() const noexcept {
            return minOrMaxKey(true);
        };

        //returns a maximal key
        vector<Value> maxKey() const noexcept {
            return minOrMaxKey(false);
        };

        Value ip(const vector<Key> &i) const noexcept override {
            return {};
        }
    };
}

#endif