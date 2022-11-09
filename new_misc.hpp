#ifndef G4M_NEW_MISC_H
#define G4M_NEW_MISC_H

#include <map>
#include <vector>

#include <cmath>
#include <algorithm>
#include <functional>
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
    class Vipol {
    public:
        virtual ~Vipol() = default;

        virtual VAL ip(const IDX &i) const noexcept = 0;

        VAL operator()(const IDX &i) const noexcept {
            return ip(i);
        };

        // string representation
        virtual string str() const noexcept = 0;

        // print to a stream
        friend ostream &operator<<(ostream &os, const Vipol &obj) {
            os << obj.str();
            return os;
        }

        virtual Vipol &operator+=(VAL) noexcept = 0;

        virtual Vipol &operator*=(VAL) noexcept = 0;
    };

    template<floating_point Key, floating_point Value>
    class Ipol : public Vipol<Key, Value> {
    private:
        // find min or max key
        Key minOrMaxKey(const bool min_flag) const noexcept {
            if (data.empty())
                return {};

            return (min_flag ? data.begin()->first : data.rbegin()->first);
        }

        // find min or max value
        Value minOrMaxValue(const bool min_flag) const noexcept {
            if (data.empty())
                return {};

            const auto values = data | rv::values;

            return (min_flag ? ranges::min(values) : ranges::max(values));
        }

        // find min or max value < or <= x
        Value
        minOrMaxValueLessOrNotGreater(const Key x, const bool min_flag, const bool strict) const noexcept {
            auto it = (strict ? data.lower_bound(x) : data.upper_bound(x));

            if (it != data.begin())
                --it;

            if (strict ? it->first >= x : it->first > x)
                return {};

            const auto subRange = ranges::subrange(data.begin(), ++it) | rv::values;

            return (min_flag ? ranges::min(subRange) : ranges::max(subRange));
        }

        // find min or max value > or >= x
        Value
        minOrMaxValueGreaterOrNotLess(const Key x, const bool min_flag, const bool strict) const noexcept {
            auto it = (strict ? data.upper_bound(x) : data.lower_bound(x));

            if (it == data.end())
                return {};

            const auto subRange = ranges::subrange(it, data.end()) | rv::values;

            return (min_flag ? ranges::min(subRange) : ranges::max(subRange));
        }

        // find min or max value in (x, y) or [x, y]
        Value
        minOrMaxValueRange(const Key x, const Key y, const bool min_flag, const bool strict) const noexcept {
            auto itA = (strict ? data.upper_bound(x) : data.lower_bound(x));
            auto itB = (strict ? data.lower_bound(y) : data.upper_bound(y));

            if (itB != data.begin())
                --itB;

            if (itA == data.end() || (strict ? itB->first >= y : itB->first > y))
                return {};

            const auto subRange = ranges::subrange(itA, ++itB) | rv::values;

            if (!subRange)
                return {};

            return (min_flag ? ranges::min(subRange) : ranges::max(subRange));
        }

    public:
        // access: [myObject].data.[method]
        // other methods: https://en.cppreference.com/w/cpp/container/map
        map<Key, Value> data;

        Ipol() = default;

        // assign to map
        Ipol(const map<Key, Value> &data_map) noexcept: data{data_map} {}

        // assign to map
        Ipol &operator=(const map<Key, Value> &data_map) noexcept {
            data = data_map;
            return *this;
        }

        // add x to all
        Ipol &operator+=(const Value x) noexcept override {
            for (auto &value: data | rv::values)
                value += x;
            return *this;
        }

        // multiply all by x
        Ipol &operator*=(const Value x) noexcept override {
            for (auto &value: data | rv::values)
                value *= x;
            return *this;
        }

        string str() const noexcept override {
            string s = "Ipol data:\n";
            for (const auto &[key, value]: data)
                s += format("{}: {}\n", key, value);
            return s;
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
            return ranges::any_of(data | rv::values, [](const auto x) -> bool { return x != 0; });
        }

        // interpolate i (better for pointers, default ())
        Value ip(const Key &i) const noexcept override {
            if (data.empty())
                return {};

            const auto itUp = data.lower_bound(i);

            if (itUp == data.end())
                return data.rbegin()->second;

            if (itUp == data.begin())
                return data.begin()->second;

            const auto itLo = prev(itUp);

            const auto &[i0, v0] = *itLo;
            const auto &[i1, v1] = *itUp;
            return lerp(v0, v1, (i - i0) / (i1 - i0));  // interpolate intermediate value
        }
    };

    //Multidimensional interpolation
    template<floating_point Key, floating_point Value>
    class Ipolm : public Vipol<vector<Key>, Value> {
    private:
        //returns a minimal or maximal key
        vector<Key> minOrMaxKey(const bool min_flag) const noexcept {
            if (data.empty())
                return vector<Key>{}; // returns empty vector!

            // copy assignment of the first vector
            vector<Key> ret = data.begin()->first; // every vector must be the same size!

            if (min_flag)
                for (size_t i = 0; i < ret.size(); ++i)
                    for (const auto &key: data | rv::keys)
                        ret[i] = min(key[i], ret[i]);
            else
                for (size_t i = 0; i < ret.size(); ++i)
                    for (const auto &key: data | rv::keys)
                        ret[i] = max(key[i], ret[i]);

            return ret;
        }

    public:
        map<vector<Key>, Value> data;

        Ipolm() = default;

        // assign to map
        Ipolm(const map<vector<Key>, Value> &data_map) noexcept: data{data_map} {}

        // assign to map
        Ipolm &operator=(const map<vector<Key>, Value> &data_map) noexcept {
            data = data_map;
            return *this;
        }

        // add x to all
        Ipolm &operator+=(const Value x) noexcept override {
            for (auto &value: data | rv::values)
                value += x;
            return *this;
        }

        // multiply all by x
        Ipolm &operator*=(const Value x) noexcept override {
            for (auto &value: data | rv::values)
                value *= x;
            return *this;
        }

        string str() const noexcept override {
            string s = "Ipolm data:\n";
            for (const auto &[key, value]: data) {
                for (const auto el: key)
                    s += format("{:>2} ", el);
                s += format(":\t{}\n", value);
            }
            return s;
        }

        //returns a minimal key
        vector<Key> minKey() const noexcept {
            return minOrMaxKey(true);
        };

        //returns a maximal key
        vector<Key> maxKey() const noexcept {
            return minOrMaxKey(false);
        };

        Value ip(const vector<Key> &vec) const noexcept override {
            const size_t regions = 1 << vec.size(); // 2 ^ i.size(), but better
            // values of near points & shortest distance in region
            vector<pair<vector<Value>, Value> > y_dist(regions, {{}, -1});

            Value d{}; // distance
            size_t pos{}; // memory where to save the minimum distance
            size_t mul{}; // multiplier to get right pos
            for (Value tmp{}; const auto &[key, value]: data) {
                d = 0;
                pos = 0;
                mul = 1;
                for (size_t i = 0; i < key.size(); ++i) {
                    tmp = key[i] - vec[i];

                    if (tmp > 0)
                        pos += mul;

                    if (pos >= regions) {
                        cerr << format("out of range problem: pos = {}, regions = {}", pos, regions) << endl;
                        return Value{};
                    }

                    d += abs(tmp); // d += abs(tmp) - manhattan distance; d += tmp * tmp - geometric interpolation
                    mul <<= 1; // mul *= 2
                }
                if (auto &[y, dist] = y_dist[pos]; d <= dist || dist < 0) {
                    if (dist != d) {
                        y.clear();
                        dist = d;
                    }
                    y.push_back(value);
                }
            }

            Value ip = 0, distSum = 0;
            int64_t n = 0;
            for (const auto &[y, dist]: y_dist)
                if (dist > 0 && n == 0)
                    for (const auto y_el: y) {
                        ip += y_el / dist;
                        distSum += 1 / dist;
                    }
                else if (dist == 0) {
                    if (n == 0)
                        ip = 0;
                    for (const auto y_el: y) {
                        ip += y_el;
                        ++n;
                    }
                }

            if (n > 0)
                ip /= n;
            else if (distSum > 0)
                ip /= distSum;
            else
                ip = 0; // numeric_limits<Value>::quiet_NaN();
            return ip;
        }
    };

    template<floating_point T>
    class Fipol : public Vipol<T, T> {
    private:
        T intercept{}, zoom{};
    public:
        vector<T> data;

        Fipol() = default;

        Fipol(const vector<T> &data_vec) noexcept: data{data_vec} {}

        Fipol &operator=(const vector<T> &data_vec) noexcept {
            data = data_vec;
            return *this;
        }

        explicit Fipol(const Ipol<T, T> &ipol, const T zoom_ = 1) : zoom{zoom_}, intercept{zoom_ * ipol.minKey()} {
            const size_t n = 1 + zoom * (ceil(ipol.maxKey()) - floor(ipol.minKey()));
            data.assign(n, 0);
            for (size_t i = 0; i < data.size(); ++i)
                data[i] = ipol((intercept + i) / zoom);
        };

        // add x to all
        Fipol &operator+=(const T x) noexcept override {
            for (auto &value: data)
                value += x;
            return *this;
        }

        // multiply all by x
        Fipol &operator*=(const T x) noexcept override {
            for (auto &value: data)
                value *= x;
            return *this;
        }

        string str() const noexcept override {
            string s = "Fipol data:\n";
            for (const auto el: data)
                s += format("{}\n", el);
            return s;
        }

        T ip(const T &i) const noexcept override {
            if (data.empty())
                return T{};
            if (data.size() == 1)
                return data.front();

            const T zi = i * zoom + intercept;
            if (zi < 0)
                return data.front();
            if (zi >= data.size() - 1)
                return data.back();

            const size_t i0 = zi;
            return lerp(data[i0], data[i0 + 1], zi - i0);  // linear interpolation function
        }
    };
}

#endif