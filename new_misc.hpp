#ifndef G4M_NEW_MISC_H
#define G4M_NEW_MISC_H

#include <map>
#include <vector>

#include <cmath>
#include <algorithm>
#include <numeric>
#include <functional>
#include <ranges>

#include <iostream>
#include <string>
#include <format>

using namespace std;
namespace rv = ranges::views;

using ld = long double;

template<class T>
concept FloatingPointVector = same_as<T, vector<class T::value_type>> && floating_point<class T::value_type>;

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
        [[nodiscard]] virtual string str() const noexcept = 0;

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

        [[nodiscard]] string str() const noexcept override {
            string s = "Ipol data:\n";
            s.reserve(s.length() + 32 * data.size());
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
        [[nodiscard]] bool nonZero() const noexcept {
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

        [[nodiscard]] string str() const noexcept override {
            string s = "Ipolm data:\n";
            s.reserve(s.length() + 64 * data.size());
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
                        cerr << format("out of range problem: pos = {}, regions = {}", pos, regions)
                             << endl;  // or to logger
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

    // Fast Interpolation where the steps of the index are 1 and starting at 0
    template<floating_point T>
    class Fipol : public Vipol<T, T> {
    private:
        T intercept = 0, zoom = 1;
    public:
        vector<T> data;

        Fipol() = default;

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

        [[nodiscard]] string str() const noexcept override {
            string s = "Fipol data:\n";
            s.reserve(s.length() + 16 * data.size());
            for (const auto el: data)
                s += format("{}\n", el);
            return s;
        }

        T ip(const T &i) const noexcept override {
            if (data.empty())
                return 0;
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


    template<floating_point T>
    class Fipolm : public Vipol<vector<T>, T> {
    private:
        size_t dim = 0;
        vector<size_t> n;  // Size of array and how many dimensions does the index have
        vector<T> intercept; // If data range does not start from 0
        vector<T> zoom;

        void fillMap(const Ipolm<T, T> &t, vector<T> &key, const size_t idx_ = 0, const size_t adim = 0,
                     const size_t mul = 1) {
            for (size_t i = 0; i < n[adim]; ++i) {
                key[adim] = (intercept[adim] + i) / zoom[adim];
                if (adim + 1 < dim) {
                    fillMap(t, key, idx_ + i * mul, adim + 1, mul * n[adim]);
                } else {
                    data[idx_ + i * mul] = t(key);
                }
            }
        };

    public:
        vector<T> data;

        explicit Fipolm(const vector<size_t> &n_) : n{n_}, dim{n_.size()} {
            intercept.assign(dim, 0);
            zoom.assign(dim, 1);
            data.assign(accumulate(n.cbegin(), n.cend(), 1, multiplies<>{}), 0);
        }

        Fipolm(const Ipolm<T, T> &t, const vector<T> &zoom_) : zoom{zoom_} {
            vector<T> idxMin = t.minKey();
            vector<T> idxMax = t.maxKey();

            dim = idxMin.size();

            intercept.reserve(dim);
            for (size_t i = 0; i < dim; ++i)
                intercept.push_back(zoom[i] * idxMin[i]);

            n.reserve(dim);
            for (size_t i = 0; i < dim; ++i)
                n.push_back(1 + zoom[i] * (ceil(idxMax[i]) - floor(idxMin[i])));

            data.assign(accumulate(n.cbegin(), n.cend(), 1, multiplies<>{}), 0);

            fillMap(t, vector<T>(dim));
        }

        explicit Fipolm(const Ipolm<T, T> &t) : Fipolm(t, vector<T>{t.minKey().size(), 1}) {}

        bool insert(const vector<size_t> &indices, const T value) {
            if (indices.size() > dim) {
                throw invalid_argument{"indices size is bigger than dim"};
            }
            for (size_t i = 0; i < indices.size(); ++i)
                if (indices[i] >= n[i])
                    return false;

            size_t index = 0;
            size_t mul = 1;
            for (size_t i = 0; i < indices.size(); ++i) {
                index += indices[i] * mul;
                mul *= n[i];
            }
            data[index] = value;
            return true;
        }

        [[nodiscard]] vector<size_t> getN() const noexcept {
            return n;
        }

        // add x to all
        Fipolm &operator+=(const T x) noexcept override {
            for (auto &value: data)
                value += x;
            return *this;
        }

        // multiply all by x
        Fipolm &operator*=(const T x) noexcept override {
            for (auto &value: data)
                value *= x;
            return *this;
        }

        [[nodiscard]] string str() const noexcept override {
            string s = "Fipolm data:\n";
            s.reserve(s.length() + 16 * data.size());
            for (const auto el: data)
                s += format("{}\n", el);
            return s;
        }

        T ip(const vector<T> &i) const noexcept override {
            if (i.size() != dim)
                return 0;

            T k_j = clamp(i[0] * zoom[0] + intercept[0], 0., n[0] - 1.);  // TODO to check sign

            size_t sur = 1 << n.size();  // n Surrounding points 2 ^ dim
            vector<size_t> idx(sur, string::npos);  // Index for surrounding points
            vector<T> dist(sur);  // Distance of surrounding points

            idx[0] = floor(k_j);
            idx[1] = ceil(k_j);
            dist[0] = abs(k_j - idx[0]); // Manhattan distance
            dist[1] = abs(k_j - idx[1]);

            uint32_t mul = n[0];  // Array size in dimension 0

            size_t t{}, uc{}, uf{};
            T dc{}, df{};
            for (size_t j = 1; j < dim; ++j) {
                k_j = clamp(i[j] * zoom[j] + intercept[j], 0., n[j] - 1.);  // TODO to check sign
                t = 1 << j;  // 2^n points used in this dim
                uc = ceil(k_j) * mul; // Index where of grid point
                uf = floor(k_j) * mul;
                dc = abs(k_j - ceil(k_j)); // Manhattan distance
                df = abs(k_j - floor(k_j));
                for (size_t k = 0; k < t; ++k) {
                    idx[k + t] = idx[k] + uc;
                    idx[k] += uf;
                    dist[k + t] = dist[k] + dc;
                    dist[k] += df;
                }
                mul *= n[j];
            }

            T sdist = 0;
            T sval = 0;
            for (size_t j = 0; j < sur; ++j)
                if (idx[j] >= 0) {
                    if (dist[j] > 0) {
                        sval += data[idx[j]] / dist[j];
                        sdist += 1 / dist[j];
                    } else {
                        sval = data[idx[j]];
                        sdist = 1;
                        break;
                    }
                }
            return sdist > 0 ? sval / sdist : 0;
        }
    };


    template<floating_point T>
    class FFipol : public Vipol<T, T> {
    private:
        T intercept = 0, zoom = 1;
    public:
        vector<T> data;

        explicit FFipol(const Ipol<T, T> &ipol, const T zoom_ = 1, const T add = 0.5) : zoom{zoom_}, intercept{
                zoom_ * ipol.minKey() + add} {
            const size_t n = 1 + zoom * (ceil(ipol.maxKey()) - floor(ipol.minKey()));
            data.assign(n, 0);
            for (size_t i = 0; i < data.size(); ++i)
                data[i] = ipol((intercept + i) / zoom - add);
        };

        // add x to all
        FFipol &operator+=(const T x) noexcept override {
            for (auto &value: data)
                value += x;
            return *this;
        }

        // multiply all by x
        FFipol &operator*=(const T x) noexcept override {
            for (auto &value: data)
                value *= x;
            return *this;
        }

        [[nodiscard]] string str() const noexcept override {
            string s = "FFipol data:\n";
            s.reserve(s.length() + 16 * data.size());
            for (const auto el: data)
                s += format("{}\n", el);
            return s;
        }

        T ip(const T &i) const noexcept override {
            if (data.empty())
                return 0;
            if (data.size() == 1)
                return data.front();

            const T zi = i * zoom + intercept;
            if (zi < 0)
                return data.front();
            if (zi >= data.size() - 1)
                return data.back();

            const size_t i0 = zi;
            return data[i0];
        }
    };

    template<floating_point T>
    class FFipolm : public Vipol<vector<T>, T> {
    private:
        size_t dim = 0;
        vector<size_t> n;  // Size of array and how many dimensions does the index have
        vector<T> intercept; // If data range does not start from 0
        vector<T> zoom;


        void fillMap(const Ipolm<T, T> &t, vector<T> &key, const size_t idx_, const size_t adim, const size_t mul,
                     const T add) {
            for (size_t i = 0; i < n[adim]; ++i) {
                key[adim] = (intercept[adim] + i) / zoom[adim] + add;
                if (adim + 1 < dim) {
                    fillMap(key, idx_ + i * mul, adim + 1, mul * n[adim], add);
                } else {
                    data[idx_ + i * mul] = t(key);
                }
            }
        };

    public:
        vector<T> data;

        explicit FFipolm(const vector<size_t> &n_) : n{n_}, dim{n_.size()} {
            intercept.assign(dim, 0);
            zoom.assign(dim, 1);
            data.assign(accumulate(n.cbegin(), n.cend(), 1, multiplies<>{}), 0);
        }

        FFipolm(const Ipolm<T, T> &t, const vector<T> &zoom_, const T add = 0.5) : zoom{zoom_} {
            vector<T> idxMin = t.minKey();
            vector<T> idxMax = t.maxKey();

            dim = idxMin.size();

            intercept.reserve(dim);
            for (size_t i = 0; i < dim; ++i)
                intercept.push_back(zoom[i] * idxMin[i]);

            n.reserve(dim);
            for (size_t i = 0; i < dim; ++i)
                n.push_back(1 + zoom[i] * (ceil(idxMax[i]) - floor(idxMin[i])));

            data.assign(accumulate(n.cbegin(), n.cend(), 1, multiplies<>{}), 0);

            vector<T> key(dim, 0);
            fillMap(t, key, 0, 0, 1, add);
        }

        explicit FFipolm(const Ipolm<T, T> &t, const T add = 0.5) : FFipolm(t, vector<T>{t.minKey().size(), 1}, add) {}

        bool insert(const vector<size_t> &indices, const T value) {
            if (indices.size() != dim) {
                throw invalid_argument{"indices size doesn't equal dim"};
            }
            for (size_t i = 0; i < indices.size(); ++i)
                if (indices[i] >= n[i])
                    return false;

            size_t index = 0;
            size_t mul = 1;
            for (size_t i = 0; i < indices.size(); ++i) {
                index += indices[i] * mul;
                mul *= n[i];
            }
            data[index] = value;
            return true;
        }

        [[nodiscard]] vector<size_t> getN() const noexcept {
            return n;
        }

        // add x to all
        FFipolm &operator+=(const T x) noexcept override {
            for (auto &value: data)
                value += x;
            return *this;
        }

        // multiply all by x
        FFipolm &operator*=(const T x) noexcept override {
            for (auto &value: data)
                value *= x;
            return *this;
        }

        [[nodiscard]] string str() const noexcept override {
            string s = "FFipolm data:\n";
            s.reserve(s.length() + 16 * data.size());
            for (const auto el: data)
                s += format("{}\n", el);
            return s;
        }

        T ip(const vector<T> &i) const noexcept override {
            if (i.size() != dim)
                return 0;

            size_t idx = 0;
            size_t mul = 1;
            size_t k = 0;
            for (size_t j = 0; j < dim; ++j) {
                k = clamp(static_cast<size_t>(i[j] * zoom[j] + intercept[j]), 0uz, n[j] - 1uz);
                idx += k * mul;
                mul *= n[j];
            }
            return data[idx];
        }
    };
}

#endif