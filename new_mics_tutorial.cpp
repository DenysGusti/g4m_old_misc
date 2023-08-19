#include <iostream>
#include <chrono>

#include "new_misc.hpp"
#include "misc.h"

using namespace std;
using ld = long double;

using namespace g4m;

static uint64_t cc{};

void *operator new(const size_t size) {
    ++cc;
    cout << "Allocating " << size << " bytes" << endl;
    return malloc(size);
}

void operator delete(void *memory, const size_t size) noexcept {
    cout << "Freeing " << size << " bytes" << endl;
    free(memory);
}

void operator delete(void *memory) noexcept {
    free(memory);
}

class Timer {
private:
    chrono::high_resolution_clock::time_point startTimepoint;

public:
    Timer() : startTimepoint{chrono::high_resolution_clock::now()} {}

    void stop() {
        auto endTimepoint = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(endTimepoint - startTimepoint);
        cout << "elapsed time: " << duration << endl;
    }

    ~Timer() {
        stop();
    }
};

int main() {
    Ipol<ld> ipol;
    ipol.data = {{1, 2},
                 {3, 4},
                 {5, 6},
                 {7, 8},
                 {9, 10}};
    cout << ipol.str();

    ipol.data[11] = 12;
    ipol.data.insert({13, 14});
    ipol.data.emplace(15, 16);

    ipol.data.erase(15);
    erase_if(ipol.data,
             [](const auto &item) -> bool {
                 const auto &[key, value] = item;
                 return (key + value) > 25;
             });
    cout << ipol;

    cout << ipol.data.size() << '\n'
         << ipol.data.max_size() << '\n'
         << ipol.data.count(1) << '\n'
         << boolalpha << ipol.data.empty() << ' ' << ipol.data.contains(1) << ' '
         << (ipol.data.find(11) != ipol.data.end()) << ' ' << (ipol.data > map<ld, ld>{{1000, -1}}) << '\n';

    map tmp = {pair{1.L, 2.L}};
    swap(ipol.data, tmp);
    cout << ipol;
    swap(ipol.data, tmp);

    cout << ipol.data.lower_bound(5)->second << ' ' << ipol.data.upper_bound(5)->second << '\n'
         << ipol.data[9] << ' ' << ipol.data.at(7) << '\n';

    ipol.data[6] = ipol(6);
    cout << ipol
         << ipol.nonZero() << '\n';

    cout << ipol.minKey() << ' ' << ipol.maxKey() << '\n'
         << ipol.minValue() << ' ' << ipol.maxValue() << '\n'
         << ipol.maxValueGreater(5) << ' ' << ipol.maxValueNotLess(5) << '\n'
         << ipol.maxValueLess(5) << ' ' << ipol.maxValueNotGreater(5) << '\n'
         << ipol.minValueGreater(5) << ' ' << ipol.minValueNotLess(5) << '\n'
         << ipol.minValueLess(5) << ' ' << ipol.minValueNotGreater(5) << '\n'
         << ipol.maxValueRangeNotStrict(5, 7) << ' ' << ipol.maxValueRangeStrict(5, 7) << '\n'
         << ipol.minValueRangeNotStrict(-10, 10) << ' ' << ipol.minValueRangeStrict(8, 9) << '\n';

//    ipol.data.clear();

    IpolM<ld> test;
    test.data = {{{1, 2}, 3},
                 {{4, 5}, 6}};

    for (const auto &[k, v]: test.data)
        cout << k.size() << endl;

    cout << test;
    cout << test.minKey().front() << ' ' << test.maxKey().front() << endl;

    IpolM<ld> vd;
    vd.data[{0, 0, 0}] = 10;
    vd.data[{10, 5, 0}] = 20;
    cout << vd;

    cout << "Min:";
    for (const auto i: vd.minKey())
        cout << "\t" << i;
    cout << endl;
    cout << "Max:";
    for (const auto i: vd.maxKey())
        cout << "\t" << i;
    cout << endl;

    cout << vd({{10, 5, 0}}) << endl;
    cout << vd({{5, 2.5, 0}}) << endl;
    cout << vd({{7.5, 3.75, 0}}) << endl;
    cout << vd({{0, 0, 0}}) << endl;
    vd *= 2.5;
    cout << vd({{0, 0, 0}}) << endl;

    vd.data.clear();

    vd.data[{10, 10}] = 110;
    vd.data[{20, 10}] = 120;
    vd.data[{10, 20}] = 210;
    vd.data[{20, 30}] = 320;
    vd.data[{15, 30}] = 999;
    cout << vd;
    cout << "mip: " << vd({{15, 15}}) << endl;

    {
        Timer t;
        g4m::ipol<double, double> d;
        d.insert(0., 0.);
        d.insert(5., 10.);
        d.insert(15., 12.);
        d.insert(35., 13.4);
        d.insert(60., 16.2);
        d.insert(100., 20.);
        g4m::fipol<double> fd(d);

        Ipol<ld> dd;
        dd.data = {{0,   0},
                   {5,   10},
                   {15,  12},
                   {35,  13.4},
                   {60,  16.2},
                   {100, 20}};
        FIpol<ld> fdd{dd};
    }

    {
        g4m::fipol<double> fip(10);
        fip.fill(10.);
        fip.insert(4, 15.5);
        fip.insert(5, 25.5);
        fip.insert(6, 35.5);
        fip.insert(16, 35.5);
        cout << fip.g(3.) << endl;
        cout << fip.g(3.5) << endl;
        cout << fip.g(4.) << endl;
        cout << fip.g(4.5) << endl;
        cout << fip.g(5.) << endl;
        cout << fip.g(16.) << endl;
        cout << fip.g((unsigned int) 3) << endl;
        cout << fip.g((unsigned int) 4) << endl;
        cout << fip.g((unsigned int) 16) << endl;

        fip *= 2.5;
        cout << fip.g(4.5) << endl;
        cout << fip.g((unsigned int) 4) << endl;

        FIpol<ld> fip1;
        fip1.data = {10, 10, 10, 10, 15.5, 25.5, 35.5};
        cout << format("\n{}\n{}\n{}\n{}\n{}\n{}\n\n", fip1(3), fip1(3.5), fip1(4), fip1(4.5), fip1(5), fip1(16));
        fip1 *= 2.5;
        cout << fip1;
        cout << format("\n{}\n{}\n{}\n{}\n{}\n{}\n\n", fip1(3), fip1(3.5), fip1(4), fip1(4.5), fip1(5), fip1(16));
    }
    cout << cc << '\n';
    return 0;
}

/*
Min:	0	0	0
Max:	10	5	0
20
15
17.5
10
25
mip: 171.429
*/
