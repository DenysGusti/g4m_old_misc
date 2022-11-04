#include <iostream>

#include "new_misc.hpp"
#include "misc.h"

using namespace std;
using ld = long double;

int main() {
    g4m::Ipol<ld, ld> ipol;
    ipol = {{1, 2},
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

    ipol.data.clear();

    g4m::Ipolm<ld, ld> test = {{{{1, 2}, 3}, {{4, 5}, 6}}};
    for (const auto &[k, v]: test.data)
        cout << k.size() << endl;

    cout << test;
    cout << test.minKey().front() << ' ' << test.maxKey().front() << endl;

        g4m::Ipolm<ld, ld> vd;
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

        cout << vd({10, 5, 0}) << endl;
        cout << vd({5, 2.5, 0}) << endl;
        cout << vd({7.5, 3.75, 0}) << endl;
        cout << vd({0, 0, 0}) << endl;
        vd *= 2.5;
        cout << vd({0, 0, 0}) << endl;

        vd.data.clear();

        vd.data[{10, 10}] = 110;
        vd.data[{20, 10}] = 120;
        vd.data[{10, 20}] = 210;
        vd.data[{20, 30}] = 320;
        vd.data[{15, 30}] = 999;
        cout << vd;
        cout << "mip: " << vd({15, 15}) << endl;
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
