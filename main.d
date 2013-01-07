module main;

import tvalue : tvalue1, tvalue2;
import wafom : wafom;
import sobol : Sobol, Sobols, direction_numbers;

import std.stdio;

void main()
{
    foreach (length; [4, 6, 8, 10, 12, 14, 16]){
    auto DN =
    [
        direction_numbers([1], 3, length),
        direction_numbers([1, 1], 7, length),
        direction_numbers([1, 3, 7], 11, length),
        direction_numbers([1, 1, 5], 13, length),
        direction_numbers([1, 3, 1, 1], 19, length),
        direction_numbers([1, 1, 3, 7], 25, length),
        direction_numbers([1, 3, 3, 9, 9], 37, length),
        direction_numbers([1, 3, 7, 13, 3], 59, length),
        direction_numbers([1, 1, 5, 11, 27], 47, length),
        direction_numbers([1, 3, 5, 1, 15], 61, length),
    ];
    foreach (i; 0..10)
    {
        auto sobol = Sobols(DN[0..(i + 1)]);
        auto t1 = sobol.save.tvalue1();
        "t-val = ".writeln(t1, " by Algorithm 1");
        auto t2 = sobol.save.tvalue2();
        "t-val = ".writeln(t2, " by Algorithm 2");
        if (t1 != t2)
            "    Error!".writeln();
        "wafom = ".writeln(sobol.wafom());
        writeln();
    }
    }
}

version (array) Tuple!(ulong, double) measure
(
    ulong p, ulong e, ulong m, ulong s
)(
    F!(p, e)[m][m][s] C,
    F!(p, e) delegate (ulong) phi,
    ulong delegate (F!(p, e)) ihp
)
{
    auto net = DigitalNet!(p, e, m, s)(C, phi, ihp);
    return Tuple!(ulong, double)(tvalue!(p ^^ e, m, s)(net), wafom!(n, b, m, s));
}
