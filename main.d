module main;

import tvalue : tvalue, tvalue1;
import wafom : wafom;
import sobol : Sobol, Sobols, direction_numbers;

import std.stdio;

const length = 6;

void main()
{
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
    //version (none)
    foreach (i; 0..3)
    {
        auto sobol = Sobols(DN[0..(i + 1)]);
        string sep = "";
        foreach (x; sobol.save)
        {
            foreach (e; x) e.write(",");
            writeln();
        }
        writeln();
        "t-val = ".writeln(sobol.save.tvalue1(), " by Algorithm 1");
        "t-val = ".writeln(sobol.save.tvalue(), " by Algorithm 2");
        "wafom = ".writeln(sobol.wafom());
        writeln();
    }
    version (none)
    foreach (x; Sobol(direction_numbers([1], 2, 6)))
    {
        x.write(", ");
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
