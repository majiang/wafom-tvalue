module main;

import tvalue : tvalue1, tvalue2;
import wafom : wafom;
import sobol : sobols, direction_numbers, defaultSobols;
import pointset;
import std.algorithm : min, max, reduce;

import std.stdio;
version = random_search_wafom;
version (arc)
{
    void main()
    {
        foreach (x; defaultSobols(2, 10, 10))
        {
            writeln(x, ",");
        }
    }
}
else version (tmp) 
{
import integral : integral;
immutable size_t dimension = 2;
double f(double[] x)
{
    assert (x.length == dimension);
    foreach (e; x)
    {
        assert (0 < e && e < 1);
    }
    return x[0] ^^ 2 + x[1] ^^ 3;
}
void main()
{
    immutable size_t precision = 5;
    immutable size_t lg_length = 10;
    "integral[x^2+y^3] = %f".writefln(
        integral!f(randomPoints(dimension, precision, lg_length)));
}
}
else version (random_search_wafom) void main()
{
    double[] wafoms;
    immutable n = 1000, dimension = 4, precision = 16, lg_length = 10;

    foreach (i; 0..n)
    {
        wafoms ~= dimension.randomPoints(precision, lg_length).wafom();
    }
    "after %d times of iteration, min-wafom = %f and av-wafom = %f for %d-dimensional %d-bit point set.".writefln(
        n, wafoms.reduce!min(), wafoms.reduce!q{a+b} / n, dimension, precision
    );
    "full data: %s".writefln(wafoms);
}
else version (tvaluedebug) void main()
{
    foreach (length; [4, 6, 8, 10, 12, 14, 16]) // m
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
        foreach (i; 0..10)
        {
            writeln("length = ", length, "; dimension = ", i+1);
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
else void main()
{
}
