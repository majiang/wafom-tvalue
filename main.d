module main;

import tvalue : tvalue1, tvalue2;
import wafom : wafom;
import sobol : defaultSobols;
import pointset : randomPoints;
import std.algorithm : min, max, reduce;

import std.stdio;
version = integral;
version (integral)
{
    import asianoption : default_integrand;
    import integral : integral;
    void main()
    {
        "description\tt-value\twafom\tintegral\tbasis".writeln();
        foreach (d; 16..33)
        {
            auto P = defaultSobols(4, d, d);
            "Sobol(%d)\t%d\t%f\t%.15f".writefln(
                d, P.save.tvalue1(), P.save.wafom(), integral!default_integrand(P));
        }
        version (none) foreach (d; 16..17)
        {
            foreach (t; 0..5000)
            {
                auto P = randomPoints(4, 32, d);
                "random(%d)\t%d\t%f\t%.15f\t%s".writefln(
                    d, P.save.tvalue1(), P.save.wafom(), integral!default_integrand(P), P.basis);
            }
        }
    }
}
else version (asianoption)
{
    void main()
    {
        auto P = defaultSobols(4, 16, 16);
//        P.integral_sobol();
    }
}
else version (unittest_only)
{
    void main()
    {
    }
}
else version (tvaluegeneral)
{
    void main()
    {
        size_t dimension = 2;
        size_t precision = 12;
        size_t lg_length = 10;
        foreach (i; 0..100)
        {
            randomPoints(dimension, precision, lg_length).tvalue2().writeln();
        }
    }
}
else version (arc)
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
        foreach (dimension; 1..11)
        {
            writeln("length.lg = ", length, "; dimension = ", dimension);
            auto sobol = defaultSobols(dimension, length, length);
            auto t1 = sobol.save.tvalue1();
            "algorithm 1: t-val = ".writeln(t1);
            auto t2 = sobol.save.tvalue2();
            "algorithm 2: t-val = ".writeln(t2);
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
