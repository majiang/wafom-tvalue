import std.stdio;
import std.traits : isUnsigned;
import pointset : Bisectable;

private static flist = function (size_t n)
{
    double ret[];
    double cur = 1.0;
    foreach (i; 0..n)
    {
        ret ~= cur;
        cur *= 0.5;
    }
    return ret ~ cur;
}(33);

/** Perform numerical integration of a function f: [0..1)<sup>s</sup> -> R by the point set P.

s is implicitly given as P.dimension. f must be callable with double[] and return a double.

Usage:
----------------
immutable size_t dimension;
double f(double[] x)
{
    assert (x.length == dimension);
    //...
    return something;
}
double result = integral!f(randomPoints(dimension, precision, lg_length));
----------------
*/
double integral(alias f, R)(R P)// if (isUnsigned!T)
{
    double factor = flist[P.precision];
    double shift = flist[P.precision + 1];
    double result = 0;
    foreach (x; P)
    {
        result += f(x.affine(factor, shift));
    }
    return result * flist[P.dimensionF2];
}

/// ditto
double bintegral(alias f, R)(R P) if (Bisectable!R)
{
    if (P.bisectable)
    {
        auto Q = P.bisect();
        return (Q[0].bintegral!(f, R) + Q[1].bintegral!(f, R)) * 0.5;
    }
    return P.integral!(f, R)();
}

private double[] affine(T)(T[] x, double factor, double shift) if (isUnsigned!T)
{
    auto ret = new double[x.length];
    //ret[] = x[] * factor + shift;
    foreach (i, c; x)
    {
        ret[i] = c * factor + shift;
    }
    return ret;
}
