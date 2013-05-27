import std.stdio;
import std.traits : isUnsigned;
import pointset : Bisectable;

static flist = function (size_t n)
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

/** Perform integration of a function f: [0..1)<sup>s</sup> -> R by the point set P.

* s is implicitly given as P.dimension. f must support double opCall(double[] x)).

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
double integral(alias f, T, R)(R P) if (isUnsigned!T)
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

double bintegral(alias f, T, R)(R P) if (isUnsigned!T && Bisectable!R)
{
    if (P.bisectable)
    {
        auto Q = P.bisect();
        return (Q[0].bintegral!(f, T, R) + Q[1].bintegral!(f, T, R)) * 0.5;
    }
    return P.integral!(f, T, R);
}

private double[] affine(T)(T[] x, double factor, double shift) if (isUnsigned!T)
{
    auto ret = new double[x.length];
    foreach (i, c; x)
    {
        ret[i] = c * factor + shift;
    }
    return ret;
}
