import std.stdio;
import pointset : BasisPoints;

/** Perform integration of a function f: [0..1)<sup><var>s</var></sup> -> R by the point set P.

<var>s</var> is implicitly given as P.dimension.
f must support double opCall(double[] x)).

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
double integral(alias f)(BasisPoints P)
{
    double factor = 1.0 / (1UL << P.precision);
    double shift = 0.5 * factor;
    double result = 0;
    foreach (x; P)
    {
        result += f(x.affine(factor, shift));
    }
    return result / P.length;
}

private double[] affine(ulong[] x, double factor, double shift)
{
    auto ret = new double[x.length];
    foreach (i, c; x)
    {
        ret[i] = c * factor + shift;
    }
    return ret;
}
