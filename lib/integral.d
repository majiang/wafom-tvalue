module lib.integral;

import std.stdio;
import lib.pointsettype : Bisectable, toReals;
import std.traits : ParameterTypeTuple;
import std.range : ElementType;

private static flist = function (size_t n)
{
    real ret[];
    real cur = 1.0;
    foreach (i; 0..n)
    {
        ret ~= cur;
        cur *= 0.5;
    }
    return ret ~ cur;
}(65);

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
real integral(alias f, R)(R P)// if (isUnsigned!T)
{
    real result = 0;
    foreach (x; P.toReals!real())
        result += f(x);
    return result * flist[P.dimensionF2];
}

real integralNoncentering(alias f, R)(R P)// if (isUnsigned!T)
{
    real result = 0;
    foreach (x; P.toReals!real())
        result += f(x.revertCentering(P.precision + 1));
    return result * flist[P.dimensionF2];
}
real[] revertCentering(in real[] x, size_t n)
{
    real[] ret;
    foreach (e; x)
        ret ~= e - flist[n];
    return ret;
}

/// ditto
real bintegral(alias f, R)(R P) if (Bisectable!R)
{
    if (P.bisectable)
    {
        auto Q = P.bisect();
        return (Q[0].bintegral!(f, R)() + Q[1].bintegral!(f, R)()) * 0.5;
    }
    return P.integral!(f, R)();
}

real bintegralNoncentering(alias f, R)(R P) if (Bisectable!R)
{
    if (P.bisectable)
    {
        auto Q = P.bisect();
        return (Q[0].bintegralNoncentering!(f, R)() + Q[1].bintegralNoncentering!(f, R)()) * 0.5;
    }
    return P.integralNoncentering!(f, R)();
}
