module lib.integral;

import std.stdio;
import lib.pointsettype : Bisectable, toReals;
import lib.testfunction;
import std.traits : ParameterTypeTuple;
import std.range : ElementType;
import std.typecons : Flag;

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


F signedIntegrationError(Flag!"centering" centering = Flag!"centering".yes, F, R)(R P, TestFunction!F f)
{
    return P.integral!(centering, F, R)(f) - f.I;
}

F integral(Flag!"centering" centering = Flag!"centering".yes, F, R)(R P, Function!F f)
    if (Bisectable!R)
{
    if (P.bisectable)
    {
        auto Q = P.bisect();
        return (Q[0].integral!centering(f) + Q[1].integral!centering(f)) * 0.5;
    }
    F result = 0;
    foreach (x; P.toReals!(F, centering)())
        result += f(x);
    return result * flist[P.dimensionF2];
}

auto integrals(Flag!"centering" centering = Flag!"centering".yes, F, R)(R Ps, Function!F f)
{
    return Ps.map!(P => P.integral!centering(f))();
}

auto integrationErrors(F, R)(R Ps, TestFunction!F f)
{
    return Ps.map!(P => P.signedIntegrationError!centering(f))();
}

import std.algorithm : map, reduce, min, max;
import std.conv : to;
alias minimum = reduce!min;
alias maximum = reduce!max;
//alias maxabs = reduce!((a, b) => a.abs().max(b.abs()))();

F stdDev(F, R)(R xs)
{
    auto xm = xs.mean!F();
    return xs.map!(x => x - xm).rootMeanSquare!F();
}

F mean(F, R)(R xs)
{
    return xs.sum!F() / xs.length;
}

private F meanPositive(F, R)(R xs)
{
    return xs.sumPositive!F() / xs.length;
}

F rootMeanSquare(F, R)(R xs)
{
    import std.math : sqrt;
    return xs.map!(a => a * a)().meanPositive!F().sqrt();
}

F sum(F, R)(R xs)
{
    return xs.map!(to!F)().reduce!((a, b) => a + b)();
}

private alias sumPositive = sum;

version (stand_alone_int):
final class Tf : TestFunction!real
{
    this (size_t dimensionR)
    {
        super (dimensionR * 0.5);
    }
    real opCall(in real[] x)
    {
        return x.reduce!((a, b) => a + b)();
    }
}

void main(string[] args)
{
    import std.conv : to;
    import lib.pointsettype;
    immutable dimR = 2;
    immutable dimB = args[1].to!size_t();
    auto P = randomPointSet!ubyte(Precision(8), DimensionR(dimR), DimensionF2(dimB));
    auto f = new Tf(dimR);
    "centering = Y: ".write(); P.signedIntegrationError(f).writeln();
    "centering = N: ".write(); P.signedIntegrationError!(Flag!"centering".no)(f).writeln();
}
