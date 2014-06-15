alias uint U;
import std.stdio;
import lib.pointsettype;
import lib.integral;
import lib.wafom;
import std.conv : to;
import std.string;
import std.math;
alias log2 lg;

auto f(in real[] x)
{
    import std.algorithm;
    return (x.reduce!((a, b) => a + b)() * -2).exp();
}

real I(size_t dimension)
{
    immutable I1 = (-2).expm1() / -2;
    real s = dimension;
    return I1 ^^ s;
}

auto err(T)(T P)
{
    return (P.integral!f() - P.dimensionR.I()).abs();
}

import std.algorithm, std.functional;
alias sup = pipe!(reduce!max, lg);
auto rms(T)(T xs)
{
    size_t count;
    real sumsq = 0;
    foreach (x; xs)
    {
        sumsq += x * x;
        count += 1;
    }
    return (sumsq / count).lg() / 2;
}

void main(string[] args)
{
    immutable numDS = args[1].to!size_t();
    foreach (line; stdin.byLine())
    {
        auto P = line.fromString!U();
        real[] e;
        foreach (i; 0..numDS)
            e ~= P.shiftRandomly().err();
        auto dat = e.adjoin!(sup, rms);
        "%s,%.15f,%.15f,%.15f".writefln(line.strip(), P.bimswafom().lg(), dat[0], dat[1]);
    }
}
