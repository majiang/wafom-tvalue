module mmdsmc;

alias uint U;
import std.stdio;
import lib.integral;
//import lib.testfunction;
import std.conv : to;
import std.string;
import std.math;
alias log2 lg;
import std.random;

import std.algorithm, std.functional;

auto genzFactory(F)(size_t genzIndex, size_t dimensionR, in F coeff)
{
    import lib.testfunction : _gf = genzFactory;
    import std.functional : memoize;
    alias gf = memoize!(_gf!F);
    return gf(genzIndex, dimensionR, coeff);
}

void main(string[] args)
{
    import lib.montecarlo, lib.pointsettype;
    if (args.length < 5)
        return "mmds-genz-mc genzIndex dimR [dimBmin )dimBmax [coeff]".writeln();
    immutable
        genzIndex = args[1].to!size_t(),
        dimR = args[2].to!size_t(),
        dimBmin = args[3].to!size_t(),
        dimBmax = args[4].to!size_t(),
        coeff = (5 < args.length) ? args[5].to!real() : 1;
    auto f = genzIndex.genzFactory(dimR, coeff);
    real[] sample;
    foreach (x; UniformPoints!()(DimensionR(dimR), DimensionF2(dimBmax)).toReals!real())
    //foreach (x; MonteCarloPoints!uint(Precision(32), DimensionR(dimR), DimensionF2(dimBmax)).toReals!real())
        sample ~= f(x);
    immutable v = sample.stdDev!real().lg();
    foreach (dimB; dimBmin..dimBmax)
        "%d,%.15f".writefln(dimB, v - dimB * 0.5);
}
