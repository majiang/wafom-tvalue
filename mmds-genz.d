module mmds;

alias uint U;
import std.stdio;
import lib.pointsettype;
import lib.integral;
import lib.wafom;
//import lib.testfunction;
import std.conv : to;
import std.string;
import std.math;
alias log2 lg;

import std.algorithm, std.functional;
alias sup = pipe!(map!abs, reduce!max, lg);
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

auto genzFactory(F)(size_t genzIndex, size_t dimensionR, in F coeff)
{
    import lib.testfunction : _gf = genzFactory;
    import std.functional : memoize;
    alias gf = memoize!(_gf!F);
    return gf(genzIndex, dimensionR, coeff);
}

void main(string[] args)
{
    if (args.length < 3)
        return "mmds-genz genzIndex numDS [coeff]".writeln();
    immutable
        genzIndex = args[1].to!size_t(),
        numDS = args[2].to!size_t(),
        coeff = (3 < args.length) ? args[3].to!real() : 1;
    foreach (line; stdin.byLine())
    {
        auto P = line.fromString!U();
        auto f = genzIndex.genzFactory(P.dimensionR, coeff);
        debug if (10 < P.dimensionF2) break;
        real[] e;
        foreach (i; 0..numDS)
            e ~= P.shiftRandomly().signedIntegrationError(f);
        auto dat = e.adjoin!(sup, rms);
        "%s,%.15f,%.15f,%.15f".writefln
        (line.strip(), P.bimswafom().lg(), dat[0], dat[1]);
    }
}
