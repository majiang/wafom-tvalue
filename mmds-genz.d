module mmds;

alias uint U;
import std.stdio;
import lib.pointsettype;
import lib.integral;
import lib.wafom;
import lib.testfunction;
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

void main(string[] args)
{
    immutable numDS = args[2].to!size_t();
    auto f = (args.length < 4) ?
        args[1].to!size_t().genzFactory!real(4) :
        args[1].to!size_t().genzFactory!real(4, args[3].to!real)
    ;
    stderr.writeln(f.I);
    foreach (line; stdin.byLine())
    {
        auto P = line.fromString!U();
        debug if (10 < P.dimensionF2) break;
        real[] e;
        foreach (i; 0..numDS)
            e ~= P.shiftRandomly().signedIntegrationError(f);
        auto dat = e.adjoin!(sup, rms);
        "%s,%.15f,%.15f,%.15f".writefln
        (line.strip(), P.bimswafom().lg(), dat[0], dat[1]);
    }
}
