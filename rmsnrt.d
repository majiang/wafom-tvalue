import lib.pointsettype;
import lib.wafom;
import std.conv;
import std.math : lg = log2;
import std.stdio;

auto ifany(T, T defaultValue)(T[] arr, size_t idx)
{
    if (arr.length <= idx)
        return defaultValue;
    return arr[idx];
}


void main(string[] args)
{
    if (args.length != 5)
        return stderr.writeln("precision(<=32) dimF2 dimR generate");
    auto buf = args[1..$];
    immutable precision = buf[0].to!size_t().min(32);
    immutable dimensionF2 = buf[1].to!size_t();
    immutable dimensionR = buf[2].to!size_t();
    immutable count = buf[3].to!size_t();
    foreach (i; 0..count)
    {
        auto P = randomPointSet!uint(Precision(precision), DimensionR(dimensionR), DimensionF2(dimensionF2));
        "%s,%.15f".writefln(P, P.bimsnrtwafom().lg());
    }
}
