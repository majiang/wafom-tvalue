module main;
import std.stdio, std.array, std.string, std.conv;
import lib.pointset : randomVector, nonshiftedRandomBasisPoints;

void main()
{
    auto buf = readln().strip().split();
    immutable precision = buf[0].to!size_t();
    immutable dimensionF2 = buf[1].to!size_t();
    immutable dimensionR = buf[2].to!size_t();
    auto P = nonshiftedRandomBasisPoints!ubyte(precision, dimensionR, dimensionF2);
    auto v = randomVector!ubyte(precision, dimensionR);
    foreach (X; P)
    {
        foreach (x; X)
            x.write(",");
        writeln();
    }
    writeln();
    foreach (X; P * v)
    {
        foreach (x; X)
            x.write(",");
        writeln();
    }
    writeln();
}
