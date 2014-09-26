module main;
import std.stdio, std.array, std.string, std.conv;
import lib.pointsettype;

/** Short sample of the library.

*/

void main(string[] args)
{
    auto buf = args[1..$];
    if (len(buf) != 4)
        return "random-generate-many precision dimB dimR count".writeln();
    immutable precision = buf[0].to!size_t();
    immutable dimensionF2 = buf[1].to!size_t();
    immutable dimensionR = buf[2].to!size_t();
    immutable count = buf[3].to!size_t();
    foreach (i; 0..count)
        randomPointSet!uint(Precision(precision), DimensionR(dimensionR), DimensionF2(dimensionF2)).writeln();
}
