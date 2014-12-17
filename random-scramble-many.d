module main;
import std.stdio, std.array, std.string, std.conv;
import lib.pointsettype;
alias U = uint;

void main(string[] args)
{
    auto buf = args[1..$];
    if (buf.length != 4)
        return "random-scramble-many precision dimB dimR count".writeln();
    immutable precision = buf[0].to!size_t();
    immutable dimensionF2 = buf[1].to!size_t();
    immutable dimensionR = buf[2].to!size_t();
    immutable count = buf[3].to!size_t();
    auto P = source(Precision(precision), DimensionR(dimensionR), DimensionF2(dimensionF2));
    foreach (i; 0..count)
    	P.scrambleRandomly()[0].writeln();
}

auto source(Precision precision, DimensionR dimensionR, DimensionF2 dimensionF2)
{
	import std.algorithm;
	immutable filename = "nx-low-s%02d-all.csv".format(dimensionR.s);
	foreach (P; File(filename).byLine().map!(fromString!U)())
		if (P.precision == precision.n && P.dimensionR == dimensionR.s && P.dimensionF2 == dimensionF2.m)
			return P;
	assert (false, "no source point set to scramble: %s".format(filename));
}
