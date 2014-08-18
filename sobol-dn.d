module sobol;
import std.stdio, std.conv, std.algorithm;
import lib.pointsettype, lib.sobol;

void main(string[] args)
{
	if (args.length < 3)
		return "sobol-dn dimR dimF2...".writeln();
	immutable dimR = args[1].to!size_t();
	foreach (dimB; args[2..$].map!(to!size_t)())
	defaultSobols!uint(DimensionR(dimR), Precision(32), DimensionF2(dimB)).writeln();
}
