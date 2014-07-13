import annealing;
import lib.wafom, lib.pointsettype;
alias U = uint;

import std.stdio;
import std.math : lg = log2;

auto second(alias f, T)(T x)
{
	return tuple(x[0], f(x[1]));
}

void main(string[] args)
{
	if (args.length < 7)
		return usage();
	auto dimR = args[1].to!size_t();
	auto dimBmin = args[2].to!size_t();
	auto dimBmax = args[3].to!size_t();
	auto step = args[4].to!size_t();
	auto Tinit = args[5].to!double();
	auto Tend = args[6].to!double();
	size_t duplication = 1;
	if (args.length == 8)
		duplication = args[7].to!size_t();
	foreach(dimB; dimBmin..dimBmax)
	foreach (i; 0..duplication)
		"%s,%.15f".writefln(step.anneal!
			(bimsnrtwafom,
			(() => randomPointSet!U(Precision(32), DimensionR(dimR), DimensionF2(dimB))),
			neighborPointSet!1,
			KirkpatrickAcceptance)(Tinit, Tend).second!lg().expand);
}

void usage()
{
	"stat_anneal dimR dimB dimB step from to [duplication]".writeln();
}
