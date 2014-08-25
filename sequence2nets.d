import ui.input, lib.pointsettype;
import std.conv : to;
import std.algorithm : map;
import std.stdio;

void main(string[] args)
{
	auto P = getDigitalNet!uint();
	if (args.length == 1)
		return "sequence2nets dimB...".writeln();
	foreach (m; args[1..$].map!(to!size_t)())
		P.changedTo(DimensionF2(m)).writeln();
}
