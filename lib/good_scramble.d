import lib.pointsettype;

import std.algorithm : map, reduce, schwartzSort;
import std.array : array;

auto find_scramble(alias criterion, U)(ShiftedBasisPoints!U P, size_t count)
{
	return randomScramblesFor(P, count)
		.schwartzSort!(a => criterion(P * a))().front;
}

version (stand_alone_scramble) void main()
{
	import std.stdio;
	import std.string : strip;
	import lib.newsobol;
	import lib.wafom;

	while (true)
	{
		auto buf = readln().strip().split().map!(to!size_t)();
		auto P = sobolPointSet!uint(Precision(buf[0]), DimensionR(buf[1]), DimensionF2(buf[2]));
		"initial state: %f".writefln(P.bipmswafom());
		"scrambled state: %f".writefln((P * P.find_scramble!bipmswafom(1000)).bipmswafom());
	}
}