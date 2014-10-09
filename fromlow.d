import lib.pointsettype : Precision, DimensionR, DimensionF2, randomPointSet, ShiftedBasisPoints;
import lib.wafom;
import lib.scramble : deepcopy, shuffleBasis;

auto best(alias generator, alias criterion)(size_t n)
{
	auto ret = [generator()];
	auto bp = criterion(ret[0]);
	foreach (i; 1..n)
	{
		ret ~= generator();
		auto cs = criterion(ret[1]);
		if (bp < cs)
			ret = ret[0..1];
		else
			ret = ret[1..2];
	}
	return ret[0];
}

import std.algorithm : map;
import std.conv : to;
import std.stdio;
import std.math : lg = log2;
alias U = uint;

void main(string[] args)
{
	if (args.length < 6)
		return "fromlow precision dimR[will-be-doubled] dimB n[to find half] m[to find composition]".writeln();
	auto buf = args[1..$].map!(to!size_t)();
	immutable precision = buf[0], dimR = buf[1], dimB = buf[2], n = buf[3], m = buf[4];
	auto P = best!
	(
		() => randomPointSet!U(Precision(precision), DimensionR(dimR), DimensionF2(dimB)),
		biwafom
	)(n);
	"%s,%.15f".writefln(P, P.biwafom().lg());
	auto basis = P.basis.deepcopy();
	auto Q = best!
	(
		() => ShiftedBasisPoints!U(basis.shuffleBasis(), Precision(precision))
		    & ShiftedBasisPoints!U(basis.shuffleBasis(), Precision(precision)),
		biwafom
	)(m);
	"%s,%.15f".writefln(Q, Q.biwafom().lg());
}
