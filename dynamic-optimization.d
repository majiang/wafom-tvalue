module dynamic;
import lib.approximate;
import lib.pointsettype;

import std.typecons : Tuple, tuple;
alias U = uint;
alias DN = Tuple!(real, U[][]);

auto ddup(T)(in T[][] arr)
{
	T[][] ret;
	foreach (row; arr)
		ret ~= row.dup;
	return ret;
}

void main(string[] args)
{
	import std.stdio;
	import std.array;
	import std.algorithm : sort;
	import std.conv : to;
	import std.container : heapify;
	import std.math : lg = log2;

	args.popFront();
	if (args.length < 6)
		return stderr.writeln(
"dynamic-optimization s m n v p q
  s : dimension over [0..1)
  m : dimension over {0, 1}
  n : the number of precision digits
  v : target volatility (stdandard deviation / average)
  p : search population
  q : accept population");
	auto dimR = DimensionR(args[0].to!size_t());
	auto dimB = DimensionF2(args[1].to!size_t());
	auto prec = Precision(args[2].to!size_t());
	immutable target_volatility = args[3].to!real();
	immutable search_size = args[4].to!size_t();
	immutable find_size = args[5].to!size_t();

	immutable c = binarysearch(prec, dimR, dimB, target_volatility);
	stderr.writefln("c = %.4f", c);

	auto H = (new DN[find_size]).heapify(0);
	version (NX)
		auto NXS = nxPointSet!U(prec, dimR, dimB);
	foreach (i; 0..search_size)
	{
		version (NX)
			auto P = NXS.scrambleRandomly()[0];
		else
			auto P = randomPointSet!U(prec, dimR, dimB);
		auto w = P.WAFOM(c);
		if (H.length < find_size)
		{
			H.insert(w.tuple(P.basis.ddup));
			continue;
		}
		if (w < H.front[0])
			H.replaceFront(w.tuple(P.basis.ddup));
	}
	foreach (dn; H.release().sort())
		"%s,%.15f,%.15f".writefln(ShiftedBasisPoints!U(dn[1], prec), c, dn[0].lg());
}

version (NX)
auto nxPointSet(U)(Precision prec, DimensionR dimR, DimensionF2 dimB)
{
	import std.stdio, std.string, std.algorithm;
	try
		foreach (P; File("nx-low-s%02d-all.csv".format(dimR.s)).byLine().map!(fromString!U)())
			if (P.precision == prec.n &&
				P.dimensionR == dimR.s &&
				P.dimensionF2 == dimB.m)
				return P;
	catch
		assert (false, "Place a file `nx-low-s%02d-all.csv` in the same directory with the program.".format(dimR.s));
	assert (false, "The file `nx-low-%s02d-all.csv` must have a line for m = %d".format(dimR.s, dimB.m));
}
