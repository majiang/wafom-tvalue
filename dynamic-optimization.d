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
	foreach (i; 0..search_size)
	{
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
