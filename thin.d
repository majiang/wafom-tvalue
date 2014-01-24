import std.range;
import std.typecons : tuple;
import std.stdio;
import std.algorithm : map, sort;

/** Point set generation with threshold.

*/
auto generation(alias generate, alias criterion, T, Args...)(in T threshold, Args args)
{
	static struct R
	{
		T threshold;
		Args args;
		void popFront(){}
		enum empty = false;
		@property auto front()
		{
			while (true)
			{
				auto P = generate(args);
				auto q = criterion(P);
				if (q < threshold)
					return tuple(q, P);
			}
		}
	}
	return R(threshold, args);
}


/** Generate and select sample to plot in a graph.

*/
auto thin
	(alias generate, alias criterion, T, Args...)
	(size_t sample, size_t population, T threshold, Args args)
{
	auto buf =
		generation!(generate, criterion)(threshold, args).
		takeExactly(population).
		array();
	auto idx = buf.zip(buf.length.iota()).map!(a => tuple(a[0][0], a[1]))().array();
	sort(idx);
	auto p = idx.progress(sample);
	return p.map!(a => buf[a[1]][1])();
}

/** Binary search for decreasing function. */
real binarySearch(alias f)(real a, real b)
{
	auto eps = 1e-10;
	real m = (a + b) / 2;
	if (f(m) == 0 || b - a < eps)
		return m;
	if (0 < f(m))
		return binarySearch!f(m, b);
	return binarySearch!f(a, m);
}

/** select sample so that the contiguous pair has difference greater or equal to gap.
*/
auto largegap(T)(T idx, real gap)
{
	auto buf = idx[0 .. 1];
	foreach (i; idx[1 .. $])
		if (gap <= i[0] - buf[$-1][0])
			buf ~= i;
	return buf;
}

auto progress(T)(T idx, ptrdiff_t length)
{
	real gap = binarySearch!
		(g => cast(ptrdiff_t)(largegap(idx, g).length) - length)
		(0, idx[$-1][0] - idx[0][0]);
	return largegap(idx, gap);
}


import lib.pointset : nonshiftedRandomBasisPoints;
import std.math : lg = log2;
import std.conv : to;

version (prep)
void main(string[] buf)
{
	import std.container : heapify;
	if (buf.length != 8)
		return writeln("Usage:\n    thin precision dimR dimB sample < population sample < population");
	auto uints = buf[1..$].map!(to!size_t)();
	size_t
		_precision = uints[0],
		_dimR = uints[1],
		_dimB = uints[2],
		_sample = uints[5],
		_population = uints[6];
	auto pq = (new double[_sample + 1]).heapify(0);
	import lib.wafom : criterion = bipmswafom;
	import std.math : lg = log2;
	foreach (i; 0.._population)
		pq.conditionalInsert(nonshiftedRandomBasisPoints!uint(_precision, _dimR, _dimB).criterion().lg());
	auto rel = pq.release();
	auto ret = rel.front;
	rel.popFront();
	ret += rel.front;
	ret *= 0.5;
	"thin-gen %s %s %s %s %s %e".writefln(
		buf[1],
		buf[2],
		buf[3],
		buf[4],
		buf[5],
		ret
	);
}
else
void main(string[] buf)
{
	buf.popFront();
	if (buf.length != 6)
		return writeln("Usage:\n    thin precision dimR dimB sample < population threshold");
	//auto buf = readln().strip().split();
	auto uints = buf[0..5].map!(to!size_t)();
	size_t
		_precision = uints[0],
		_dimR = uints[1],
		_dimB = uints[2],
		sample = uints[3],
		population = uints[4];
	auto threshold = buf[5].to!real();
	import lib.wafom : bipmswafom;
	foreach (t; thin!(
		((precision, dimR, dimB) => nonshiftedRandomBasisPoints!uint(precision, dimR, dimB)),
		(P => P.bipmswafom().lg())
		)(sample, population, threshold, _precision, _dimR, _dimB))
		t.writeln();
}
