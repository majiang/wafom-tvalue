import std.range;
import std.typecons : tuple;
import std.stdio;
import std.algorithm : map, sort;

auto generation(alias generate, alias criterion, T)(in T threshold)
{
	static struct R
	{
		T threshold;
		void popFront(){}
		enum empty = false;
		@property auto front()
		{
			while (true)
			{
				auto P = generate();
				auto q = criterion(P);
				if (q < threshold)
					return tuple(q, P);
			}
		}
	}
	return R(threshold);
}


/** Select sample to plot in a graph.

*/
auto thin
	(alias generate, alias criterion, T)
	(size_t sample, size_t population, T threshold)
{
	auto buf =
		generation!(generate, criterion)(threshold).
		takeExactly(sample).
		array();
	auto idx = buf.zip(buf.length.iota()).map!(a => tuple(a[0][0], a[1]))().array();
	sort(idx);
	auto p = idx.progress(population);
	return p.map!(a => buf[a[1]][0])();
}

/** Binary search for decreasing function */
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

auto largegap(T)(T idx, real gap)
{
	auto buf = idx[0 .. 1];
	idx.popFront();
	foreach (i; idx)
	{
		if (gap <= i[0] - buf[$-1][0])
			buf ~= i;
	}
	return buf;
}

auto progress(T)(T idx, ptrdiff_t length)
{
	real gap = binarySearch!
		(g => cast(ptrdiff_t)(largegap(idx, g).length) - length)
		(0, idx[$-1][0] - idx[0][0]);
	return largegap(idx, gap);
}

void main()
{
	import std.math : lg = log2;
	enum
		precision = 32,
		dimR = 4,
		dimB = 12;
	import lib.wafom : bipmswafom;
	import lib.pointset : nonshiftedRandomBasisPoints;
	foreach (t; thin!(
		(() => nonshiftedRandomBasisPoints!uint(precision, dimR, dimB)),
		(P => P.bipmswafom().lg())
		)(10000, 1000, -14.0))
		t.writeln();
}
