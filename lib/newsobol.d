module lib.sobol;

import lib.pointsettype;

import exception = std.exception;
import std.algorithm : map;
import std.conv : to;
import std.array : split, array, front;
debug import std.stdio;

auto sobolPointSet(U)(in Precision precision, in DimensionR dimensionR, in DimensionF2 dimensionF2)
{
	return ShiftedBasisPoints!U(sobolBasis!U(precision, dimensionR, dimensionF2), precision);
}

auto sobolBasis(U)(in Precision precision, in DimensionR dimensionR, in DimensionF2 dimensionF2)
{
	return
		dimensionR.s.initialDirectionNumbers()
		.map!(x => x.split(",").generateDirectionNumbers!U(dimensionF2.m).shift!U(precision.n)())()
		.array().transpose();
}

version (stand_alone_sobol) void main()
{
	import std.stdio;
	import std.string : strip;
	while (true)
	{
		"?".write();
		auto buf = readln().strip().split().map!(to!size_t)();
		sobolBasis!ulong(Precision(buf[0]), DimensionR(buf[1]), DimensionF2(buf[2])).writeln();
	}
}

private:

string[] initialDirectionNumbers(in size_t dimensionR)
{
	//debug "in/out: initial DN".writeln();
	return import ("sobol.csv").split()[0 .. dimensionR];
}

auto generateDirectionNumbers(U)(string[] buf, size_t length)
{
	//debug "in: generate DN".writeln();
	immutable pp = buf[0].to!ulong();
	auto ret = buf[1..$].map!(to!U)().array();
	auto il = ret.length;
	ret.length = length;
	auto d = pp.degree();
	foreach (i; il..length)
	{
		foreach (j; 0..d)
		{
			if (pp >> j & 1)
				ret[i] ^= ret[i - d + j];
			ret[i] <<= 1;
		}
		ret[i] ^= ret[i - d];
	}
	//debug "out: generate DN".writeln();
	return ret;
}

auto shift(U)(U[] x, size_t n)
{
	exception.enforce(x.length <= n, "Sobol sequence needs precision at least dim[F2].");
	debug "in: shift(%s, %d)".writefln(x, n);
	foreach (i, ref c; x)
	{
		debug "shift x[%d]: %d <<= %d - 1 - %d".writefln(i, c, n, i);
		c <<= n - 1 - i;
		debug readln();
	}
	debug "out: shift".writeln();
	return x;
}

auto transpose(T)(T[][] xss)
{
	//debug "in: transpose".writeln();
	auto ret = xss.front.map!(x => [x]).array();
	foreach (xs; xss[1..$])
		foreach (i, x; xs)
			ret[i] ~= x;
	//debug "out: transpose".writeln();
	return ret;
}

size_t degree(ulong polynomial)
in
{
    assert (polynomial);
}
body
{
	debug "in: degree".writeln();
	if (polynomial == 1)
		return 0;
	return (polynomial >> 1).degree() + 1;
}
