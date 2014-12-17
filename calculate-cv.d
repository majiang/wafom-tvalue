module cv.calculate;

import lib.pointsettype;

void main(string[] args)
{
	import std.stdio;
	import std.conv : to;
	if (args.length < 5)
		return stderr.writeln("calculate dimR dimB precision pdiff");
	immutable
		s = args[1].to!size_t(),
		m = args[2].to!size_t(),
		n = args[3].to!size_t(),
		c = args[4].to!real();
	"%.15e".writefln(cv(DimensionR(s), DimensionF2(m), Precision(n), c));
}

F cv(F)(in DimensionR s, in DimensionF2 m, in Precision n, in F c)
{
	return volatility0!F(m) * volatility1(s, n, c);
}

import std.math;

F volatility0(F)(in DimensionF2 m)
{
	F v = 2;
	return (v ^^ m.m - 1).sqrt();
}

auto volatility1(F)(in DimensionR s, in Precision n, in F c)
{
	import std.math;
	F numerator = 1, denominator = 1;
	auto eps = 2 ^^ (c - 1);
	static assert (is (typeof (eps) == F));
	foreach (j; 0..n.n)
	{
		eps *= 0.5;
		numerator *= 1 + eps * eps;
		denominator *= 1 + eps;
	}
	numerator = numerator ^^ s.s - 1;
	denominator = denominator ^^ s.s - 1;
	return numerator.sqrt() / denominator;
}
