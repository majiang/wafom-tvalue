module cv.bsearch;
import lib.pointsettype;

void main(string[] args)
{
	import std.stdio;
	import std.conv : to;
	if (args.length < 5)
		return stderr.writeln("bsearch-cv s m n v");
	immutable
		s = DimensionR(args[1].to!size_t()),
		m = DimensionF2(args[2].to!size_t()),
		n = Precision(args[3].to!size_t()),
		v = args[4].to!real();
	assert (4 <= s.s && s.s <= 16, "only 4 <= s <= 16 is supported");
	"%.15f".writefln(bs(s, m, n, v));
}

auto bs(F)(in DimensionR s, in DimensionF2 m, in Precision n, in F v, in F left=-1, in F right=1)
{
	import std.string, std.exception;
	debug
	{
		import std.string, std.stdio;
		stderr.writefln("bs(s,%d,n,v,%f,%f)", m.m, left, right);
	}
	enforce(left <= right, "bs(s, m, n, v, %f > %f)".format(left, right));
	enum EPS = 1e-10;
	immutable
		lv = cv!F(s, m, n, left),
		rv = cv!F(s, m, n, right),
		mid = (left + right) * 0.5,
		mv = cv!F(s, m, n, mid);
	if (right - left < EPS)
		return mid;
	enforce(lv >= rv, "bs: f(%f) = %f < %f = f(%f)".format(left, lv, rv, right));
	if (v < rv)
		try
			return bs(s, m, n, v, right, right+1);
		catch
			return bs(s, DimensionF2(m.m-1), n, v);
	if (lv < v)
		try
			return bs(s, m, n, v, left-1, right);
		catch
			return bs(s, DimensionF2(m.m+1), n, v);
	if (mv < v)
		return bs(s, m, n, v, left, mid);
	if (v < mv)
		return bs(s, m, n, v, mid, right);
	return mid;
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
