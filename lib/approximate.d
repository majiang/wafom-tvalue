module lib.approximate;
import lib.integral, lib.pointsettype;
import lib.testfunction : Exponential;
import std.functional : memoize;

auto _exponential(F)(in F c, size_t dimR)
{
	auto a = -(2.0 ^^ c);
	return new Exponential!F(a, dimR);
}
alias exponential(F) = memoize!(_exponential!F);

auto WAFOM(F, R)(R P, in F c)
	if (Bisectable!R)
{
	return P.signedIntegrationError(exponential!F(c, P.dimensionR));
}

auto volatility(F)(Precision n, DimensionR s, DimensionF2 m, in F c)
{
	import std.math;
	F numerator = 1, denominator = 1;
	foreach (h; 0..m.m)
	{
		numerator /= 4;
		denominator /= 2;
	}
	numerator = denominator - numerator;
	immutable ret = numerator.sqrt() / denominator;
	numerator = 1;
	denominator = 1;
	F eps = 2.0 ^^ (c - 1);
	foreach (j; 0..n.n)
	{
		eps *= 0.5;
		numerator *= 1 + eps * eps;
		denominator *= 1 + eps;
	}
	numerator = numerator ^^ s.s;
	denominator = denominator ^^ s.s;
	numerator -= 1;
	denominator -= 1;
	if (numerator == 0)
		return F.nan;
	return ret * numerator.sqrt() / denominator;
}

auto binarysearch(F, F EPS = 0.0001)(in Precision n, in DimensionR s, in DimensionF2 m, in F targetVolatility, in F left = 0, in F right = 1)
in
{
	assert (left <= right);
	assert (volatility!F(n, s, m, left) >= volatility!F(n, s, m, right));
}
body
{
	if (right - left < EPS)
		return (left + right) * 0.5;
	immutable lv = volatility!F(n, s, m, left);
	immutable rv = volatility!F(n, s, m, right);
	if (lv < targetVolatility) // move left
		return binarysearch(n, s, m, targetVolatility, left * 3 - right * 2, left);
	if (rv > targetVolatility) // move right
		return binarysearch(n, s, m, targetVolatility, right, right * 3 - left * 2);
	immutable middle = (left + right) * 0.5;
	immutable mv = volatility!F(n, s, m, middle);
	if (mv < targetVolatility)
		return binarysearch(n, s, m, targetVolatility, left, middle);
	if (mv > targetVolatility)
		return binarysearch(n, s, m, targetVolatility, middle, right);
	assert (false);
}


version (standalone_approximate) void main()
{
	import std.stdio;
	auto P = randomPointSet!uint(Precision(32), DimensionR(4), DimensionF2((8)));
	"%(%.4f\n%)".writefln([
/*		volatility(Precision(32), DimensionR(5), DimensionF2(16), 2.0),
		volatility(Precision(32), DimensionR(5), DimensionF2(10), 1.0),
		volatility(Precision(32), DimensionR(9), DimensionF2(12), 0.0),
		volatility(Precision(32), DimensionR(9), DimensionF2(8), -1.0),
		volatility(Precision(32), DimensionR(14), DimensionF2(8), -2.0),//*/
		binarysearch(Precision(32), DimensionR(4), DimensionF2(12), 1.0),
		binarysearch(Precision(32), DimensionR(6), DimensionF2(12), 1.0),
		binarysearch(Precision(32), DimensionR(9), DimensionF2(12), 1.0),
		binarysearch(Precision(32), DimensionR(9), DimensionF2(8), 1.0),
		binarysearch(Precision(32), DimensionR(14), DimensionF2(8), 1.0),//*/
			]);
	//foreach (real c; [-1.0, 0, 1, 2])
	//	P.WAFOM(c).writeln();
}
