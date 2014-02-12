import std.math : exp;

auto pxe(real x)
{
	return (-x).exp();
}

void main(string[] args)
{
	import std.stdio;
	import std.conv : to;
	import std.string : strip;
	foreach (line; args[1..$])
		"%(%.15f,%)".writefln(line.strip().to!size_t().walsh!mb().bitrev());
		//line.strip().to!size_t().walsh!trivial().bitrev().writeln();
		//line.strip().to!size_t().walsh!pxe().bitrev().writeln();
}

auto trivial(real x)
{
	return x < 0.25 ? 1 : 0;
}

auto mb(real x)
{
	return (1 - x).bm();
}

auto bm(real x)
{
	import std.math;
	return (x.log() * -2).sqrt().exp();
}

auto walsh(alias f)(size_t v)
{
	auto ret = v.discretize!f();
	fwt(ret);
	return ret;
}

private:

auto bitrev(T)(T[] x)
{
	import core.bitop;
	auto n = x.length;
	auto ret = new T[n];
	size_t u;
	while (n)
	{
		n <<= 1;
		u += 1;
	}
	foreach (uint i, ref c; ret)
	{
		c = x[i.bitswap() >> u];
		debug
		{
			import std.stdio : writefln;
			import std.conv : text;
			text("%0", 32-u, "b -> %0", 32-u, "b").writefln(i, i.bitswap() >> u);
		}
	}
	return ret;
}

auto discretize(alias f)(size_t v)
{
	import std.algorithm : map;
	import std.array : array;

	immutable s = 1U << v;
	real[] ret;
	foreach (k; 0..s)
		ret ~= integral!f(k, s);
	return ret;
}

auto integral(alias f)(size_t k, size_t s)
{
	immutable p = 1000;
	real ret = 0;
	foreach (i; 0..p)
		ret += f(((i + 0.5L) / p + k) / s);
	return ret / p;
}

import std.traits : isFloatingPoint;

void fwt(T)(T[] x)
    if (isFloatingPoint!T)
{
    import std.math : sqrt;
    fwt_rec(x);
    auto v = (1.0 / x.length);
    foreach (ref e; x)
        e *= v;
}

/// walsh transform
void fwt_rec(T)(T[] x)
    if (isFloatingPoint!T || isIntegral!T)
{
    auto n = x.length;
    if (x.length == 1)
        return;
    assert (!(n & 1));
    n >>= 1;
    foreach (i; 0..n)
    {
    	immutable j = i + n;
        x[i] += x[j];
        x[j] *= -2;
        x[j] += x[i];
    }
    fwt_rec(x[0..n]);
    fwt_rec(x[n..$]);
}
