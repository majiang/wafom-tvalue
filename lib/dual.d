module lib.dual;

import lib.pointsettype;
import std.bigint;
import std.algorithm;
import std.stdio;

version (yoshiki_net)
void main()
{
	foreach (m; 1..10)
	{
		auto b = new uint[][] (m, 1);
		uint v = 1;
		foreach (ref x; b)
		{
			x[0] |= v;
			v <<= 1;
		}
		foreach_reverse (ref x; b)
		{
			x[0] |= v;
//			v <<= 1;
		}
		auto P = ShiftedBasisPoints!uint(b, Precision(2 * m));
		auto dP = P.dualNet();
		assert (P.precision == dP.precision);
		assert (P.dimensionR == dP.dimensionR);
		assert (P.dimensionF2 == dP.dimensionF2);
		writefln("P.mu2 = %(%d, %).", P.mu!2().array().sort());
		writefln("dP.mu2 = %(%d, %).", dP.mu!2().array().sort());
	}

}

version (standalone_dual)
void main()
{
	auto P = randomPointSet!ubyte(Precision(4), DimensionR(2), DimensionF2(4));
	stderr.writefln("%(%(%04b %)\n%)", P.basis);
	auto b = P.encodeBasis();
	stderr.writefln("%(%08b %)", b.map!(x => x.toLong()));
	auto Q = b.decodeBasis!ubyte(DimensionR(2), Precision(4));
	stderr.writefln("%(%(%04b %)\n%)", Q.basis);
	auto c = Q.encodeBasis();
	stderr.writefln("%(%08b %)", c.map!(x => x.toLong()));
	assert (b == c);
	auto d = c.toEchelon();
	stderr.writefln("%(%08b %)", d.map!(x => x.toLong()));
	auto R = d.decodeBasis!ubyte(DimensionR(2), Precision(4));
	stderr.writefln("P = %((%(%d,%))%|, %)", P.array().sort()); 
	stderr.writefln("R = %((%(%d,%))%|, %)", R.array().sort()); 
	auto dP = P.dualNet();
	stderr.writefln("%(%(%04b %)\n%)", dP.basis);
	auto ddP = dP.dualNet();
	stderr.writefln("%(%(%04b %)\n%)", ddP.basis);
	stderr.writefln("dP = %((%(%d,%))%|, %)", dP.array().sort()); 
	stderr.writefln("ddP = %((%(%d,%))%|, %)", ddP.array().sort()); 
	stderr.writefln("dP.mu2 = %(%d, %)", dP.mu!2().array().sort());
}

auto mu(size_t alpha, U)(ShiftedBasisPoints!U P)
{
	stderr.writefln("This is impure.");
	return P.map!(k => k.map!(k => mu!alpha(k, Precision(P.precision)))().reduce!("a+b")());
}

// mu-alpha weight for dualed net.
auto mu(size_t alpha, U)(U k, Precision n)
	if (!alpha)
{
	stderr.writefln("This is impure.");
	return 0;
}

auto mu(size_t alpha, U)(U k, Precision n)
	if (alpha)
{
	stderr.writefln("This is impure.");
	if (!k)
		return 0;
	immutable U l = k & -k;
	import core.bitop : bsf;
	return n.n-l.bsf() + (k ^ l).mu!(alpha - 1, U)(n);
}

BigInt[] encodeBasis(U)(ShiftedBasisPoints!U P)
{
	return P.basis.map!(x => encodeVector!U(x.dup))().array();
}

ShiftedBasisPoints!U decodeBasis(U)(BigInt[] p, DimensionR s, Precision n)
{
	return ShiftedBasisPoints!U(p.map!(x => x.decodeVector!U(s.s))().array(), n);
}

import std.typecons : Tuple;
alias BitInfo = Tuple!(BigInt, "free", BigInt, "bound");

auto dualNet(U)(ShiftedBasisPoints!U P)
{
	return P.encodeBasis().dual(P.dimensionR * P.precision).decodeBasis!U(DimensionR(P.dimensionR), Precision(P.precision));
}

auto dual(BigInt[] b, size_t l)
{
	auto separated = b.toEchelon.separateBits(l);
	auto freeBits = separated.free;
	BigInt[] d;
	while (freeBits)
	{
		d ~= freeBits & (-freeBits);
		freeBits ^= d[$-1];
	}
	assert (d.length + b.length == l);
	foreach (ref x; d)
	{
		BigInt y;
		foreach (z; b)
			if (x & z)
				y |= z;
		x |= y & separated.bound;
	}
	return d;
}

auto separateBits(BigInt[] b, size_t l)
{
	BigInt freeBits, boundBits;
	BigInt x = 1;
	x <<= l - 1;
	while (x)
	{
		if (b.length && (b[0] & x))
		{
			boundBits |= x;
			b = b[1..$];
		}
		else
			freeBits |= x;
		x >>= 1;
	}
	return BitInfo(freeBits, boundBits);
}

auto toEchelon(BigInt[] b)
{
	immutable m = b.length;
	foreach (i; 0..m)
	{
		sort!"a > b"(b);
		foreach (j; 0..m)
		{
			if (i != j)
				b[j] = b[j].min(b[j] ^ b[i]);
		}
	}
	return b;
}

auto encodeVector(U)(U[] v)
{
	immutable s = v.length;
	BigInt r = 0;
	BigInt x = 1;
	while (v.any())
	{
		foreach (i; 0..s)
		{
			if (v[i] & 1)
				r |= x;
			v[i] >>= 1;
			x <<= 1;
		}
	}
	return r;
}

auto decodeVector(U)(BigInt v, size_t s)
{
	auto ret = new U[s];
	U x = 1;
	while (v)
	{
		foreach (i; 0..s)
		{
			if (v & 1)
				ret[i] |= x;
			v >>= 1;
		}
		x <<= 1;
	}
	return ret;
}
