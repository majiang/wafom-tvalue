module lib.pointsettype;

import exception = std.exception;
import conv = std.conv;
import random = std.random;
import std.traits : isUnsigned;
debug import std.stdio;

/// Safe types for point set property.
struct DimensionF2
{
	immutable size_t m;
}
/// ditto
struct DimensionR
{
	immutable size_t s;
	this (size_t s)
	{
		exception.enforce(s);
		this.s = s;
	}
}
/// ditto
struct Precision
{
	immutable size_t n;
	this (size_t n)
	{
		exception.enforce(n <= 64);
		this.n = n;
	}
}
/// Check bisectability.
enum Bisectable(T) = is (typeof ({
	T P;
	if (P.bisectable)
		T[2] x = P.bisect();
}));
/// Check the point set property.
enum isPointSet(T) = is (typeof ({
	T P;
	DimensionF2 dimensionF2 = P.dimensionF2;
	DimensionR dimensionR = P.dimensionR;
	Precision precision = P.precision;
	foreach (x; P)
	{}
}));


/// Standard point set type.
struct ShiftedBasisPoints(U, Size = size_t)
	if (isUnsigned!U)
{
	alias U ComponentType;
	DimensionF2 dimensionF2;///
	DimensionR dimensionR;///
	Precision precision;///
	immutable Size length;
	private
	{
		Size position;
		const U[][] basis;
		const U[] shift;
		U[] current;
	}

	///
	this (in U[][] basis, Precision precision, in U[] shift)
	in
	{
		assert (precision.n);
		assert (precision.n <= (U.sizeof << 3));
	}
	body
	{
		this.dimensionF2 = DimensionF2(basis.length);
		this.dimensionR = DimensionR(shift.length);
		this.precision = precision;
		this.length = (cast(typeof (this.length)) 1) << (this.dimensionF2.m);
		this.position = 0;
		foreach (b; basis)
			assert (this.dimensionR.s == b.length);
		this.basis = basis;
		this.shift = shift;
		this.current = this.shift.dup;
	}
	///
	this (in U[][] basis, Precision precision)
	{
		exception.enforce(basis.length);
		this (basis, precision, new U[basis[0].length]);
	}
	@property U[] front() const/// Input range primitives.
	{
		return current.dup;
	}
	@property bool empty() const/// ditto
	{
		return position == length;
	}
	void popFront()/// ditto
	{
		exception.enforce(!empty);
		position += 1;
		if (this.empty)
			return;
		current[] ^= basis[position.bottom_zeros()][];
	}
	string toString()/// Serialize into ASCII. The inverse of fromString!U.
	{
		string ret = conv.text(precision.n, " ", dimensionF2.m, " ", dimensionR.s);
		foreach (l; basis)
			foreach (x; l)
				ret ~= conv.text(" ", x);
		foreach (x; shift)
			ret ~= conv.text(" ", x);
		return ret;
	}
	typeof (this) opBinary(string op)(in U[] shift)
		if (op == "+")
	{
		return this.shiftBy(shift);
	}
	typeof (this) opBinary(string op)(in U[] vector)
		if (op == "*")
	{
		return this.extendBy(vector);
	}
	typeof (this) opBinary(string op)(in U[][] vectors)
		if (op == "*")
	{
		return this.extendBy(vectors);
	}
	typeof (this) opBinary(string op)(in int amount)
		if (op == "<<" || op == ">>")
	{
		return this.changedTo(Precision(precision.n + (op == "<<" ? 1 : -1) * amount));
	}
}

/// Check if P is bisectable.
@property bool bisectable(U)(in ShiftedBasisPoints!U P)
{
	if (P.dimensionR.s == 1)
		return 1 < P.basis.length;
	return 10 < P.basis.length;
}

/// Bisect P if possible; otherwise throw exception.
ShiftedBasisPoints!U[2] bisect(U)(in ShiftedBasisPoints!U P)
{
	exception.enforce(P.bisectable);
	auto former = ShiftedBasisPoints!U(P.basis[1..$], P.precision, P.shift);
	return [former, former.shifted(basis[0])];
}

/// Shift P by a vector.
ShiftedBasisPoints!U shiftBy(U)(in ShiftedBasisPoints!U P, in U[] shift)
{
	auto new_shift = P.shift.dup;
	new_shift[] ^= shift[];
	return ShiftedBasisPoints!U(P.basis, P.precision, new_shift);
}

/// Extend P by a vector.
ShiftedBasisPoints!U extendBy(U)(in ShiftedBasisPoints!U P, in U[] vector)
{
	return ShiftedBasisPoints!U(P.basis ~ vector, P.precision, P.shift);
}

/// Extend P by vectors.
ShiftedBasisPoints!U extendBy(U)(in ShiftedBasisPoints!U P, in U[][] vectors)
{
	return ShiftedBasisPoints!U(P.basis ~ vectors, P.precision, P.shift);
}

/// Change precision of P.
ShiftedBasisPoints!U changedTo(U)(in ShiftedBasisPoints!U P, Precision precision)
{
	if (P.precision.n == precision.n)
		return P;
	U[][] basis;
	foreach (base; P.basis)
		basis ~= base.dup;
	U[] shift = P.shift.dup;
	if (P.precision.n < precision.n)
	{
		immutable precisionIncrement = precision.n - P.precision.n;
		foreach (base; basis)
			foreach (ref b; base)
				b <<= precisionIncrement;
		foreach (ref s; shift)
			s <<= precisionIncrement;
	}
	else
	{
		immutable precisionDecrement = P.precision.n - precision.n;
		foreach (base; basis)
			foreach (ref b; base)
				b >>= precisionDecrement;
		foreach (ref s; shift)
			s >>= precisionDecrement;
	}
	return ShiftedBasisPoints!U(basis, precision, shift);
}

/// Change dimensionF2 of P.
ShiftedBasisPoints!U changedTo(U)(in ShiftedBasisPoints!U P, DimensionF2 dimensionF2)
{
	return ShiftedBasisPoints!U(basis[0..dimensionF2.m], P.precision, P.shift);
}

/// ditto
ShiftedBasisPoints!U dimensionF2ShrinkBy(U)(in ShiftedBasisPoints!U P, size_t decrement)
{
	return P.changedTo!(DimensionF2(P.dimensionF2.m - decrement));
}

import std.math;
/// Convert to array of (0..1)<sup>s</sup>.
auto toReals(F, U)(ShiftedBasisPoints!U P)
{
	static struct R
	{
		F factor, shift;
		F[] current;
		ShiftedBasisPoints!U P;
		this (ShiftedBasisPoints!U P)
		{
			this.P = P;
			this.factor = 0.5 ^^ (cast(F)(P.precision.n));
			this.shift = this.factor * 0.5;
			foreach (x; P.front)
				current ~= x * factor + shift;
		}
		@property auto empty() const
		{
			return P.empty;
		}
		@property auto front() const
		{
			return current;
		}
		void popFront()
		{
			P.popFront();
			foreach (i, x; P.front)
				current[i] = x * factor + shift;
		}
	}
	return R(P);
}

private U uniform_number(U)(in Precision precision)
	if (isUnsigned!U)
{
	return random.uniform!("[]", U, U)(U.min, U.max) >> ((U.sizeof << 3) - precision.n);
}

U[] uniform_vector(U)(in Precision precision, in DimensionR dimensionR)
{
	U[] ret;
	foreach (i; 0..dimensionR.s)
		ret ~= precision.uniform_number!U();
	return ret;
}

U[][] uniform_vectors(U)(in Precision precision, in DimensionR dimensionR, in DimensionF2 dimensionF2)
{
	U[][] ret;
	foreach (i; 0..dimensionF2.m)
		ret ~= precision.uniform_vector!U(dimensionR);
	return ret;
}

/// Generate a point set by uniformly and randomly selecting its basis.
auto randomPointSet(U)(in Precision precision, in DimensionR dimensionR, in DimensionF2 dimensionF2)
{
	return ShiftedBasisPoints!U(uniform_vectors!U(precision, dimensionR, dimensionF2), precision);
}

/// Generate a vector to shift P.
auto randomShiftFor(U)(ShiftedBasisPoints!U P)
{
	return uniform_vector!U(P.precision, P.dimensionR);
}

/// Generate vectors to shift P.
auto randomShiftsFor(U)(ShiftedBasisPoints!U P, size_t numShift)
{
	// Wraps by inappropriate use of DimensionF2.
	return uniform_vectors!U(P.precision, P.dimensionR, DimensionF2(numShift));
}

/// Shift P by uniformly and randomly selected vector.
auto shiftRandomly(U)(ShiftedBasisPoints!U P)
{
	return P + randomShiftFor(P);
}

/// Scramble P by uniformly and randomly selected nonsingular lower triangular matrix.
/// Time complexity is O(n<sup>2</sup>sm).
/// Currently one cannot retain the scrambling matrix from the result; to fix.
auto scrambleRandomly(U)(ShiftedBasisPoints!U P)
{
	auto basis = new U[][](P.dimensionF2.m, P.dimensionR.s);
	foreach (i; 0..P.dimensionF2.m)
		basis[i][] = P.basis[i][];
	foreach (i; 0..P.dimensionR.s)
		foreach (j; 1..P.precision.n)
			foreach (k; 0..j)
				if (coin())
					foreach (l; 0..P.dimensionF2.m)
						if (basis[l][i] >> j & 1)
							basis[l][i] ^= (1 << k);
	return ShiftedBasisPoints!U(basis, P.precision);
}

/// Guess precision of a point set from its basis.
auto guessPrecision(U)(U[][] basis)
	if (isUnsigned!U)
{
	T x;
	foreach (base; basis)
		foreach (b; base)
			x |= b;
	size_t n;
	while (x)
	{
		n += 1;
		x >>= 1;
	}
	return Precision(n);
}

/// Deserialize from ASCII. The inverse of toString.
auto fromString(U)(const(char)[] line)
	if (isUnsigned!U)
{
	import std.string, std.array, std.algorithm;
	import std.conv : to;
	auto buf = line.strip().findSplitBefore(",")[0].split();
	immutable n = buf.front.to!size_t(); buf.popFront();
	immutable m = buf.front.to!size_t(); buf.popFront();
	immutable s = buf.front.to!size_t(); buf.popFront();
	auto basis = new U[][](m, s);
	foreach (i; 0..m)
		foreach (j; 0..s)
		{
			basis[i][j] = buf.front.to!U();
			buf.popFront();
		}
	if (buf.empty)
		return ShiftedBasisPoints!U(basis, Precision(n));
	exception.enforce(buf.length == s);
	auto shift = new U[s];
	foreach (j; 0..s)
		shift[j] = buf[j].to!U();
	return ShiftedBasisPoints!U(basis, Precision(n), shift);
}


version (stand_alone)
void main()
{
	import std.array;
	import std.stdio;

	auto P = randomPointSet!(ubyte)
		(Precision(8), DimensionR(1), DimensionF2(4));
	auto Q = P.toString().fromString!ubyte();
	foreach (x; P)
	{
		auto y = Q.front;
		if (!Q.empty)
			Q.popFront();
		"%s %s %s".writefln(x, x == y ? "==" : "!=", y);
	}
	auto R = P.shiftRandomly();
	auto S = P.extendBy(R.shift);
	"%((%(%s %))%| %)".writefln(R.array());writeln();
	"%((%(%s %))%| %)".writefln(S.array());
	auto T = ShiftedBasisPoints!ubyte(P.basis, P.precision).scrambleRandomly();
	foreach (x; ShiftedBasisPoints!ubyte(P.basis, P.precision))
	{
		"%(%08b %) -> %(%08b %)".writefln(x, T.front);
		if (!T.empty)
			T.popFront();
	}
}

private auto bottom_zeros(Size)(Size x)
{
	assert (x);
	size_t ret;
	while ((x & 1) == 0)
	{
		x >>= 1;
		ret += 1;
	}
	return ret;
}

private auto coin()
{
	return random.uniform(0, 2) == 1;
}
