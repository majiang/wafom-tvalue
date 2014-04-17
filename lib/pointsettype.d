module lib.pointsettype;

import exception = std.exception;
import conv = std.conv;
import random = std.random;
import std.range : ElementType;
import std.traits : isUnsigned;
import std.array : array, empty, front, popFront;
import std.typecons : tuple;
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
	size_t
		dimensionF2 = P.dimensionF2,
		dimensionR = P.dimensionR,
		precision = P.precision;
	foreach (x; P)
	{}
}));

static assert (isPointSet!(ShiftedBasisPoints!ubyte));
static assert (isPointSet!(ShiftedBasisPoints!ushort));
static assert (isPointSet!(ShiftedBasisPoints!uint));
static assert (isPointSet!(ShiftedBasisPoints!ulong));
static assert (Bisectable!(ShiftedBasisPoints!ubyte));
static assert (Bisectable!(ShiftedBasisPoints!ushort));
static assert (Bisectable!(ShiftedBasisPoints!uint));
static assert (Bisectable!(ShiftedBasisPoints!ulong));


/// Standard point set type.
struct ShiftedBasisPoints(U, Size = size_t)
	if (isUnsigned!U)
{
	alias U ComponentType;
	immutable size_t dimensionF2;///
	immutable size_t dimensionR;///
	immutable size_t precision;///
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
		this.dimensionF2 = basis.length;
		this.dimensionR = shift.length;
		this.precision = precision.n;
		this.length = (cast(typeof (this.length)) 1) << (this.dimensionF2);
		this.position = 0;
		foreach (b; basis)
			assert (this.dimensionR == b.length);
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
		string ret = conv.text(precision, " ", dimensionF2, " ", dimensionR);
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
	typeof (this) opBinary(string op, T)(T scramble)
		if (op == "*" && is (ElementType!T == bool))
	{
		return this.scrambleBy(scramble);
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
	@property bool bisectable()
	{
		if (this.dimensionR == 1)
			return 1 < this.basis.length;
		return 10 < this.basis.length;
	}
	typeof (this)[2] bisect()
	{
		exception.enforce(this.bisectable);
		auto former = typeof (this)(this.basis[1..$], Precision(this.precision), this.shift);
		return [former, former + basis[0]];
	}
}
version (none)
/// Check if P is bisectable.
@property bool bisectable(U)(in ShiftedBasisPoints!U P)
{
	if (P.dimensionR == 1)
		return 1 < P.basis.length;
	return 10 < P.basis.length;
}

version (none)
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
	return ShiftedBasisPoints!U(P.basis, Precision(P.precision), new_shift);
}

/// Scramble P by a nonsingular lower triangular matrix.
auto scrambleBy(U, S)(ShiftedBasisPoints!U P, S scramble)
	if (is (ElementType!S == bool))
{
	auto basis = new U[][](P.dimensionF2, P.dimensionR);
	foreach (i; 0..P.dimensionF2)
		basis[i][] = P.basis[i][];
	foreach (i; 0..P.dimensionR)
		foreach (j; 1..P.precision)
			foreach (k; 0..j)
			{
				if (scramble.front)
					foreach (l; 0..P.dimensionF2)
						if (basis[l][i] >> j & 1)
							basis[l][i] ^= (1 << k);
				scramble.popFront();
			}
	return ShiftedBasisPoints!U(basis, Precision(P.precision));
}

/** Scramble P by uniformly and randomly selected nonsingular lower triangular matrix.

complexity:
O(n<sup>2</sup>sm) time and O(n(n+m)s) space.

notes:
Though scramble is returned tupled with P scrambled, recovering P is a complexed operation.
*/
auto scrambleRandomly(U)(ShiftedBasisPoints!U P)
{
	auto scramble = randomScrambleFor(P).array(); // freeze the scramble.
	return tuple(P.scrambleBy(scramble), scramble);
}

/// Extend P by a vector or vectors.
ShiftedBasisPoints!U extendBy(U, T)(in ShiftedBasisPoints!U P, in T vector)
	if (is (T == U[]) || is (T == U[][]))
{
	return ShiftedBasisPoints!U(P.basis ~ vector, P.precision, P.shift);
}

// / Extend P by vectors.
//ShiftedBasisPoints!U extendBy(U)(in ShiftedBasisPoints!U P, in U[][] vectors)
//{
//	return ShiftedBasisPoints!U(P.basis ~ vectors, P.precision, P.shift);
//}

/// Change precision of P.
ShiftedBasisPoints!U changedTo(U)(in ShiftedBasisPoints!U P, Precision precision)
{
	if (P.precision == precision.n)
		return P;
	U[][] basis;
	foreach (base; P.basis)
		basis ~= base.dup;
	U[] shift = P.shift.dup;
	if (P.precision < precision.n)
	{
		immutable precisionIncrement = precision.n - P.precision;
		foreach (base; basis)
			foreach (ref b; base)
				b <<= precisionIncrement;
		foreach (ref s; shift)
			s <<= precisionIncrement;
	}
	else
	{
		immutable precisionDecrement = P.precision - precision.n;
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
	return ShiftedBasisPoints!U(basis[0..dimensionF2], P.precision, P.shift);
}

/// ditto
ShiftedBasisPoints!U dimensionF2ShrinkBy(U)(in ShiftedBasisPoints!U P, size_t decrement)
{
	return P.changedTo!(DimensionF2(P.dimensionF2 - decrement));
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
			this.factor = 0.5 ^^ (cast(F)(P.precision));
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


/// Generate a matrix to scramble P.
auto randomScrambleFor(U)(ShiftedBasisPoints!U P)
{
	return (P.dimensionR * (P.precision * (P.precision - 1) / 2)).coins();
}

/// Generate matrices to scramble P.
auto randomScramblesFor(U)(ShiftedBasisPoints!U P, size_t numScramble)
{
	bool[][] ret;
	foreach (i; 0..numScramble)
		ret ~= (P.dimensionR * (P.precision * (P.precision - 1) / 2)).coins().array();
	return ret;
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

private auto coins(size_t count = 0)
{
	static struct R
	{
	private:
		immutable size_t count;
		size_t position;
		size_t rest;
		size_t content;
		bool _empty;
		auto refresh()
		{
			content = random.uniform!("[]", size_t, size_t)(0, size_t.max);
			rest = size_t.sizeof << 3;
		}
		this (size_t count)
		{
			this.count = count;
			refresh();
		}
	public:
		void popFront()
		{
			exception.enforce(!empty);
			position += 1;
			_empty = position == count;
			if (_empty)
				return;
			rest -= 1;
			content >>= 1;
			if (!rest)
				refresh;
		}
		@property bool empty()
		{
			return _empty;
		}
		@property bool front()
		{
			exception.enforce(!empty);
			return content & 1;
		}
	}
	return R(count);
}

version (stand_alone_pointset):

import std.stdio;

void writePoints(R)(R P)
{
	if (P.dimensionR == 1)
		foreach (x; P)
			conv.text("%0", P.precision, "b ").writef(x[0]);
	else
		foreach (x; P)
			conv.text("%(%0", P.precision, "b %)").writefln(x);
	writeln();
}

void writeReals(R)(R P)
{
	if (P.dimensionR == 1)
		foreach (x; P.toReals())
			"%.5f".writef(x[0]);
	else
		foreach (x; P.toReals())
			"%(%.5f %)".writefln(x);
}

void main()
{
	foreach (m; [1, 2])
	{
		auto P = randomPointSet!ubyte
			(Precision(8), DimensionR(m), DimensionF2(4));
		P.writePoints();
		P.toString().fromString!ubyte().writePoints();
		auto Q = P.shiftRandomly();
		Q.writePoints();
		Q.toString().fromString!ubyte().writePoints();
		auto R = P.scrambleRandomly()[0];
		R.writePoints();
		R.toString().fromString!ubyte().writePoints();
	}
}
