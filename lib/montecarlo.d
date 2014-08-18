import lib.pointsettype;

struct MonteCarloPoints(U, Size = size_t)
	if (isUnsigned!U)
{
	import std.exception, std.random;
	alias U ComponentType;
	immutable size_t dimensionF2;
	immutable size_t dimensionR;
	immutable size_t precision;
	immutable Size length;
	enum bisectable = false;
	typeof (this)[2] bisect()
	{
		assert (false);
	}
	private
	{
		Size position;
		U[] current;
		enum maxPrecision = U.sizeof << 3;
	}
	this (Precision precision, DimensionR dimensionR, DimensionF2 dimensionF2)
	{
		Size length = 1;
		foreach (i; 0..dimensionF2.m)
			length <<= 1;
		this.dimensionF2 = dimensionF2.m;
		this.dimensionR = dimensionR.s;
		this.precision = precision.n;
		this.length = length;
		current.length = this.dimensionR;
		foreach (ref c; current)
			c = uniform!U() >> (maxPrecision - this.precision);
	}
	@property U[] front() const
	{
		return current.dup;
	}
	@property bool empty() const
	{
		return position == length;
	}
	void popFront()
	{
		exception.enforce(!empty);
		position += 1;
		if (this.empty)
			return;
		foreach (ref c; current)
			c = uniform!U() >> (maxPrecision - precision);
	}
}

void main()
{
	import std.stdio;
	import lib.testfunction;
	import lib.pointsettype;
	import lib.integral;
	auto P = MonteCarloPoints!uint(Precision(32), DimensionR(4), DimensionF2(16));
	P.signedIntegrationError(genzFactory!real(1, 4)).writeln();
}
