module lib.montecarlo;

import lib.pointsettype;

struct UniformPoints(Size = size_t)
{
	import std.random : uniform;
	immutable size_t dimensionF2;
	immutable size_t dimensionR;
	immutable Size _length;
	enum bisectable = false;
	typeof (this)[2] bisect()
	{
		assert (false);
	}

	this (DimensionR dimensionR, DimensionF2 dimensionF2)
	{
		this.dimensionF2 = dimensionF2.m;
		this.dimensionR = dimensionR.s;
		Size length = 1;
		length <<= dimensionF2.m;
		this._length = length;
	}
	auto toReals(F)()
	{
		struct UniformRealsImpl
		{
			Size length;
			F[] current;
			this (Size _length, size_t dimensionR)
			{
				this.length = _length;
				foreach (i; 0..dimensionR)
					this.current ~= 0.uniform!("[)", F, F)(1);
			}
			@property bool empty()
			{
				return length == 0;
			}
			void popFront()
			{
				assert (!empty);
				length -= 1;
				foreach (ref x; current)
					x = 0.uniform!("[)", F, F)(1);
			}
			@property F[] front()
			{
				return current.dup;
			}
		}
		return UniformRealsImpl(_length, dimensionR);
	}
}

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

version (standalone_montecarlo) void main()
{
	import std.stdio;
	import lib.testfunction;
	import lib.pointsettype;
	import lib.integral;
	auto P = MonteCarloPoints!uint(Precision(32), DimensionR(4), DimensionF2(16));
	P.signedIntegrationError(genzFactory!real(1, 4)).writeln();
}
