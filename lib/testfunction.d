module lib.testfunction;
debug import std.stdio;

import std.math, std.algorithm, std.range;
alias reduce!((a, b) => a + b) sumation;
alias reduce!((a, b) => a * b) product;

interface Function(F)
{
	abstract F opCall(in F[] x);
}

abstract class TestFunction(F) : Function!F
{
	immutable F I;
	this (F I)
	{
		this.I = I;
	}
	override string toString()
	{
		assert (false);
	}
}

interface Validator(F)
{
	bool validate(in F[] x);
}

class LengthValidator(F) : Validator!F
{
	this (size_t length)
	{
		this.length = length;
	}
	bool validate(in F[] x)
	{
		return x.length == length;
	}
	private immutable size_t length;
}

abstract class TestFunctionWithValidation(F) : TestFunction!F
{
	this (Validator!F validator, F I)
	{
		this.validator = validator;
		super (I);
	}
	F opCall(in F[] x)
	in
	{
		assert (validator.validate(x));
	}
	body
	{
		return this.opCallImpl(x);
	}
	abstract F opCallImpl(in F[] x);
protected:
	Validator!F validator;
}

abstract class TestFunctionWithDimensionRValidation(F) : TestFunctionWithValidation!F
{
	this (size_t dimensionR, F I)
	{
		super (new LengthValidator!F(dimensionR), I);
	}
}

final class ExponentialWithInheritance(F) : TestFunctionWithDimensionRValidation!F
{
	this (F a, size_t dimensionR)
	{
		this.a = a;
		super (dimensionR, (expm1(a) / a) ^^ dimensionR);
	}
	override F opCallImpl(in F[] x)
	{
		return (x.sumation() * a).exp();
	}
private:
	F a;
}

private auto genzToString(F)(size_t index, F[] a)
{
	import std.string;
	return "%d %d %(%.15f %)".format(index, a.length, a);
}
private auto appendGenzParameter(F)(string x, F[] u ...)
{
	import std.string;
	return "%s %(%.15f %)".format(x, u);
}

final class Oscillatory(F) : TestFunction!F
{
	this (F[] a, F u)
	{
		this.u = u;
		this.a = a;
		immutable s = a.length;
		immutable F t0 = 2 * PI * u - (s & 3) * PI_2;
		F mI = 0;
		foreach (i; 0..(1u << s))
		{
			F t = t0;
			auto dir = 1;
			foreach (j, c; a)
				if (i >> j & 1)
					t += c;
				else
					dir *= -1;
			mI += t.cos() * dir;
		}
		super (mI / a.product());
	}
	F opCall(in F[] x)
	in
	{
		assert (x.length == a.length);
	}
	body
	{
		F t = 2 * PI * u;
		foreach (i, c; x)
			t += a[i] * c;
		return t.cos();
	}
	override string toString()
	{
		return 0.genzToString(a).appendGenzParameter(u);
	}
private:
	F u;
	F[] a;
}


final class ProductPeak(F) : TestFunction!F
{
	this (F[] a, F[] u)
	in
	{
		assert (u.length == a.length);
	}
	body
	{
		this.a = a;
		this.u = u;
		F mI = 1;
		foreach (i, c; a)
			mI *= c * (atan(c * u[i]) - atan(c * u[i] - c));
		super (mI);
	}
	F opCall(in F[] x)
	in
	{
		assert (x.length == a.length);
	}
	body
	{
		F ret = 1;
        foreach (i, c; x)
            ret /= 1 / (a[i] * a[i]) + (c - u[i]) * (c - u[i]);
        return ret;
	}
	override string toString()
	{
		return 1.genzToString(a).appendGenzParameter(u);
	}
private:
	F[] u, a;
}


final class CornerPeak(F) : TestFunction!F
{
	this (F[] a)
	{
		this.a = a;
		F mI = 0;
		auto s = a.length;
        foreach (i; 0..(1<<s))
        {
            F y = 1;
            F sign = 1;
            foreach (j, c; a)
            {
                y += (i >> j & 1) * c;
                sign *= (i >> j & 1) ? -1 : 1;
            }
            mI += sign / y;
        }
        foreach (i, c; a)
        	mI /= (i + 1) * c;
        super (mI);
 	}
	F opCall(in F[] x)
	in
	{
		assert (x.length == a.length);
	}
	body
	{
		F ret = 1;
		foreach (i, c; x)
			ret += c * a[i];
		return ret ^^ -(a.length + 1.0);
	}
	override string toString()
	{
		return 2.genzToString(a).appendGenzParameter!F();
	}
private:
	F[] a;
}


final class Gaussian(F) : TestFunction!F
{
	this (F[] a, F[] u)
	in
	{
		assert (u.length == a.length);
	}
	body
	{
		this.a = a;
		this.u = u;
        import std.mathspecial : phi = normalDistribution;
		F mI = PI.sqrt() ^^ a.length;
		foreach (i, c; a)
			mI *= (phi((1 - u[i]) * SQRT2 * c) - phi(-u[i] * SQRT2 * c)) / c;
		super (mI);
	}
	F opCall(in F[] x)
	in
	{
		assert (x.length == a.length);
	}
	body
	{
		F exponent = 0;
		foreach (i, c; x)
			exponent += a[i] * a[i] * (c - u[i]) * (c - u[i]);
		return (-exponent).exp();
	}
	override string toString()
	{
		return 3.genzToString(a).appendGenzParameter(u);
	}
private:
	F[] a, u;
}


final class Continuous(F) : TestFunction!F
{
	this (F[] a, F[] u)
	in
	{
		assert (u.length == a.length);
	}
	body
	{
		this.a = a;
		this.u = u;
		F mI = 1;
		foreach (i, c; a)
			mI *= -(expm1(-c * u[i]) + expm1(c * (u[i] - 1))) / c;
		super (mI);
	}
	F opCall(in F[] x)
	in
	{
		assert (x.length == a.length);
	}
	body
	{
		F exponent = 0;
		foreach (i, c; x)
			exponent -= a[i] * abs(c - u[i]);
		return exponent.exp();
	}
	override string toString()
	{
		return 4.genzToString(a).appendGenzParameter(u);
	}
private:
	F[] a, u;
}


final class Discontinuous(F) : TestFunction!F
{
	this (F[] a, F[] u)
	in
	{
		assert (u.length == a.length);
	}
	body
	{
		this.a = a;
		this.u = u;
		F mI = 1;
		foreach (i, c; a)
			mI *= expm1(c * u[i]) / c;
		super (mI);
	}
	F opCall(in F[] x)
	in
	{
		assert (x.length == a.length);
	}
	body
	{
		foreach (i, c; x)
			if (c > u[i])
				return 0;
		F exponent = 0;
		foreach (i, c; x)
			exponent += a[i] * c;
		return exponent.exp();
	}
	override string toString()
	{
		return 5.genzToString(a).appendGenzParameter(u);
	}
private:
	F[] a, u;
}


TestFunction!F stringToGenz(F)(string x)
{
	import std.string : strip;
	import std.conv : to;
	auto buf = x.strip().split();
	immutable size_t
		index = buf[0].to!size_t(),
		dimensionR = buf[1].to!size_t();
	auto a = buf[2..$].map!(to!F)().array();
	// index dimR a... u...
	if (index == 0)
		return new Oscillatory!F(a[0..dimensionR], a[dimensionR]);
	if (index == 2)
		return new CornerPeak!F(a);
	auto u = a[dimensionR..$];
	a = a[0..dimensionR];
	if (index == 1)
		return new ProductPeak!F(a, u);
	if (index == 3)
		return new Gaussian!F(a, u);
	if (index == 4)
		return new Continuous!F(a, u);
	if (index == 5)
		return new Discontinuous!F(a, u);
	assert (false);
}


TestFunction!F genzFactory(F)(size_t index, size_t dimensionR, F difficulty = 1)
{
	import std.random;
	auto rg = Random();
	version (dynamic)
		rg.seed(unpredictableSeed);
	else
	{
		enum constantSeed = 20140807;
		rg.seed(constantSeed);
	}

	auto a = index.genzParameter!F(dimensionR);
	foreach (ref ai; a)
		ai *= difficulty;

	if (index == 0)
		return new Oscillatory!F(a, uniform(0.0, 1.0, rg));
	if (index == 2)
		return new CornerPeak!F(a);
	F[] u;
	foreach (i; 0..dimensionR)
		u ~= uniform(0.0, 1.0);
	if (index == 1)
		return new ProductPeak!F(a, u);
	if (index == 3)
		return new Gaussian!F(a, u);
	if (index == 4)
		return new Continuous!F(a, u);
	foreach (i; 2..dimensionR)
		u[i] = 1;
	if (index == 5)
		return new Discontinuous!F(a, u);
	assert (false);
}

private template genzParameter(F)
{
	F[] _gp = [0.9, 0.725, 0.185, 0.703, 0.204, 0.43];
	auto genzParameter(size_t index, size_t dimensionR)
	out (result)
	{
		import std.stdio;
		assert (result.reduce!((a, b) => a + b)().approxEqual(_gp[index] * dimensionR));
	}
	body
	{
		auto x = _gp[index] * 2 / 3;
		immutable F d = x / (dimensionR - 1);
		F[] a;
		foreach (i; 0..dimensionR)
		{
			a ~= x;
			x += d;
		}
		return a;
	}
}


private F choose(F)(size_t n, size_t r)
{
	if (r == 0)
		return 1;
	return n.choose!F(r - 1) * (n - r + 1) / r;
}

F getIntegralPoweredSumWhichWeWantToBePrivate(F)(in size_t u, in size_t s)
in
{
	assert (s);
}
body
{
	import std.range, std.conv, std.functional;
	alias f = memoize!getIntegralPoweredSumWhichWeWantToBePrivate;
	if (s == 1)
		return 1 / (u + 1).to!F();
	return (u + 1).iota().map!(i => f(u - i, s - 1) * u.choose!F(i) / (i + 1))().sumation();
}

final class PoweredSum(F) : TestFunction!F
{
	this (size_t exponent, size_t dimensionR)
	{
		this.exponent = exponent;
		this.dimensionR = dimensionR;
		super (getIntegralPoweredSumWhichWeWantToBePrivate!F(exponent, dimensionR));
	}
	F opCall(in F[] x)
	in
	{
		assert (x.length == dimensionR);
	}
	body
	{
		return x.sumation() ^^ exponent;
	}
private:
	F exponent;
	size_t dimensionR;
}


final class Exponential(F) : TestFunction!F
{
	this (F a, size_t dimensionR)
	{
		this.a = a;
		this.dimensionR = dimensionR;
		F mI = 0;
		F ct = 1;
		auto i = 1;
		while (mI != mI + ct)
		{
			mI += ct;
			ct *= a;
			i += 1;
			ct /= i;
		}
		super (mI ^^ dimensionR);
	}
	F opCall(in F[] x)
	in
	{
		assert (x.length == dimensionR);
	}
	body
	{
		return (x.sumation() * a).exp();
	}
private:
	F a;
	size_t dimensionR;
}


final class ProductOfCosine(F) : TestFunction!F
{
	this (F a, size_t dimensionR)
	{
		this.a = a;
		this.dimensionR = dimensionR;
		super ((a.sin() / a) ^^ dimensionR);
	}
	F opCall(in F[] x)
	in
	{
		assert (x.length == dimensionR);
	}
	body
	{
		return x.map!(x => (x * a).cos())().product();
	}
private:
	F a;
	size_t dimensionR;
}


final class CosineOfSum(F) : TestFunction!F
{
	this (F a, size_t dimensionR)
	{
		this.a = a;
		this.dimensionR = dimensionR;
		super ((dimensionR + 1).iota().map!(i => ((dimensionR & 3) * -PI_2 + a * i).cos() * (-1.0) ^^ (dimensionR - i) * choose!F(dimensionR, i))().sumation() / a ^^ dimensionR);
	}
	F opCall(in F[] x)
	in
	{
		assert (x.length == dimensionR);
	}
	body
	{
		return (x.sumation() * a).cos();
	}
private:
	F a;
	size_t dimensionR;
}


//final class CornerGaussian(F) : TestFunction!F{}


final class Triangle(F) : TestFunction!F
{
	this (size_t dimensionR)
	{
		this.dimensionR = dimensionR;
		F mI = 0.5;
		super (mI ^^ dimensionR);
	}
	F opCall(in F[] x)
	in
	{
		assert (x.length == dimensionR);
	}
	body
	{
		return x.map!(x => (3*x).min((3*x-2).abs()))().product();
	}
private:
	size_t dimensionR;
}


final class Rectangle(F) : TestFunction!F
{
	this (size_t dimensionR)
	{
		this.dimensionR = dimensionR;
		F mI = 3;
		super (1 / mI ^^ dimensionR);
	}
	F opCall(in F[] x)
	in
	{
		assert (x.length == dimensionR);
	}
	body
	{
        return x.map!((F x) => (((1 <= 3 * x) && (3 * x < 2)) ? -1.0 : 1.0))().product();
	}
private:
	size_t dimensionR;
}


private F[] InverseNormalDistribution(F)(in F[] x)
{
	auto ret = new F[x.length];
	version (box_muller)
	{
		debug exception.enforce(!(x.length & 1));
		immutable n = x.length >> 1;
		foreach (i; 0..n)
		{
			immutable
				r = sqrt(-2 * log(1 - x[i])),
				t = expi(2 * PI * x[i+n]);
			ret[i] = r * t.re;
			ret[i + n] = r * t.im;
		}
	}
	else
	{
		import std.mathspecial;
		foreach (i, c; x)
			ret[i] = normalDistributionInverse(c);
	}
	return ret;
}

class AsianOption(F) : Function!F
{
	this (F r, F T, int s, F P0, F sigma, F k)
	{
		immutable t = T / s;
		this.stdv = sigma * sqrt(t);
		this.A = exp((r - sigma * sigma / 2) * t);
		this.p0ds = P0 / s;
		this.k = k;
		this.reduction = exp(-r * T);
	}
	F opCall(in F[] x)
	{
		const ys = InverseNormalDistribution!F(x);
		F v = 0;
		foreach (y; ys)
		{
			v += 1;
			v *= A * exp(stdv * y);
		}
		v *= p0ds;
		v -= k;
		if (v > 0)
			return v * reduction;
		return 0;
	}
private:
	immutable F A, stdv, p0ds, k, reduction;
}

/// QMC integrand to integrate exp(-A(x, x)) from -infty to b.
class MultivariateNormalProbability(F) : Function!F
{
	// see http://goo.gl/Vf3w2j
	import std.mathspecial : normalDistribution, normalDistributionInverse;

	/// function [0..1]^s -> R
	F opCall(in F[] xraw)
	in
	{
		assert (xraw.length == s);
		foreach (x; xraw)
			assert (0 <= x && x <= 1);
	}
	body
	{
		F[] x;
		e(xraw, x); // x[i] = E[i](xraw[i] * F[i](b[i])).
		F exponent = 0;
		foreach (i, xi; x)
			foreach (j, xj; x)
			{
				if (i == j)
				{
					exponent += xi * A[i][i] * xi;
					break;
				}
				exponent += xi * 2 * A[i][j] * xj;
			}
		return exp(-exponent) * prodfb;
	}
	private void e(in F[] x, out F[] y)
	{
		y.length = this.s;
		foreach (i, xi; x)
			y[i] = ey[i] * normalDistributionInverse(ex[i] * xi);
	}
	this (F[][] A, F[] b, F dd = 1)
	in
	{
		//assert (0 <= dd && dd <= 1);
	}
	body
	{
		F[] diag;
		foreach (i, row; A)
			diag ~= row[i];
		foreach (i, row; A)
			foreach (j, a; row[0..i])
				{
					immutable t = a.abs * dd;
					diag[i] -= t;
					diag[j] -= t;
				}
		this (A, diag, b);
	}
private:
	this (F[][] A, F[] diag, F[] b)
	in
	{
		immutable n = diag.length;
		assert (b.length == n);
		assert (A.length == n);
		foreach (i, row; A)
			assert (row.length == i + 1);
		foreach (c; diag)
			assert (0 < c);
	}
	body
	{
		debug
		{
			stderr.writeln("A, C, b:");
			stderr.writefln("%(%(%.3f %)\n%)\n", A);
			stderr.writefln("%(%.3f %)", diag);
			stderr.writefln("%(%.3f %)", b);
		}
		this.s = diag.length;
		this.A.length = this.s;
		foreach (i, ref row; this.A)
		{
			row = A[i].dup;
			row[i] -= diag[i];
		}
		foreach (i, c; diag)
		{
			ex ~= normalDistribution(sqrt(2 * c) * b[i]);
			ey ~= 1 / sqrt(2 * c);
		}
		debug
		{
			stderr.writefln("F(b), sx, sy:\n%(%(%.3f %)\n%)", [ex, ey]);
		}
		F fb = 1;
		foreach (i, c; diag)
			fb *= sqrt(PI / c) * ex[i];
		this.prodfb = fb;
	}
	immutable size_t s;
	F[][] A;
	immutable(F)[] ex, ey;
public:
	immutable F prodfb;
	// E(F(b)x) = ey * Phi^{-1}(ex*x)
}

version (stand_alone) void main()
{
	import lib.pointsettype, std.stdio, lib.integral;
	enum dimR = 4;
	version (asian_option)
		Function!real[] functions;
	else
		TestFunction!real[] functions;
	version (genz)
		foreach (i; 0..6)
			functions ~= genzFactory!real(i, dimR);
	else version (asian_option)
	{
		functions ~= new AsianOption!real(0.0295588022, 1.0 / 3.0, dimR, 40, 0.2, 35);
	}
	else
	{
		functions ~= new PoweredSum!real(6, dimR);
		"PS".writeln();
		functions ~= new Exponential!real(-2, dimR);
		"Ex".writeln();
		functions ~= new ProductOfCosine!real(2, dimR);
		"PC".writeln();
		functions ~= new CosineOfSum!real(2, dimR);
		"CS".writeln();
		functions ~= new Triangle!real(dimR);
		"Tr".writeln();
		functions ~= new Rectangle!real(dimR);
		"Rc".writeln();
		functions ~= new ExponentialWithInheritance!real(-2, 4);
		functions ~= new ExponentialWithInheritance!real(-2, 6);
	}
	foreach (f; functions)
	{
		version (asian_option){}
		else
			f.I.write(",");
		foreach (dimB; [8, 9, 10, 11, 12, 13, 14])
		{
			auto P = randomPointSet!uint(Precision(32), DimensionR(dimR), DimensionF2(dimB));
			version (asian_option)
				P.integral(f).write(",");
			else
				P.signedIntegrationError(f).write(",");
		}
		writeln();
	}
}
