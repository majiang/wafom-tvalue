module genz;

import std.stdio;
import std.math;

/// simple exponential function: integral is (e-1)^s.
template exponential(size_t s)
{
    double f(double[] x)
    in
    {
        assert (x.length == s);
    }
    body
    {
        auto e = 0.0;
        foreach (c; x)
            e += c;
        return exp(e);
    }
    enum I = (E - 1) ^^ s;
}

template negexp(size_t s)
{
    double f(double[] x)
    in
    {
        assert (x.length == s);
    }
    body
    {
        auto e = 0.0;
        foreach (c; x)
            e += c;
        return exp(s - e);
    }
    enum I = (E - 1) ^^ s;
}

template oscillatoryC(alias u, alias a)
{
    alias oscillatory!() g;
    double f(double[] x)
    {
        return g.f(x, u, a);
    }
    double memoI;
    @property double I()
    {
        if (memoI.isNaN)
            return memoI = g.I(u, a);
        return memoI;
    }
}

template oscillatory()
{
    double f(double[] x, double u, double[] a)
    {
        auto t = 2 * PI * u;
        foreach (i, c; x)
            t += a[i] * c;
        return t.cos();
    }
    double I(double u, double[] a)
    {
        import std.algorithm : reduce;
        auto s = a.length;
        auto t0 = 2 * PI * u - (s & 3) * PI_2;
        double ret = 0;
        foreach (i; 0..(1u << s))
        {
            auto t = t0;
            auto dir = 1;
            foreach (j, c; a)
                if (i >> j & 1)
                    t += c;
                else
                    dir *= -1;
            ret += t.cos() * dir;
        }
        return ret / a.reduce!((p, q) => p * q);
    }
}

version (none) unittest
{
    alias prodpeak!() pp;
    foreach (a; [[1.0], [0.7, 1.3], [0.5, 1.0, 1.5], [0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7]])
        pp.I([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], a).writeln();
}

template prodpeakC(alias u, alias a)
{
    alias prodpeak!() g;
    double f(double[] x)
    {
        return g.f(x, u, a);
    }
    double memoI;
    @property double I()
    {
        if (memoI.isNaN)
            return memoI = g.I(u, a);
        return memoI;
    }
}

template prodpeak()
{
    double f(double[] x, double[] u, double[] a)
    {
        auto ret = 1.0;
        foreach (i, c; x)
            ret /= 1 / (a[i] * a[i]) + (c - u[i]) * (c - u[i]);
        return ret;
    }
    double I(double[] u, double[] a)
    {
        auto ret = 1.0;
        foreach (i, c; a)
            ret *= c * (atan(c * u[i]) - atan(c * u[i] - c));
        return ret;
    }
}

unittest
{
    alias cornpeak!() cp;
    foreach (a; [
        [1.0, 2.0],
        [1.5, 2.0, 2.5],
        [1.0, 1.5, 2.0, 2.5],
        [0.5, 1.0, 1.5, 2.0, 2.5],
        [0.7, 1.0, 1.3, 1.6, 1.9, 2.2]
    ])
        cp.I(a).writeln();
}

template cornpeakC(alias a)
{
    alias cornpeak!() g;
    double f(double[] x)
    {
        return g.f(x, a);
    }
    double memoI;
    @property double I()
    {
        if (memoI.isNaN)
            return memoI = g.I(a);
        return memoI;
    }
}

template cornpeak()
{
    int oddbit(ulong x)
    {
        x = (x >> 32) ^ (x & 0x00000000FFFFFFFF);
        x = (x >> 16) ^ (x & 0x000000000000FFFF);
        x = (x >> 8) ^ (x & 0x00000000000000FF);
        x = (x >> 4) ^ (x & 0x000000000000000F);
        x = (x >> 2) ^ (x & 0x0000000000000003);
        return ((x >> 1) ^ x)
            & 1; // for range propagation
    }
    double f(double[] x, double[] a)
    {
        auto ret = 1.0;
        foreach (i, c; x)
            ret += c * a[i];
        return ret ^^ -(a.length + 1.0);
    }
    double I(double[] a)
    {
        auto ret = 0.0;
        auto s = a.length;
        foreach (i; 0..(1<<s))
        {
            auto cur_denom = 1.0;
            foreach (j, c; a)
                cur_denom += (i >> j & 1) * c;
            ret += (-1) ^^ (i.oddbit()) / cur_denom;
        }
        foreach (i, c; a)
            ret /= (i + 1) * c;
        return ret;
    }
}

template gaussianC(alias a, alias u)
{
    alias gaussian!() g;
    double f(double[] x)
    {
        return g.f(x, a, u);
    }
    double memoI;
    @property double I()
    {
        if (memoI.isNaN)
            return memoI = g.I(a, u);
        return memoI;
    }
}


template gaussian()
{
    double f(double[] x, double[] a, double[] u)
    {
        auto exponent = 0.0;
        foreach (i, c; x)
            exponent += a[i] * a[i] * (c - u[i]) * (c - u[i]);
        return exp(-exponent);
    }
    double I(double[] a, double[] u)
    {
        import std.mathspecial : phi = normalDistribution;
        auto ret = PI ^^ (a.length * 0.5);
        foreach (i, c; a)
            ret *= (phi((1.0-u[i]) * SQRT2 * c) - phi(-u[i] * SQRT2 * c)) / c;
        return ret;
    }
}

template continuousC(alias a, alias u)
{
    alias continuous!() g;
    double f(double[] x)
    {
        return g.f(x, a, u);
    }
    double memoI;
    @property double I()
    {
        if (memoI.isNaN)
            return memoI = g.I(a, u);
        return memoI;
    }
}

template continuous()
{
    double f(double[] x, double[] a, double[] u)
    {
        auto exponent = 0.0;
        foreach (i, c; x)
            exponent += a[i] * abs(c - u[i]);
        return exp(-exponent);
    }
    double I(double[] a, double[] u)
    {
        auto ret = 1.0;
        foreach (i, c; a)
            ret *= -(expm1(-c * u[i]) +expm1(c * (u[i] - 1))) / c;
        return ret;
    }
}

template discontiC(alias a, alias u)
{
    alias g = disconti!();
    double f(double[] x)
    {
        return g.f(x, a, u);
    }
    double memoI;
    @property double I()
    {
        if (memoI.isNaN)
            return memoI = g.I(a, u);
        return memoI;
    }
}

template disconti()
{
    double f(double[] x, double[] a, double[] u)
    {
        foreach (i, z; x)
            if (z > u[i])
                return 0;
        auto exponent = 0.0;
        foreach (i, c; x)
            exponent += a[i] * c;
        return exp(exponent);
    }
    double I(double[] a, double[] u)
    {
        auto ret = 1.0;
        foreach (i, c; a)
            ret *= expm1(c * u[i]) / c;
        return ret;
    }
}

version (none):
import lib.pointset : ShiftedBasisPoints;
alias ShiftedBasisPoints!uint PS;
auto save(PS P)
{
    return PS(P.basis, P.precision);
}
void main()
{
    import lib.pointset : nonshiftedRandomBasisPoints;
    alias nonshiftedRandomBasisPoints!uint genPS;

    import std.random : uniform;
    import std.stdio;

    alias oscillatory!() f1;
    alias prodpeak!() f2;
    alias cornpeak!() f3;
    alias gaussian!() f4;
    alias continuous!() f5;
    alias disconti!() f6;

    immutable size_t precision = 32;
    immutable size_t dimF2 = 14;
    immutable average = 0x1p-14;

    double[] a, u;
    {
        alias f0 = oscillatoryC!(0.5, [0.7, 1.3]);
        auto P = genPS(precision, 2, dimF2);
        import lib.integral : bintegral;
        "I[f] = %.10f; P[f] = %.10f".writefln(f0.I, P.bintegral!(f0.f));
    }
    foreach (i; 1..11)
    {
        auto P = genPS(precision, i, dimF2);
        a ~= uniform(0.2, 1.8);
        u ~= uniform(0.1, 0.9);
        double I;

        I = 0.0;
        foreach (x; P.save)
            I += f1.f(x.scaling(), u[0], a);
        I *= average;
        "I[f] = %.10f; P[f] = %.10f".writefln(I, f1.I(u[0], a));

        I = 0.0;
        foreach (x; P.save)
            I += f2.f(x.scaling(), u, a);
        I *= average;
        "I[f] = %.10f; P[f] = %.10f".writefln(I, f2.I(u, a));

        I = 0.0;
        foreach (x; P.save)
            I += f3.f(x.scaling(), a);
        I *= average;
        "I[f] = %.10f; P[f] = %.10f".writefln(I, f3.I(a));

        I = 0.0;
        foreach (x; P.save)
            I += f4.f(x.scaling(), a, u);
        I *= average;
        "I[f] = %.10f; P[f] = %.10f".writefln(I, f4.I(a, u));

        I = 0.0;
        foreach (x; P.save)
            I += f5.f(x.scaling(), a, u);
        I *= average;
        "I[f] = %.10f; P[f] = %.10f".writefln(I, f5.I(a, u));

        I = 0.0;
        foreach (x; P.save)
            I += f6.f(x.scaling(), a, u);
        I *= average;
        "I[f] = %.10f; P[f] = %.10f".writefln(I, f6.I(a, u));

        writeln();
    }
}

auto scaling(uint[] xs)
{
    immutable centering = 0x1p-33;
    immutable scaling = 0x1p-32;
    double[] ret;
    foreach (x; xs)
        ret ~= x * scaling + centering;
    return ret;
}

