module wafom;

import std.functional : memoize;
import std.traits : isUnsigned, ReturnType, ParameterTypeTuple;
import std.math;
import std.conv : text, to;

debug import std.stdio;
debug = speedup;
import pointset : Bisectable, nonshiftedRandomBasisPoints;

private
size_t mu_star(T)(T x, in size_t precision) if (isUnsigned!T)
{
    if (x == 0)
        return 0;
    size_t ret = precision + 1;
    while (x)
    {
        ret -= 1;
        x >>= 1;
    }
    return ret;
}

private double[] _ff_(size_t exponent)(in size_t precision)
{
    static assert (exponent < 2);
    static if (exponent == 0)
        auto ret = [precision * 0.5 + 1];
    else
        auto ret = [1.5 - 0.5 ^^ (precision + 1)];
    foreach (i; 1..(precision + 1))
    {
        static if (exponent == 0)
            ret ~= i * 0.5;
        else
            ret ~= 1.5 * (1 - 0.5 ^^ i);
    }
    return ret;
}


/** Apply bisect algorithm for function f: R -> double.

Params:
f = A function R -> double. f must satisfy the invariant f(P.bisect()[0]) + f(P.bisect()[1]) == f(P) * 2 whenever P is bisectable.
R = input type of f. R must be Bisectable (i.e., has property bisectable and method bisect).

TODO: generalize to f: R -> T
*/
template Bisect(alias f, R)
{
    static assert (Bisectable!R);
    double Bisect(R P)
    {
        if (!P.bisectable)
            return f(P);
        auto Q = P.bisect();
        return (Bisect!(f, R)(Q[0]) + Bisect!(f, R)(Q[1])) * 0.5;
    }
}


static immutable string scale_and_return = q{
    foreach (i; 0..P.dimensionF2)
        ret *= 0.5;
    return ret - 1;
};

/** Compute NRT WAFOM of point set P. */
double nrt(size_t exponent, R)(R P)
{
    import std.algorithm : reduce;
    auto f =P.precision._ff_!exponent();
    double ret = reduce!((ret, B) => ret + reduce!((a, b) => a * f[b.mu_star(P.precision)])(1.0, B))(0.0, P);
    mixin (scale_and_return);
}

/** Compute wafom of a general quasi-Monte Carlo point set.

* Algorithm:
* Equation (4.2) of wafom-arxiv.<ul>
*   <li>input P = (x[n] : 0 <= n < b<sup>m</sup>) in [0..1)<sup>s</sup>.</li>
*   <li>return sum[x in P](prod[i in s][j in m](1+(-1)<sup>x[i,j]</sup>2<sup>-j-1</sup>)-1) / |P|.</li>
* </ul>

* Params:
* P = an array-of-integral-type range which has attributes precision and dimensionF2. If v is an output vector, v*2<sup>-precision</sup> is an element of [0..2<sup>m</sup>)<sup>s</sup>.

Remarks:
Using double, precision > 54 is no better than precision = 54.
*/
double rapid_dick(size_t exponent, R)(R P)
{
    double ret = 0;
    foreach (B; P)
    {
        double cur = 1;
        foreach (l; B)
            cur *= l.rapid_wafom_factor!exponent(P.precision);
        ret += cur;//reduce!((cur, l) => cur * (l.rapid_wafom_factor!exponent(P.precision)))(1.0, B);
    }
    mixin (scale_and_return);
}

/// ditto
double slow_dick(size_t exponent, R)(R P)
{
    double ret = 0;
    auto f = factors(P.precision, exponent);
    foreach (B; P)
    {
        double cur = 1;
        foreach (l; B)
            foreach (j, c; f)
                cur *= c[(l >> j) & 1];
        ret += cur;
    }
    mixin (scale_and_return);
}

version (verbose) unittest
{
    import pointset : randomPoints;
    "dimF2,algorithm,dick,dick2,nrt,nrt2".writeln();
    foreach (i; 8..13)
    {
        auto P = randomPoints!uint(4, 32, i);
        "%d,%s,%.15e,%.15e,%.15e,%.15e,".writefln
            (i, "orig", P.biwafom(), P.bimswafom(), P.binrtwafom(), P.bimsnrtwafom());
        "%d,%s,%.15e,%.15e,%.15e,%.15e,".writefln
            (i, "fast", P.rapid_dick!1(), P.rapid_dick!2().sqrt(), P.nrt!0(), P.nrt!1().sqrt());
        "%d,%s,%.15e,%.15e,".writefln
            (i, "slow", P.slow_dick!1(), P.slow_dick!2().sqrt());
    }
}

/// WAFOM and other figure of merit computed by Bisect algorithm.
template biwafom(R)
{
    auto biwafom(R P)
    {
        return Bisect!(slow_dick!(1, R))(P);
    }
}
template prwafom(R) /// ditto
{
    auto prwafom(R P)
    {
        return precise_dick!(1, R)(P);
    }
}
template bimswafom(R) /// ditto
{
    auto bimswafom(R P)
    {
        return Bisect!(rapid_dick!(2, R))(P).sqrt();
    }
}
template binrtwafom(R) /// ditto
{
    auto binrtwafom(R P)
    {
        return Bisect!(nrt!(0, R))(P);
    }
}
template bimsnrtwafom(R) /// ditto
{
    auto bimsnrtwafom(R P)
    {
        return Bisect!(nrt!(1, R))(P).sqrt();
    }
}

double rapid_wafom_factor(size_t exponent)(ulong x, ptrdiff_t precision)
{
    debug auto memo = memoize!(get_memo_factor!exponent)(); // just memoize for debug, to compile faster.
    else static memo = get_memo_factor!exponent(); // CTFE for release, to execute faster.
    double ret = 1;
    while (0 < precision)
    {
        ret *= memo[precision - 1][x & 255];
        precision -= 8;
        x >>= 8;
    }
    return ret;
}


double[256][64] get_memo_factor(size_t exponent)()
{
    import std.algorithm : min, max;
    double[256][64] ret;
    auto f = _factors(64, exponent);
    foreach (i; 0..64)
    {
        foreach (j; 0..2 << min(i, 7))
        {
            ret[i][j] = 1;
            foreach (k, c; f[$-(i+1)..$-max(0, i-7)])
            {
                ret[i][j] *= c[(j >> k) & 1];
            }
        }
    }
    return ret;
}

double[2][] _factors(size_t precision, size_t exponent = 1)
{
    assert (precision);
    auto ret = new double[2][precision];
    double recip = 1.0;
    foreach (i; 0..exponent) recip /= 2;
    foreach (i; 0..2)
    {
        ret[$-1][i] = recip; // Mr. Yoshiki suggests recip * recip
        foreach_reverse (j; 1..precision)
        {
            ret[j-1][i] = ret[j][i] * recip;
        }
    }
    foreach (j; 0..precision)
    {
        ret[j][0] += 1;
        ret[j][1] *= -1;
        ret[j][1] += 1;
    }
    return ret;
}
alias memoize!_factors factors;

unittest
{
    debug (verbose) "factors(64) = ".writeln();
    foreach (x; factors(64, 1))
        debug (verbose)
            if (!(x[0] == 1.0 && x[1] == 1.0))
                x.writeln();
}

private auto dick_wep(R)(R P)
{
    Polynomial ret;
    foreach (B; P)
    {
        auto cur = Polynomial();
        foreach (l; B)
            foreach (j; 0..P.precision)
                cur *= (P.precision - j) * ((l >> j & 1) ? -1 : 1);
        ret += cur;
    }
    return ret >> P.dimensionF2;
}

/// WAFOM computed from dick weight enumerator polynomial.
double precise_dick(size_t exponent, R)(R P)
{
    return P.dick_wep().substinv!(1 << exponent)();
}

string dick_weight_enumerator_polynomial_csv(R)(R P)
{
    return P.dick_wep().toCSV();
}

/// string representation for dick weight enumerator polynomial.
string dick_weight_enumerator_polynomial(R)(R P)
{
    return P.dick_wep().toString();
}

import std.bigint;
import std.algorithm : max;
private double toDouble(BigInt x)
{
    double f = 1;
    while (long.max < x)
    {
        f *= 2;
        x >>= 1;
    }
    return x.toLong() * f;
}

private struct Polynomial // WAFOM-specified polynomial. NOT FOR GENERAL USE.
{
    BigInt[] coef = [BigInt(1)];
    invariant()
    {
        assert (1 <= this.coef.length);
        assert (this.coef[0] == 1);
    }
    Polynomial opOpAssign(string op)(Polynomial other) if (op == "+")
    {
        this.coef.length = this.coef.length.max(other.coef.length);
        foreach (i, c; other.coef)
            if (i)
                this.coef[i] += c;
        return this;
    }
    Polynomial opBinary(string op)(size_t rhs)
        if (op == ">>")
    {
        Polynomial tmp;
        tmp.coef = this.coef.dup;
        tmp >>= rhs;
        return tmp;
    }
    Polynomial opOpAssign(string op)(size_t rhs)
        if (op == ">>")
    {
        foreach (ref c; coef[1..$])
        {
            assert (c == (c >> rhs) << rhs);
            c >>= rhs;
        }
        return this;
    }
    Polynomial opOpAssign(string op)(ptrdiff_t rhs)
        if (op == "*")
        in { assert (rhs); } body
    {
        immutable bool negative = rhs < 0;
        immutable size_t position = negative ? -rhs : rhs;
        immutable old_length = this.coef.length;
        this.coef.length += position;
        if (negative)
            foreach_reverse (i; 0..old_length)
                this.coef[i + position] -= this.coef[i];
        else
            foreach_reverse (i; 0..old_length)
                this.coef[i + position] += this.coef[i];
        return this;
    }
    string toString()
    {
        auto ret = "1";
        foreach (i, c; this.coef)
            if (i && c)
                ret ~= " + " ~ c.to!string() ~ "x^" ~ i.to!string();
        return ret;
    }
    string toCSV()
    {
        string ret;
        foreach (c; this.coef)
            ret ~= "," ~ c.to!string();
        return ret;
    }
    double substinv(size_t x)()
    {
        BigInt invsubst = 0;
        foreach (i, c; this.coef)
        {
            if (i == 0) continue;
            invsubst += c;
            invsubst *= x;
        }
        assert (0 < invsubst);
        double ret = invsubst.toDouble();
        foreach (i; 0..this.coef.length)
            ret /= x;
        return ret;
    }
    deprecated double subst(double x, size_t scale)
    {
        double ret = 0;// = x;
        foreach_reverse (i, c; coef)
        {
            if (i == 0) break;
            ret += c.toDouble();
            ret *= x;
        }
        foreach (i; 0..scale)
        {
            ret *= 0.5;
        }
        return ret;
    }
}

Polynomial fromCSV(string line)
{
    import std.array : split;
    BigInt[] coef;
    foreach (w; line.split(","))
        coef ~= BigInt(w);
    return Polynomial(coef);
}

version (none) unittest
{
    Polynomial().toString().writeln();
    auto x = Polynomial(3, false);
    auto y = Polynomial(1, false);
    x.toString().writeln(" = ", x.subst(0.5, 0));
    y.toString().writeln(" = ", y.subst(0.5, 0));
    (x * y).toString().writeln(" = ", (x * y).subst(0.5, 0));
}

version (verbose) unittest
{
    foreach (d; 8..17)
    {
        "P.dimensionF2 = %d".writefln(d);
        foreach (i; 0..100)
        {
            auto P = nonshiftedRandomBasisPoints!uint(32, 4, d);
            "%.15f,%.15f,%.15f".writefln(P.wafom(), P.biwafom(), P.nrtwafom());
        }
    }
}

deprecated double wafom(R)(R P)
{
    assert (P.precision);
    double ret = 0;
    debug (speedup) auto f = factors(P.precision, 1);
    foreach (B; P)
    {
        double cur = 1;
        debug (speedup) double cur_backup = 1;
        foreach (l; B)
        {
            cur *= l.wafom_factor(P.precision);
            debug (speedup) foreach (j, c; f) cur_backup *= c[(l >> j) & 1];
        }
        debug (speedup) auto diff = cur - cur_backup;
        debug (speedup) {assert (diff * diff < 1e-10, text("diff = ", diff));}
        ret += cur;
    }
    foreach (i; 0..P.dimensionF2)
        ret *= 0.5;
    return ret - 1;
}

version (verbose) unittest
{
    import pointset : randomPoints;
    foreach (i; 0..10)
    {
        auto P = randomPoints!uint(4, 32, 10);
        auto w = P.wafom();
        auto m = P.mswafom();
        debug (verbose) "wafom = %.15f; mswfm = %.15f".writefln(w, m);
    }
}

deprecated double wafom_factor(ulong x, ptrdiff_t precision)
{
    debug auto memo = memoize!get_memo(); 
    else static memo = get_memo(); 
    double ret = 1;
    while (0 < precision)
    {
        ret *= memo[precision - 1][x & 255];
        precision -= 8;
        x >>= 8;
    }
    return ret;
}

deprecated double[256][64] get_memo()
{
    import std.algorithm : min, max;
    double[256][64] ret;
    auto f = _factors(64);
    foreach (i; 0..64)
        foreach (j; 0..2 << min(i, 7))
        {
            ret[i][j] = 1;
            foreach (k, c; f[$-(i+1)..$-max(0, i-7)])
                ret[i][j] *= c[(j >> k) & 1];
        }
    return ret;
}
