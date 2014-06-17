module lib.wafom;

import std.functional : memoize;
import std.traits : isUnsigned;
import std.math;
import std.conv : to;

import lib.pointsettype : Bisectable, isPointSet;

debug import std.stdio;
//debug import lib.pointset : nonshiftedRandomBasisPoints;

template biwafom(R) if (Bisectable!R) /// Bisect Dick WAFOM
{
    auto biwafom(R P)
    {
        return Bisect!(rapid_dick!(1, R))(P);
    }
}
template bipwafom(R) if (Bisectable!R)
{
    auto bipwafom(R P)
    {
        return Bisect!(rapid_proper_dick!(1, R))(P);
    }
}
template bipmswafom(R) if (Bisectable!R)
{
    auto bipmswafom(R P)
    {
        return Bisect!(rapid_proper_dick!(2, R))(P).sqrt();
    }
}
template prwafom(R) /// Precise Dick WAFOM
{
    auto prwafom(R P)
    {
        return precise_dick!(1, R)(P);
    }
}
template bimswafom(R) if (Bisectable!R) /// Bisect root mean square Dick WAFOM
{
    auto bimswafom(R P)
    {
        return Bisect!(rapid_dick!(2, R))(P).sqrt();
    }
}
template binrtwafom(R) if (Bisectable!R) /// Bisect NRT WAFOM
{
    auto binrtwafom(R P)
    {
        return Bisect!(nrt!(1, R))(P);
    }
}
template bimsnrtwafom(R) if (Bisectable!R) /// Bisect root mean square NRT WAFOM
{
    auto bimsnrtwafom(R P)
    {
        return Bisect!(nrt!(2, R))(P).sqrt();
    }
}


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

private real[] _ff_(size_t exponent)(in size_t precision)
{
    static assert (exponent == 1 || exponent == 2);
    static if (exponent == 0)
        real[] ret = [precision * 0.5 + 1];
    else
        real[] ret = [1.5 - 0.5 ^^ (precision + 1)];
    foreach (i; 1..(precision + 1))
    {
        static if (exponent == 1)
            ret ~= i * 0.5;
        static if (exponent == 2)
            ret ~= 1.5 * (1 - 0.5 ^^ i);
    }
    return ret;
}


/** Apply bisect algorithm for function f: R -> real.

Params:
f = A function R -> real. f must satisfy the invariant f(P.bisect()[0]) + f(P.bisect()[1]) == f(P) * 2 whenever P is bisectable.
R = input type of f. R must be Bisectable (i.e., has property bisectable and method bisect).

TODO: generalize to f: R -> T
*/
template Bisect(alias f, R) if (Bisectable!R)
{
    real Bisect(R P)
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

/** Compute NRT WAFOM of a digital net. */
real nrt(size_t exponent, R)(R P)
{
    import std.algorithm : reduce;
    auto f =P.precision._ff_!exponent();
    real ret = reduce!((ret, B) => ret + reduce!((a, b) => a * f[b.mu_star(P.precision)])(1.0, B))(0.0, P);
    mixin (scale_and_return);
}

/** Compute Dick WAFOM of a digital net.

* Algorithm:
* Equation (4.2) of wafom-arxiv.<ul>
*   <li>input P = (x[n] : 0 <= n < b<sup>m</sup>) in [0..1)<sup>s</sup>.</li>
*   <li>return sum[x in P](prod[i in s][j in m](1+(-1)<sup>x[i,j]</sup>2<sup>-j-1</sup>)-1) / |P|.</li>
* </ul>

* Params:
* P = an array-of-integral-type range which has attributes precision and dimensionF2. If v is an output vector, v*2<sup>-precision</sup> is an element of [0..2<sup>m</sup>)<sup>s</sup>.

Remarks:
Using real, precision > 54 is no better than precision = 54.
*/
real rapid_dick(size_t exponent, R)(R P)
{
    real ret = 0;
    foreach (B; P)
    {
        real cur = 1;
        foreach (l; B)
            cur *= l.rapid_wafom_factor!exponent(P.precision);
        ret += cur;//reduce!((cur, l) => cur * (l.rapid_wafom_factor!exponent(P.precision)))(1.0, B);
    }
    mixin (scale_and_return);
}
real rapid_proper_dick(size_t exponent, R)(R P)
{
    real ret = 0;
    foreach (B; P)
    {
        real cur = 1;
        foreach (l; B)
            cur *= l.rapid_wafom_factor!(exponent, 1)(P.precision);
        ret += cur;//reduce!((cur, l) => cur * (l.rapid_wafom_factor!exponent(P.precision)))(1.0, B);
    }
    mixin (scale_and_return);
}

/// ditto
real slow_dick(size_t exponent, R)(R P)
{
    real ret = 0;
    auto f = factors(P.precision, exponent);
    foreach (B; P)
    {
        real cur = 1;
        foreach (l; B)
            foreach (j, c; f)
                cur *= c[(l >> j) & 1];
        ret += cur;
    }
    mixin (scale_and_return);
}

version (verbose) unittest
{
    import lib.pointset : randomPoints;
    "dimF2,algorithm,dick,dick2,nrt,nrt2".writeln();
    foreach (i; 8..13)
    {
        auto P = randomPoints!uint(4, 32, i);
        "%d,%s,%.15e,%.15e,%.15e,%.15e,".writefln
            (i, "orig", P.biwafom(), P.bimswafom(), P.binrtwafom(), P.bimsnrtwafom());
        "%d,%s,%.15e,%.15e,%.15e,%.15e,".writefln
            (i, "fast", P.rapid_dick!1(), P.rapid_dick!2().sqrt(), P.nrt!1(), P.nrt!2().sqrt());
        "%d,%s,%.15e,%.15e,".writefln
            (i, "slow", P.slow_dick!1(), P.slow_dick!2().sqrt());
    }
}


real rapid_wafom_factor(size_t exponent, size_t offset = 0)(ulong x, ptrdiff_t precision)
{
    debug auto memo = memoize!(get_memo_factor!exponent)(offset); // just memoize for debug, to compile faster.
    else static memo = get_memo_factor!exponent(offset); // CTFE for release, to execute faster.
    real ret = 1;
    while (0 < precision)
    {
        ret *= memo[precision - 1][x & 255];
        precision -= 8;
        x >>= 8;
    }
    return ret;
}


real[256][64] get_memo_factor(size_t exponent)(size_t offset)
{
    import std.algorithm : min, max;
    real[256][64] ret;
    auto f = _factors(64, exponent, offset);
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

real[2][] _factors(size_t precision, size_t exponent = 1, size_t offset = 0)
{
    assert (precision);
    auto ret = new real[2][precision];
    real recip = 1.0;
    foreach (i; 0..exponent) recip /= 2;
    foreach (i; 0..2)
    {
        if (offset)
            ret[$-1][i] = recip * recip; // Mr. Yoshiki suggests recip * recip
        else
            ret[$-1][i] = recip;
        foreach_reverse (j; 1..precision)
            ret[j-1][i] = ret[j][i] * recip;
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
real precise_dick(size_t exponent, R)(R P)
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
private real toDouble(BigInt x)
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
    size_t max_length = 2147483647;
    this (BigInt[] coef)
    {
        this.coef = coef;
    }
    this (size_t max_length)
    {
        this.max_length = max_length;
    }
    invariant()
    {
        assert (1 <= this.coef.length);
        assert (this.coef[0] == 1);
    }
    Polynomial opOpAssign(string op)(Polynomial other) if (op == "+")
    {
        this.coef.length = this.coef.length.max(other.coef.length).min(max_length);
        foreach (i, c; other.coef)
            if (0 < i && i < max_length)
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
        auto old_length = this.coef.length;
        this.coef.length = (this.coef.length + position).min(max_length);
        old_length = old_length.min(max_length - position);
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
