import std.algorithm : map, reduce;

import std.exception : enforce;
import std.traits : isUnsigned, ReturnType;
import std.bigint;

/** given a function f and range of ranges P,
return sum[X in P] prod[x in X] f(x).

Does not work with multiplicand_memoized, which is important for Dick-type WAFOMs.
Use for NRT-type WAFOMs.
*/
auto mapprodsum(alias f, R)(R P)
{
    return P.map!(
        X => X.map!f().reduce!((a, b) => a * b)()
    )().reduce!((a, b) => a + b)();
}

string MixinMPS(string variableName, string multiplicand, string Identity, string additionalArguments)
{
    return "auto " ~ variableName ~ " = P.map!(X => X.map!(x => x.multiplicand_memoized!" ~
        "(" ~ multiplicand ~ ", " ~ Identity ~ ")(" ~ additionalArguments ~ ")"
    ~ ").reduce!((a, b) => a * b)())().reduce!((a, b) => a + b)();";
}

// the name P cannot be changed without modification in MixinMPS.
auto preciseDickYoshikiWafom(R)(R P, ptrdiff_t power)
{
    mixin ("ret".MixinMPS(
        "preciseDickYoshikiWeight", "BigFloat(BigInt(1))", "power, P.precision"));
    ret >>= P.dimensionF2;
    ret -= BigFloat(BigInt(1) << ret.exponent, ret.exponent);
    if (power == -2)
        return ret.sqrt();
    else
        return ret;
}

auto standardDickYoshikiWafom(R)(R P, ptrdiff_t power)
{
import std.math : sqrt;
    mixin ("ret".MixinMPS(
        "standardDickYoshikiWeight", "1.0", "power, P.precision"));
    ret  = ret * (0.5 ^^ P.dimensionF2) - 1;
    if (power == -2)
        return ret.sqrt();
    else
        return ret;
}

version (stub)
auto standardNRTWafom(R)(R P, ptrdiff_t power)
{
    auto ret = P.mapprodsum!(/*
        proper
    */)(/*
        implementation
    */);
}

debug void main()
{
import lib.pointset : nonshiftedRandomBasisPoints;
import lib.wafom;
import std.stdio;

    auto P = nonshiftedRandomBasisPoints!uint(16, 4, 8);
    "running main...".writeln();
    auto
        wo1 = P.bipwafom(),
        wo2 = P.bipmswafom();
    "normal algorithm finished.".writeln();
    auto
        ws1 = P.standardDickYoshikiWafom(-1),
        ws2 = P.standardDickYoshikiWafom(-2);
    "generalized algorithm finished.".writeln();
    auto
        wp1 = P.preciseDickYoshikiWafom(-1),
        wp2 = P.preciseDickYoshikiWafom(-2);
    "%.15e ?= %.15e ?= %s".writefln(wo1, ws1, wp1);
    "%.15e ?= %.15e ?= %s".writefln(wo2, ws2, wp2);
}

/** Dick weight proposed by Yoshiki: 1 + (-1)^bit 2^(power*(position+1)).

for WAFOM set power = -1,
for RMS-WAFOM set power = -2.
*/
auto preciseDickYoshikiWeight(size_t position, bool bit, ptrdiff_t power)
{
    immutable positive = 0 < power;
    immutable apow = positive ? power : -power;
    position += 1;

    BigInt ret = 1;
    ret <<= position * apow;
    ret += bit ? -1 : 1;
    return BigFloat(ret, positive ? 0 : position * apow);
}

/// ditto
auto standardDickYoshikiWeight(size_t position, bool bit, ptrdiff_t power)
{
    return 1.0 + (bit ? -1.0 : 1.0) * (2.0 ^^ (power * (1.0 + position)));
}

enum stride = 8; /// The stride of memoizing.
/** given a weight f: (position, bit, power) => G[x],
a sequence of bits U, integer power and precision,
return prod[position] f(position, bit[position], power)(2)
where position of the bit is counted 1-ORIGIN from the LEFT.

Constraints:
G is an abelian group.
*/
auto multiplicand_memoized(alias weight, alias identity, U)(U x, ptrdiff_t power=1, size_t precision=32)
    if (isUnsigned!U)
{
    enum mask = (1 << stride) - 1;

    auto memo = multiplicand_memoizer!(weight, identity, stride)(power, precision);
    auto ret = identity;
    while (memo.length)
    {
        ret *= memo[0][x & mask];
        memo = memo[1..$];
        x >>= stride;
    }
    return ret;
}

/// helper function for multiplicand_memoized.
auto multiplicand_memoizer(alias weight, alias identity, size_t stride)(ptrdiff_t power=1, size_t precision=32)
{
import std.traits : ReturnType;
import std.typecons : Tuple, tuple;
    static ReturnType!weight[][][Tuple!(typeof (power), typeof (precision))] mm;
    if (auto p = tuple(power, precision) in mm)
        return *p;

    auto q = precision / stride;
    auto r = precision % stride;
    auto ret = new ReturnType!weight[][q + (r ? 1 : 0)];
    foreach (i; 0..q)
        foreach (j; 0..1<<stride)
        {
            ret[i] ~= identity;
            foreach (k; 0..stride)
                ret[i][j] *= weight(precision - (stride * i + k), j >> k & 1, power);
        }/*
    foreach (i; [q])*/
    if (r)
        foreach (j; 0..1<<r)
        {
            ret[q] ~= identity;
            foreach (k; 0..r)
                ret[q][j] *= weight(precision - (stride * q + k), j >> k & 1, power);
        }
    return ret;
}


/// Artitrary precision arithmetic; division is not implemented.
struct BigFloat
{
    BigInt significand;
    ptrdiff_t exponent;

    BigFloat opBinary(string op)(BigFloat other)
        if (op == "+" || op == "-" || op == "*")
    {
        auto ret = this;
        mixin ("ret "~op~"= other;");
        return ret;
    }
    BigFloat opBinary(string op)(size_t rhs)
        if (op == ">>" || op == "<<")
    {
        auto ret = this;
        mixin ("ret "~op~"= rhs;");
        return ret;
    }

    BigFloat opOpAssign(string op)(BigFloat other)
        if (op == "+" || op == "-")
    {
        this.sameGap(other);
        mixin ("this.significand "~op~"= other.significand;");
        return this;
    }
    BigFloat opOpAssign(string op)(BigFloat other)
        if (op == "*")
    {
        this.significand *= other.significand;
        //if (maxexponent == 0 || this.exponent + other.exponent < maxexponent)
        this.exponent += other.exponent;
        //else
        //    this.significand >>= other.exponent;
        return this;
    }
//    BigFloat opOpAssign(string op)(BigFloat other)
//        if (op == "/")
//    {
//        this.significand <<= other.exponent;
//        this.significand /= other.significand;
//        return this;
//    }
    BigFloat opOpAssign(string op)(size_t rhs)
        if (op == ">>")
    {
        this.exponent += rhs;
        return this;
    }
    BigFloat opOpAssign(string op)(size_t rhs)
        if (op == "<<")
    {
        this.exponent -= rhs;
        return this;
    }

    BigFloat opUnary(string op)()
        if (op == "+" || op == "-")
    {
        mixin ("return BigFloat(" ~ op ~ "this.significand, this.exponent);");
    }

import std.format : FormatSpec;
    void toString(scope void delegate(const(char)[]) sink, string formatString) const
    {
    import std.conv : to;
        this.significand.toString(sink, formatString);
        if (this.exponent)
        {
            sink(" \\times 2^{");
            sink((-this.exponent).to!string());
            sink("}");
        }
    }

    BigFloat sqrt()
    {
        if (this.significand < 0)
            return -((-this).sqrt());
        auto sqr = this.significand << this.exponent;
        auto ret = BigInt(1) << ((sqr.uintLength) << 4); // initial value
        while (true)
        {
            auto oret = ret;
            ret = (sqr + oret * oret) / (2 * oret);
            if (ret == oret)
                return BigFloat(ret, this.exponent);
        }
        assert (false);
    }

private:
    void sameGap(BigFloat other)
    {
    import std.exception : enforce;
        enforce(this.exponent == other.exponent);
    }
}

version (verbose) unittest
{
import std.stdio;
    auto
        x = BigFloat(BigInt("314159265358979323846264338327") << 20, 40),
        y = BigFloat(BigInt("271828182845904523536028747135") << 20, 40);
    "%s + %s = %s".writefln(x, y, x + y);
    "%s - %s = %s".writefln(x, y, x - y);
    "%s * %s = %s".writefln(x, y, x * y);
//    "%s / %s = %s".writefln(x, y, x / y);
    readln();
}
