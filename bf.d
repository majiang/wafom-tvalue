debug import std.stdio;
import std.algorithm : map, reduce;
import std.exception : enforce;
import std.traits : isUnsigned, ReturnType;
import std.bigint;

/// the upper bound of minimum Dick-Yoshiki weight.
auto DickYoshikiMaxMinWeight(R)(R P)
{
    immutable
        q = P.dimensionF2 / P.dimensionR,
        r = P.dimensionF2 % P.dimensionR;

    return (r + 1) * (q + 2) + P.dimensionR * q * (q + 3) / 2;
}


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

/// alternative for mapprodsum to inject mixin
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
    if (power == -4)
        return ret.qrrt();
    if (power == -3)
        return ret.cbrt();
    if (power == -2)
        return ret.sqrt();
    return ret;
}

auto lowerOnlyDickYoshikiWEP(R)(R P)
{
    return P.polynomialDickYoshikiWafom(P.DickYoshikiMaxMinWeight());
}

auto polynomialDickYoshikiWafom(R)(R P, size_t maxdeg=size_t.max)
{
import std.algorithm : joiner;
    return P.map!(X => Polynomial!BigInt(maxdeg).reduce!(
    (a, b) => a * b)
    (X.map!(x => x.polynomialDickYoshikiWeight(P.precision))().joiner())).reduce!(
    (a, b) => a + b)() >> P.dimensionF2;
}

string to_string(double from)
{
    if (from == 0)
        return "0";
    if (from < 0)
        return "-" ~ (-from).to_string();
    auto pow = 0;
    while (from < 1)
    {
        pow -= 1;
        from *= 10;
    }
    while (10 <= from)
    {
        pow += 1;
        from /= 10;
    }
    return from.normalized_to_string() ~ "e" ~ (pow < 0 ? "-" : "+") ~ (pow < 0 ? -pow : pow).int_to_string();
}

string int_to_string(int from)
in
{
    assert (0 <= from);
}
body
{
    if (from < 10)
        return [cast(immutable char)(from + '0')];
    return (from / 10).int_to_string() ~ (from % 10 + '0');
}

string normalized_to_string(double from)
in
{
    assert (1 <= from && from < 10);
}
body
{
    string ret;
    foreach (i; 0..20)
    {
        foreach_reverse (j; 0..10)
        {
            if (from < j)
                continue;
            ret ~= cast(char)('0' + j);
            from -= j;
            from *= 10;
            break;
        }
    }
    return ret[0..1] ~ "." ~ ret[1..$];
}

auto extendedDickYoshikiWafom(double lg_scale, R)(R P, ptrdiff_t power)
{
    mixin ("ret".MixinMPS(
        "extendedDickYoshikiWeight!(" ~ lg_scale.to_string() ~ ")", "1.0", "power, P.precision"));
    ret = ret * (0.5 ^^ (P.dimensionF2)) - 1;
    if (power < -1)
        return ret ^^ (-1.0 / power);
    return ret;
}

auto standardDickYoshikiWafom(R)(R P, ptrdiff_t power)
{
    mixin ("ret".MixinMPS(
        "standardDickYoshikiWeight", "1.0", "power, P.precision"));
    ret  = ret * (0.5 ^^ P.dimensionF2) - 1;
    if (power < -1)
        return ret ^^ (-1.0 / power);
    return ret;
}

auto preciseNRTWafom(R)(R P, ptrdif_t power)
{
}

auto standardNRTWafom(R)(R P, in ptrdiff_t power)
{
import std.math : sqrt;

    immutable precision = P.precision;
    auto ret = P.mapprodsum!(x => x.standardNRTMultiplicand(power, precision));
    ret = ret * (0.5 ^^ P.dimensionF2) - 1;
    if (power < -1)
        return ret ^^ (-1.0 / power);
    return ret;
}

import std.datetime : Duration, Clock, SysTime;

string toMyString(SysTime t)
{
    auto buf = t.toISOString();
    return
        buf[0..4] ~ '-' ~
        buf[4..6] ~ '-' ~
        buf[6..8] ~ ' ' ~
        buf[9..11] ~ ':' ~
        buf[11..13] ~ ':' ~
        buf[13..18];
}


auto timeit(alias before, alias f, bParams, T...)(bParams bparams, Duration minTime, size_t minCount, T additionalParams)
{

import std.typecons : tuple;
    auto param = before(bparams.expand);
    auto start = Clock.currTime();
    foreach (i; 0..minCount)
        f(param, additionalParams);
    auto consumed = Clock.currTime() - start;
    if (minTime <= consumed)
        return tuple(consumed, minCount);
    return bparams.timeit!(before, f)(minTime, minCount * 2, additionalParams);
}

auto timeit(alias f, R)(R P, Duration minTime, size_t minCount)
{
import std.typecons : tuple;
    auto start = Clock.currTime();
    foreach (i; 0..minCount)
        f(P);
    auto consumed = Clock.currTime() - start;
    if (minTime <= consumed)
        return tuple(consumed, minCount);
    return P.timeit!f(minTime, minCount * 2);
}


void main()
{
    import ui.input : getDigitalNets;
    import lib.wafom : bipwafom, bipmswafom, biwafom, bimswafom;
    import std.stdio : stderr, stdout, writef, writefln, write, writeln, stdin, File;
    import std.datetime : seconds, minutes;
    import std.conv : to;
    import std.array : split;
    import std.string : strip, format;
    import std.algorithm : map;

    import lib.pointset : nonshiftedRandomBasisPoints;
    alias nonshiftedRandomBasisPoints!uint generate;

    immutable precision = 32;
    foreach (line; stdin.byLine())
    {
        auto input = line.strip().split().map!(to!size_t);
        immutable dimR = input[0];
        immutable dimF2 = input[1];
        immutable count = input[2];
        auto output = File("comparison-s%02d-m%02d-%dps.csv".format(dimR, dimF2, count), "w");
        foreach (i; 0..count)
        {
            auto P = generate(precision, dimR, dimF2);
            output.writefln(
                "%s,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,",
                P,
                P.extendedDickYoshikiWafom!(-1.0)(-1),
                P.extendedDickYoshikiWafom!(-1.0)(-2),
                P.extendedDickYoshikiWafom!(-0.5)(-1),
                P.extendedDickYoshikiWafom!(-0.5)(-2),
                P.extendedDickYoshikiWafom!(0.0)(-1),
                P.extendedDickYoshikiWafom!(0.0)(-2),
                P.extendedDickYoshikiWafom!(0.5)(-1),
                P.extendedDickYoshikiWafom!(0.5)(-2),
                P.extendedDickYoshikiWafom!(1.0)(-1),
                P.extendedDickYoshikiWafom!(1.0)(-2),
                P.extendedDickYoshikiWafom!(2.0)(-1),
                P.extendedDickYoshikiWafom!(2.0)(-2)
            );
        }
    }

    version (none) foreach (P; getDigitalNets!uint())
    {
        "(m,s)=(%d,%d): ".writef(P.dimensionF2, P.dimensionR);
        {
            auto b1 = P.timeit!bipwafom(1.seconds(), 1);
            "W1=%e, ".writef(cast(double) b1[0].total!"msecs" / b1[1]);
        }
        {
            auto b2 = P.timeit!bipmswafom(1.seconds(), 1);
            "W2=%e, ".writef(cast(double) b2[0].total!"msecs" / b2[1]);
        }
        stdout.flush();
        {
            auto bx = P.timeit!lowerOnlyDickYoshikiWEP(1.minutes(), 2);
            "WEP=%e. [ms]".writefln(cast(double) bx[0].total!"msecs" / bx[1]);
        }
        stdout.flush();
    }
}
void transform()
{
import ui.input : getDigitalNets;
import lib.wafom : bipwafom, bipmswafom;
import std.math : lg = log2;
import std.stdio : writefln;

    foreach (P; getDigitalNets!uint())
    {
        auto w1 = P.bipwafom().lg();
        auto w2 = P.bipmswafom().lg();
        auto wep = P.lowerOnlyDickYoshikiWEP();
        "%s,%.20f,%.20f,%s".writefln(P, w1, w2, wep);
    }
}

version (none) void main()
{
import lib.pointset : nonshiftedRandomBasisPoints, ShiftedBasisPoints;
alias ShiftedBasisPoints!uint PS;
alias nonshiftedRandomBasisPoints!uint generate;
import lib.wafom;
import std.stdio;
import std.datetime : Clock, dur;

import std.typecons : tuple;
alias tuple!(size_t, size_t, size_t) pscond;

import std.string : strip;
import std.array : split;
import std.conv : to;

    stderr.writeln("dimR, dimF");
    auto buf = readln().strip().split();
    immutable dimR = buf[0].to!size_t();
    immutable dimF = buf[1].to!size_t();

    auto cond = pscond(32, dimR, dimF);
    auto minTime = 1.dur!"seconds"();
    "W1: ".writeln(cond.timeit!(generate, bipwafom!PS)(minTime, 1));
    "W2: ".writeln(cond.timeit!(generate, bipmswafom!PS)(minTime, 1));
    "G1: ".writeln(cond.timeit!(generate, standardDickYoshikiWafom!PS)(minTime, 1, -1));
    "G2: ".writeln(cond.timeit!(generate, standardDickYoshikiWafom!PS)(minTime, 1, -2));
    "G3: ".writeln(cond.timeit!(generate, standardDickYoshikiWafom!PS)(minTime, 1, -3));
    "G4: ".writeln(cond.timeit!(generate, standardDickYoshikiWafom!PS)(minTime, 1, -4));
    "P1: ".writeln(cond.timeit!(generate, preciseDickYoshikiWafom!PS)(minTime, 1, -1));
    "P2: ".writeln(cond.timeit!(generate, preciseDickYoshikiWafom!PS)(minTime, 1, -2));
    "P3: ".writeln(cond.timeit!(generate, preciseDickYoshikiWafom!PS)(minTime, 1, -3));
    "P4: ".writeln(cond.timeit!(generate, preciseDickYoshikiWafom!PS)(minTime, 1, -4));
    "WP: ".writeln(cond.timeit!(generate, lowerOnlyDickYoshikiWEP!PS)(minTime, 1));
}

size_t mu_star(U)(U x, size_t precision)
{
    size_t ret = x ? precision + 1 : 0;
    while (x)
    {
        ret -= 1;
        x >>= 1;
    }
}

auto preciseNRTMultiplicand(U)(in U x, in ptrdiff_t power, in size_t precision)
{
}

auto standardNRTMultiplicand(U)(in U x, in ptrdiff_t power, in size_t precision)
{
    immutable weight = x.mu_star(precision);
    if (power == -1)
    {
        if (weight == 0)
            return 1 + 0.5 * precision;
        return weight * 0.5;
    }
    immutable b = 2.0 ^^ (power + 1.0);
    if (weight == 0)
        return 1 + b * 0.5 * (1 - b ^^ precision) / (1 - b);
    return 1 + 0.5 * ((b - b ^^ weight) / (1 - b) - b ^^ precision);
}


auto polynomialDickYoshikiWeight(U)(U x, size_t precision)
{
    immutable mask = 1 << (precision - 1);
    ptrdiff_t[] ret;
    foreach (i; 2..precision+2)
    {
        ret ~= (x & mask) ? -i : +i;
        x <<= 1;
    }
    return ret;
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
alias extendedDickYoshikiWeight!0.0 standardDickYoshikiWeight;

auto extendedDickYoshikiWeight(double lg_scale)(size_t position, bool bit, ptrdiff_t power)
{
    return 1.0 + (bit ? -1.0 : 1.0) * (2.0 ^^ (power * (1.0 + position - lg_scale)));
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

private struct Polynomial(T) // WAFOM-specified polynomial. NOT FOR GENERAL USE.
{
    T[] coef = [BigInt(1)];
    size_t max_length = size_t.max;
    this (T[] coef)
    {
        this.coef = coef;
    }
    this (size_t max_length)
    {
        this.max_length = max_length;
    }
    this (this)
    {
        coef = coef.dup;
    }
    invariant()
    {
        assert (1 <= this.coef.length);
        assert (this.coef[0] == 1);
    }
    Polynomial opBinary(string op)(Polynomial other) if (op == "+")
    {
        auto ret = this;
        ret += other;
        return ret;
    }
    Polynomial opOpAssign(string op)(Polynomial other) if (op == "+")
    {
    import std.algorithm : min, max;
        this.coef.length = this.coef.length.max(other.coef.length).min(max_length);
        foreach (i; 1..this.coef.length.min(other.coef.length))
            this.coef[i] += other.coef[i];
        return this;
    }
    Polynomial opBinary(string op)(size_t rhs)
        if (op == ">>")
    {
        Polynomial tmp = this;;
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
    Polynomial opBinary(string op)(ptrdiff_t rhs)
        if (op == "*")
    {
        auto ret = this;
        ret *= rhs;
        return ret;
    }
    Polynomial opOpAssign(string op)(ptrdiff_t rhs) /// product with ptrdiff_t rhs, which represents (1 + (sgn rhs) x ^ |rhs|).
        if (op == "*")
        in { assert (rhs); } body
    {
    import std.algorithm : min;
        immutable bool negative = rhs < 0;
        immutable size_t position = negative ? -rhs : rhs;
        if (max_length < position)
            return this;
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
    import std.conv : to;
        auto ret = "1";
        foreach (i, c; this.coef)
            if (i && c)
                ret ~= " + " ~ c.to!string() ~ "x^" ~ i.to!string();
        if (max_length != size_t.max)
            ret ~= " + ...";
        return ret;
    }
    string toCSV()
    {
    import std.conv : to;
        string ret;
        foreach (c; this.coef)
            ret ~= "," ~ c.to!string();
        return ret;
    }
    /// substitute the inverse of x.
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
    /// substitute x and scale.
    version (none) double subst(double x, size_t scale)
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


/// Arbitrary precision arithmetic; division is not implemented.
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

    double toDouble()
    {
        auto exponent = this.exponent;
        auto significand = this.significand;
        auto sign = significand < 0;
        if (sign)
            significand = -significand;
        while (significand > BigInt(1) << 63)
        {
            significand >>= 1;
            exponent -= 1;
        }
        return (sign ? -1 : +1) * significand.toLong() * (0.5 ^^ exponent);
    }

    BigFloat sqrt()
    {
        if (this.significand < 0)
            return -((-this).sqrt());
        auto sqr = this.significand;
        auto exponent = this.significand.uintLength << 5;
        sqr <<= exponent;
        exponent += this.exponent;
        if (exponent & 1)
        {
            exponent += 1;
            sqr <<= 1;
        }
        auto ret = BigInt(1) << (sqr.uintLength << 4); // initial value
        while (true)
        {
            auto oret = ret;
            ret = (sqr + oret * oret) / (2 * oret);
            if (ret == oret)
                return BigFloat(ret, exponent / 2);
        }
        assert (false);
    }

    BigFloat cbrt()
    {
        if (this.significand < 0)
            return -((-this).cbrt());
        auto cub = this.significand;
        auto exponent = this.significand.uintLength << 6; // triple the precision
        cub <<= exponent;
        exponent += this.exponent;
        while (exponent % 3)
        {
            exponent += 1;
            cub <<= 1;
        }
        auto ret = BigInt(1) << (cub.uintLength * 8 / 3); // initial value
        while (true)
        {
            auto oret = ret;
            ret = (2 * ret + cub / (ret * ret)) / 3;
            if (ret == oret)
                return BigFloat(ret, exponent / 3);
        }
        assert (false);
    }
    BigFloat qrrt()
    {
        return this.sqrt().sqrt();
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


template ulonger(size_t n)
{
    static if (1 < n)
        struct ulonger
{
    alias ulonger!n ul;
    alias ulonger!(n - 1) ui;
    alias ulonger!(n - 2) us;
    enum mbits = 32 << n;
    enum hbits = mbits >> 1;
    enum qbits = mbits >> 2;
    private ui[2] w;
    this (ui[2] w) { this.w = w; }
    this (ui x, ui y) { this.w = [x, y]; }
    this (ui x) { this.w[0] = x; }
    static if (2 < n)
        this (ulong x) { this.w[0] = ui(x); }
    this (ui[4] ww)
    {
        this.w = [ww[0], ww[2] + (ww[3] << qbits)];
        this += ul(ww[1]) << qbits;
    }
    ui[4] half()
    {
        static if (2 < n)
            return [ui(w[0].w[0]), ui(w[0].w[1]), ui(w[1].w[0]), ui(w[1].w[1])];
        else
            return [w[0] & 0xFFFF_FFFFUL, w[0] >> 32, w[1] & 0xFFFF_FFFFUL, w[1] >> 32];
    }
    ul opUnary(string op)()
        if (op == "+" || op == "~")
    {
        mixin ("return ul("~op~"w[0], "~op~"w[1]);");
    }
    ul opUnary(string op)()
        if (op == "-")
    {
        return ++(~this);
    }
    ref ul opUnary(string op)()
        if (op == "++" || op == "--")
    {
        mixin (op ~ "w[0];");
        if (op == "++" && w[0].isZero())
            ++w[1];
        if (op == "--" && (w[0] + 1).isZero())
            --w[1];
        return this;
    }
    ptrdiff_t opCmp(ul other)
    {
        auto ret = this.w[1].opCmp(other.w[1]);
        return ret ? ret : this.w[0].opCmp(other.w[0]);
    }
    ref ul opOpAssign(string op)(ptrdiff_t amount)
        if (op == "<<" || op == ">>")
    {
        static if (op == "<<") enum from = 0, to = 1, opop = ">>";
        static if (op == ">>") enum from = 1, to = 0, opop = "<<";
        if (amount < 0)
            mixin ("this "~opop~"= -amount;");
        if (mbits <= amount)
        {
            this = 0;
            return this;
        }
        if (hbits <= amount)
        {
            this.w[to] = this.w[from];
            this.w[from] = 0;
            amount -= hbits;
        }
        if (amount == 0)
            return this;
        mixin ("this.w[to] "~op~"= amount;");
        mixin ("this.w[to] |= this.w[from]" ~ opop ~ "(hbits - amount);");
        mixin ("this.w[from] "~op~"= amount;");
        return this;
    }
    ref ul opOpAssign(string op)(ul other)
        if (op == "+")
    {
        if (w[0] + other.w[0] < w[0])
            ++w[1];
        w[0] += other.w[0];
        w[1] += other.w[1];
        return this;
    }
    ref ul opOpAssign(string op)(ul other)
        if (op == "-")
    {
        this += -other;
        return this;
    }
    ref ul opAssign(ulong x)
    {
        w[0] = x;
        w[1] = 0;
        return this;
    }
    ref ul opOpAssign(string op)(ul other)
        if (op == "*")
    {
        auto tw = this.half(), ow = other.half();
        this = 0;
        foreach (i, c; tw)
            foreach (j, d; ow)
            {
                auto pos = i + j;
                if (pos < 4){
                    this += ul(c * d) << (pos * qbits);
                }
            }
        return this;
    }
    ref ul opOpAssign(string op)(ul other)
        if (op == "|" || op == "&" || op == "^")
    {
        foreach (i; 0..2)
            mixin ("w[i] "~op~"= other.w[i];");
        return this;
    }
    ul opBinary(string op)(ul other)
        if (op == "+" || op == "-" || op == "*" || op == "|" || op == "&" || op == "^")
    {
        auto ret = this;
        mixin ("ret "~op~"= other;");
        return ret;
    }
    ul opBinary(string op)(ptrdiff_t other)
        if (op == "<<" || op == ">>")
    {
        auto ret = this;
        mixin ("ret "~op~"= other;");
        return ret;
    }
    bool isZero()
    {
        return w[0].isZero() && w[1].isZero();
    }
    ref ul opOpAssign(string op)(ulong rhs)
        if (op == "-")
    {
        auto before = w[0];
        w[0] -= rhs;
        if (before < w[0])
            --w[1];
        return this;
    }
    ref ul opOpAssign(string op)(ulong rhs)
        if (op == "+")
    {
        auto before = w[0];
        w[0] += rhs;
        if (w[0] < before)
            ++w[1];
        return this;
    }
    ul opBinary(string op)(ulong rhs)
        if (op == "+" || op == "-")
    {
        auto ret = this;
        mixin ("ret "~op~"= rhs;");
        return ret;
    }
    ubyte lastDigit()
    {
        ubyte ret = (w[0].lastDigit() + w[1].lastDigit() * 6) % 10;
        return ret;
    }
    auto div10()
    {
        auto h = this.half();
        static if (n == -3){
        writeln(this.toNondecimalString(), " - half -> ");
        write("[");
        foreach (i; 0..4)
        {
            write(h[i].toNondecimalString(), ", ");
        }
        writeln("] ->");}
        foreach_reverse (i; 0..4)
        {
            if (h[i].isZero())
                continue;
            auto r = h[i].lastDigit();
            h[i] = h[i].div10();
            if (i)
            {
                static if (n == 2)
                    h[i - 1] += cast(ulong)r << qbits;
                else
                    h[i - 1] += ui(r) << qbits;
            }
        }
        static if (n == -3){
        write("[");
        foreach (i; 0..4)
        {
            write(h[i].toNondecimalString(), ", ");
        }
        writeln("]");}
        return ul(h);
    }
    string toString()
    {
        if (this.isZero)
            return "0";
        string ret = "";
        auto x = this;
        while (!x.isZero)
        {
            auto d = x.lastDigit();
            ret = cast(char)('0' + d) ~ ret;
            x = (x - d).div10();
        }
        return ret;
    }
    string toNondecimalString()
    {
    import std.conv : text;
        static if (n == 2)
            return text(w);
        else
            return text("[", w[0].toNondecimalString(), ", ", w[1].toNondecimalString(), "]");
    }
}
else
{
    static if (n == 0)
        alias uint ulonger;
    static if (n == 1)
        alias ulong ulonger;
}
}


auto random(size_t n)()
{
    static if (1 < n)
    {
        return ulonger!n(random!(n-1)(), random!(n-1)());
    }
    else
    {
        import std.random : uniform;
        static if (n == 0) return 0.uniform!"[]"(uint.max);
        static if (n == 1) return 0.uniform!"[]"(ulong.max);
    }
}

ptrdiff_t opCmp(ulong x, ulong y)
{
    if (x > y) return +1;
    if (x < y) return -1;
    return 0;
}

string toNondecimalString(ulong x)
{
import std.conv : text;
    return text(x);
}

bool isZero(ulong x)
{
    return x == 0;
}

ubyte lastDigit(ulong x)
{
    return x % 10;
}

ulong div10(ulong x)
{
    return x / 10;
}

debug (byPython) unittest
{
    auto x = random!3();
    auto y = random!3();
    auto z = random!3();
    "print (%s + %s - %s) & ((1L << 256) - 1)".writefln(x, y, x + y);
    "print (%s + %s - %s) & ((1L << 256) - 1)".writefln(y, z, y + z);
    "print (%s + %s - %s) & ((1L << 256) - 1)".writefln(z, x, z + x);
    "print (%s - %s - %s) & ((1L << 256) - 1)".writefln(x, y, x - y);
    "print (%s - %s - %s) & ((1L << 256) - 1)".writefln(y, z, y - z);
    "print (%s - %s - %s) & ((1L << 256) - 1)".writefln(z, x, z - x);
    "print (%s * %s - %s) & ((1L << 256) - 1)".writefln(x, y, x * y);
    "print (%s * %s - %s) & ((1L << 256) - 1)".writefln(y, z, y * z);
    "print (%s * %s - %s) & ((1L << 256) - 1)".writefln(z, x, z * x);
}

version (none)
void main()
{
    auto z = ulonger!2(2, 1); //z.writeln();
    auto zz = ulonger!3(z, ulonger!2(1, 2)); zz.writeln();
    auto zzz = ulonger!4(zz, ulonger!3(ulonger!2(1, 2), z)); //zzz.writeln();
    writeln();
    auto x = ulonger!2(0x4000_8000_3000_6000, 0xFFFF_FFFF_FFFF_FFFF);
    "%s + %s = %s".writefln(x, z, x + z);
    "%s * %s = %s".writefln(x, z, x * z);
    writeln();
    auto y = ulonger!2(0x0000_0000_0001_0001, 0x0000_0000_0000_0000);
    "%s * %s = ".writef(y, y);
    y *= y; y.writeln();
    "%s * %s = ".writef(y, y);
    y *= y; y.writeln();
    "%s * %s = ".writef(y, y);
    y *= y; y.writeln();
}
