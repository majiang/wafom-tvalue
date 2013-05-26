module pointset;

debug import std.stdio;

import std.exception : enforce;
import std.traits : isUnsigned;
import std.algorithm : min;

import std.random : uniform;
public import sobol : defaultSobols;
import graycode;
import std.conv : to;

/** Test whether the type has method bisect.

A point-set type T is Bisectable if it has bisect method and bisectable property.
*/
template Bisectable(T)
{
    enum Bisectable = __traits(hasMember, T, "bisect") && __traits(hasMember, T, "bisectable");
}

/** Digital Shifted Net over F_2.

Bugs:
If 64 <= dimensionF2, position and length being ulong, foreach (x; P) gives wrong output for p.  Always bisect in such cases, though it's rare in practice.
Examples:
-----
// create non-shifted digital net
auto P = ShiftedBasisPoints(randomVectors!uint(precision, dimensionR, dimensionF2), precision);
assert (P.dimensionF2 == dimensionF2);
assert (P.dimensionR == dimensionR);
assert (P.precision == precision);
foreach (j; P.popFront) assert (j == 0);
foreach (x; P) assert (x.length == dimensionR);
-----
*/
struct ShiftedBasisPoints(T) if (isUnsigned!T)
{
    immutable size_t dimensionF2;
    immutable size_t dimensionR;
    immutable size_t precision;
    immutable ulong length;

    private ulong position;
    private const T[][] basis; /// baisis[i][j] = (i-th vector of basis)'s __j-th__ component
    private const T[] shifter;
    private T[] current;

    this (in T[][] basis, in size_t precision, in T[] shifter)
    in
    {
        enforce(precision <= (T.sizeof << 3));
    }
    body
    {
        this.dimensionF2 = basis.length;
        this.dimensionR = shifter.length;
        assert (precision);
        this.precision = precision;
        this.length = 1UL << this.dimensionF2;

        this.position = 0;
        foreach (b; basis)
            enforce(this.dimensionR == b.length);
        this.basis = basis;
        this.shifter = shifter;
        this.current = this.shifter.dup;
    }
    this (in T[][] basis, in size_t precision)
    {
        enforce(basis.length);
        this (basis, precision, new T[basis[0].length]);
    }

    /// Input range primitives.
    @property T[] front() const
    {
        return current.dup;
    }
    /// ditto
    @property bool empty() const
    {
        return position == length;
    }
    /// ditto
    void popFront()
    {
        enforce(!empty);
        position += 1;
        if (this.empty)
            return;
        current[] ^= basis[position.bottom_zeros()][];
    }

    private alias ShiftedBasisPoints!T SBP;

    /// bisectability
    @property bool bisectable() const
    {
        return 0 < basis.length;
    }
    private alias SBP[2] PSBP;
    /** Return a two-element array of ShiftedBasisPoints.

    Outputs of elements joined together is equal as multiset to output of this.*/
    PSBP bisect() const
    {
        enforce(bisectable);
        assert (precision);
        auto former = SBP(basis[1..$], precision, shifter);
        return [former, former.shift(basis[0])];
        //assert (former.precision);
        //auto ret = PSBP(former, former.shift(basis[0]));
        //assert (ret.former.precision);
        //return ret;
    }

    /// ShiftedBasisPoints with its outputs digital-shifted.
    SBP shift(in T[] further_shifter) const
    {
        return SBP(basis, precision, shifter.XOR(further_shifter));
    }
    /// ShiftedBasisPoints with its outputs bit-shifted.
    SBP opBinary(string op)(in int amount) //const
    {
        if (amount == 0) return this;
        if (amount < 0)
        {
            static if (op == "<<") return this >> -amount;
            static if (op == ">>") return this << -amount;
        }
        auto new_basis = new T[][dimensionF2];
        foreach (i; 0..dimensionF2)
            new_basis[i] = basis[i].to!(T[]);

        size_t new_precision = precision;
        static if (op == "<<") new_precision += amount;
        static if (op == ">>") new_precision -= amount.min(precision);
        auto new_shifter = shifter.to!(T[]);

        static if (op == "<<")
            enforce(new_precision <= T.sizeof << 3, "overflow: ShiftedBasisPoints <<");
        static if (op == ">>")
        {
            if (new_precision == 0)
            {
                foreach (ref l; new_basis)
                    foreach (ref x; l)
                        x = 0;
                foreach (ref x; new_shifter)
                    x = 0;
                return SBP(new_basis, new_precision, new_shifter);
            }
        }
        foreach (ref l; new_basis)
            foreach (ref x; l)
            {
                static if (op == "<<")
                    x <<= amount;
                static if (op == ">>")
                    x >>= amount;
            }
        foreach (ref x; new_shifter)
        {
            static if (op == "<<")
                x <<= amount;
            static if (op == ">>")
                x >>= amount;
        }
        return SBP(new_basis, new_precision, new_shifter);
    }
    /// utility for bit-shifts.
    SBP changePrecision(in size_t new_precision) //const
    {
        if (precision < new_precision)
            return this << (new_precision - precision);
        if (new_precision < precision)
            return this >> (precision - new_precision);
        return this;
    }
}

/// Return a T with random lower precision bits.
T randomBits(T)(in size_t precision) if (isUnsigned!T)
{
    enforce(precision <= T.sizeof << 3);
    return uniform!("[]", T, T)(T.min, T.max)
        >> ((T.sizeof << 3) - precision);
}

/// Return an array of length dimensionR, each element is precision.randomBits.
T[] randomVector(T)(in size_t precision, in size_t dimensionR) if (isUnsigned!T)
{
    T[] ret;
    foreach (i; 0..dimensionR)
        ret ~= precision.randomBits!T;
    return ret;
}

/// Return an array of length count, each element is precision.randomVector(dimensionR).
T[][] randomVectors(T)(in size_t precision, in size_t dimensionR, in size_t count) if (isUnsigned!T)
{
    T[][] ret;
    foreach (i; 0..count)
        ret ~= precision.randomVector!T(dimensionR);
    return ret;
}

/// Utility for point set generation.
ShiftedBasisPoints!T randomBasisPoints(T) (in size_t precision, in size_t dimensionR, in size_t dimensionF2, Flag!"shift" shift) if (isUnsigned!T)
{
    if (shift)
        return precision.nonshiftedRandomBasisPoints(dimensionR, dimensionF2);
    else
        return precision.shiftedRandomBasisPoints(dimensionR, dimensionF2);
}

/// ditto
ShiftedBasisPoints!T nonshiftedRandomBasisPoints(T) (in size_t precision, in size_t dimensionR, in size_t dimensionF2) if (isUnsigned!T)
{
    return ShiftedBasisPoints!T(precision.randomVectors!T(dimensionR, dimensionF2), precision);
}

/// ditto
ShiftedBasisPoints!T shiftedRandomBasisPoints(T) (in size_t precision, in size_t dimensionR, in size_t dimensionF2) if (isUnsigned!T)
{
    return ShiftedBasisPoints!T (precision.randomVectors!T(dimensionR, dimensionF2), precision, precision.randomVector!T(dimensionR));
}

import wafom : wafom;
unittest
{
    auto P = nonshiftedRandomBasisPoints!ubyte(6, 2, 6);
    "P is a SBP with wafom = ".writeln(P.wafom());
    "P.bisect[0].wafom = ".writeln(P.bisect()[0].wafom());
    "P.bisect[1].wafom = ".writeln(P.bisect()[1].wafom());
    "OK?".writeln();
    readln();
}

/// vector componentwise bitwise xor.
T[] XOR(T)(in T[] x, in T[] y) if (isUnsigned!T)
{
    enforce(x.length == y.length);
    T[] ret = new T[x.length];
    ret[] = x[] ^ y[];
    return ret;
}
/// functions for backward compatibility.
auto randomPoints(T)(in size_t dimensionR, in size_t precision, in size_t dimensionF2)
{
    return nonshiftedRandomBasisPoints!T(precision, dimensionR, dimensionF2);
}
/// ditto
ShiftedBasisPoints!ulong transposedBasisPoints(in ulong[][] basis, in size_t precision)
{
    auto new_basis = new ulong[][basis[0].length];
    foreach (i; 0..new_basis.length)
    {
        new_basis[i].length = basis.length;
        foreach (j; 0..basis.length)
        {
            new_basis[i][j] = basis[j][i];
        }
    }
    return ShiftedBasisPoints!ulong(new_basis, precision);
}


import std.array : split;
import std.typecons : Tuple, Flag;
import std.string : strip;
struct DigitalNet(T)
{
    ShiftedBasisPoints!T ps;
    double wafom;
    ulong t;
}

DigitalNet!ulong lineToBP(string line, size_t precision = size_t.max)
{
    ulong[][] basis;
    double wafom;
    ulong t;
    foreach (i, bufs; line.strip().split(",,"))
    {
        auto buf = bufs.split(",");
        if (i == 0)
        {
            t = buf[0].to!ulong();
            wafom = buf[1].to!double();
            continue;
        }
        basis.length += 1;
        foreach (s; bufs.split(","))
        {
            basis[$-1] ~= s.to!ulong();
        }
    }
    assert (basis.length);
    return DigitalNet!ulong(transposedBasisPoints(basis, precision = size_t.max ? basis.guess_precision() : precision), wafom, t);
}

import std.algorithm : max;
size_t guess_precision(ulong[][] basis)
{
    ulong x = 0;
    foreach (l; basis)
        foreach (c; l)
            x = x.max(c);
    size_t precision;
    while (x)
    {
        precision += 1;
        x >>= 1;
    }
    return precision;
}

debug unittest
{
    auto c = 0;
    foreach (x; "5,0.002124192556608,5.236969948020973,,2600265188,692020818,1829963221,894032275,1090497089,651123054,2898340559,1909687544,843513215,1542217271,39519261,3977641622,,2144888475,2941401343,1387697674,1986117176,3702571292,2647056038,3871827325,2263216594,3008901273,4224148358,3048652205,3799831373,,737302895,1233368001,1654098828,2764743171,239054234,249267380,1039474368,3378013260,2468295934,902812364,993745693,2410603677,,3726908047,3018079636,1719761848,2421945980,8259646,1793582138,3611200899,137680621,2493595579,2004711502,1809926346,2378246536\n".
        lineToBP().ps)
    {
        debug (lineToBP) "%s".writefln(x);
        c += 1;
    }
    "unittest passed with %d elements".writefln(c);
}

debug unittest
{
    auto P = randomPoints!ushort(2, 10, 5);
    auto Q = P.shift(randomVector!ushort(10, 2));
    "P =".writeln();
    int i;
    foreach (x; P)
    {
        i.writeln(" -> ", x);
        i += 1;
    }
    i = 0;
    "\nQ =".writeln();
    foreach (x; Q)
    {
        i.writeln(" -> ", x);
        i += 1;
    }
    "OK?".writeln();
    readln();
}
