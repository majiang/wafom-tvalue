module lib.pointset;

debug (working) import std.stdio;

import std.exception : enforce;
import std.traits : isUnsigned;
import std.algorithm : min;
import std.typecons : Flag;

import std.random : uniform;
public import lib.sobol : defaultSobols;
import lib.graycode;
import std.conv : to, text;

private enum BisectMin = 10;

/** Test whether the type has method bisect.

A point-set type T is Bisectable if it has bisect method and bisectable property.
*/
template Bisectable(T)
{
    enum Bisectable =
        __traits(hasMember, T, "bisect") && __traits(hasMember, T, "bisectable");
}

template isPointSet(T)
{
    import std.range : isInputRange;
    enum isPointSet = isInputRange!T && __traits(hasMember, T, "dimensionF2") && __traits(hasMember, T, "precision");
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
    alias T ComponentType;
    immutable size_t dimensionF2;
    immutable size_t dimensionR;
    immutable size_t precision;
    immutable ulong length;

    private ulong position;
    public const T[][] basis; /// baisis[i][j] = (i-th vector of basis)'s __j-th__ component
    public const T[] shift;
    private T[] current;

    this (in T[][] basis, in size_t precision, in T[] shift)
    in
    {
        enforce(precision <= (T.sizeof << 3));
    }
    body
    {
        this.dimensionF2 = basis.length;
        this.dimensionR = shift.length;
        assert (precision);
        this.precision = precision;
        this.length = 1UL << this.dimensionF2;

        this.position = 0;
        foreach (b; basis)
            enforce(this.dimensionR == b.length);
        this.basis = basis;
        this.shift = shift;
        this.current = this.shift.dup;
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
    void reset()
    {
        position = 0;
        this.current = this.shift.dup;
    }

    private alias ShiftedBasisPoints!T SBP;

    /// bisectability
    @property bool bisectable() const
    {
        if (this.dimensionR == 1)
            return 2 <= basis.length;
        return BisectMin < basis.length;
    }
    private alias SBP[2] PSBP;
    /** Return a two-element array of ShiftedBasisPoints.

    Outputs of elements joined together is equal as multiset to output of this.*/
    PSBP bisect() const
    {
        enforce(bisectable);
        auto former = SBP(basis[1..$], precision, shift);
        return [former, former.shifted(basis[0])];
    }

    /// new ShiftedBasisPoints with its outputs digital-shifted.
    SBP shifted(in T[] shift) const
    in
    {
        assert (shift.length == this.dimensionR);
    }
    body
    {
        return SBP(basis, precision, this.shift.XOR(shift));
    }
    /// ditto
    SBP opBinary(string op)(in T[] shift) if (op == "+")
    {
        return this.shifted(shift);
    }
    /// new ShiftedBasisPoints extended by a vector.
    SBP extended(in T[] vector)
    in
    {
        assert (vector.length == this.dimensionR);
    }
    body
    {
        return SBP(basis ~ vector, precision);
    }
    /// ditto
    SBP extended(in T[][] vectors)
    in
    {
        foreach (vector; vectors)
            assert (vector.length == this.dimensionR);
    }
    body
    {
        return SBP(basis ~ vectors, precision);
    }
    /// ditto
    SBP opBinary(string op)(in T[] vector) if (op == "*")
    {
        return this.extended(vector);
    }
    /// ditto
    SBP opBinary(string op)(in T[][] vectors) if (op == "*")
    {
        return this.extended(vectors);
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
        auto new_shift = shift.to!(T[]);

        static if (op == "<<")
            enforce(new_precision <= T.sizeof << 3, "overflow: ShiftedBasisPoints <<");
        static if (op == ">>")
        {
            if (new_precision == 0)
            {
                foreach (ref l; new_basis)
                    foreach (ref x; l)
                        x = 0;
                foreach (ref x; new_shift)
                    x = 0;
                return SBP(new_basis, new_precision, new_shift);
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
        foreach (ref x; new_shift)
        {
            static if (op == "<<")
                x <<= amount;
            static if (op == ">>")
                x >>= amount;
        }
        return SBP(new_basis, new_precision, new_shift);
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
    string toString()
    {
        string ret = text(precision, " ", dimensionF2, " ", dimensionR);
        foreach (l; basis)
            foreach (x; l)
                ret ~= text(" ", x);
        return ret;
    }
}
unittest ///
{
    static assert (isPointSet!(ShiftedBasisPoints!uint));
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
        return precision.nonshiftedRandomBasisPoints!T(dimensionR, dimensionF2);
    else
        return precision.shiftedRandomBasisPoints!T(dimensionR, dimensionF2);
}

/// ditto
ShiftedBasisPoints!T nonshiftedRandomBasisPoints(T) (in size_t precision, in size_t dimensionR, in size_t dimensionF2) if (isUnsigned!T)
{
    return ShiftedBasisPoints!T(precision.randomVectors!T(dimensionR, dimensionF2), precision);
}

/// ditto
ShiftedBasisPoints!T shiftedRandomBasisPoints(T) (in size_t precision, in size_t dimensionR, in size_t dimensionF2) if (isUnsigned!T)
{
    return ShiftedBasisPoints!T(precision.randomVectors!T(dimensionR, dimensionF2), precision, precision.randomVector!T(dimensionR));
}

unittest
{
    import lib.wafom : biwafom;
    import std.math : approxEqual;
    auto P = nonshiftedRandomBasisPoints!uint(32, 4, 12);
    auto
        x = P.bisect()[0].biwafom(),
        y = P.bisect()[1].biwafom(),
        z = P.biwafom();
    assert ((x + y).approxEqual(z * 2));
    debug (verbose)
    {
        "P is a SBP with wafom = ".writeln(z);
        "P.bisect[0].wafom = ".writeln(x);
        "P.bisect[1].wafom = ".writeln(y);
        "average wafom = ".writeln((x + y) * 0.5);
        "OK.".writeln();
    }
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
ShiftedBasisPoints!T transposedBasisPoints(T)(in T[][] basis, in size_t precision) if (isUnsigned!T)
{
    auto new_basis = new T[][basis[0].length];
    foreach (i; 0..new_basis.length)
    {
        new_basis[i].length = basis.length;
        foreach (j; 0..basis.length)
        {
            new_basis[i][j] = basis[j][i];
        }
    }
    return ShiftedBasisPoints!T(new_basis, precision);
}

size_t guess_precision(T)(T[][] basis) if (isUnsigned!T)
{
    import std.algorithm : max;
    T x = 0;
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


unittest
{
    auto P = randomPoints!ushort(2, 10, 5);
    auto Q = P.shifted(randomVector!ushort(10, 2));
    debug (verbose) "P =".writeln();
    int i;
    foreach (x; P)
    {
        debug (verbose) i.writeln(" -> ", x);
        i += 1;
    }
    i = 0;
    debug (verbose) "\nQ =".writeln();
    foreach (x; Q)
    {
        debug (verbose) i.writeln(" -> ", x);
        i += 1;
    }
    debug (verbose) "OK?".writeln();
    debug (verbose) readln();
}

/** construct ShiftedBasisPoints from string.

In case line is a part of CSV data, ignore the first comma and the following.
The rest is interpreted as space-separated values:
    precision dimensionF2 dimensionR <basis>
where <basis> is lexicographically ordered by dimF2 and dimR.
i.e., the first dimR element is the first vector of the basis.
*/
auto fromString(T)(const(char)[] line) if (isUnsigned!T)
{
    import std.string : strip;
    import std.array : split, front, popFront, findSplitBefore;
    auto buf = line.strip().findSplitBefore(",")[0].split();
    immutable precision = buf.front.to!size_t(); buf.popFront();
    immutable dimensionF2 = buf.front.to!size_t(); buf.popFront();
    immutable dimensionR = buf.front.to!size_t(); buf.popFront();
    auto basis = new T[][dimensionF2];
    foreach (i; 0..dimensionF2)
        foreach (j; 0..dimensionR)
        {
            basis[i] ~= buf.front.to!T();
            buf.popFront();
        }
    return ShiftedBasisPoints!T(basis, precision);
}

///
unittest
{
    assert ("32 16 1 2147486260 1073761228 536900390 268705022 134484065 67133232 33829545 17045993 8658469 4205512 2368453 1340400 789021 161910 72838 54636"
            .fromString!uint().front == [0u]);
}

struct InfoPointSet(PointSetType, InfoType)
{
    PointSetType pointSet;
    InfoType info;
    alias info this;
    string toString()
    {
        return pointSet.toString() ~ "," ~ info.toString();
    }
}

unittest
{
    import std.typecons : Tuple;
    alias InfoPointSet!(ShiftedBasisPoints!uint, Tuple!(double, "wafom", ulong, "tvalue")) IPS;
    IPS().toString().writeln();
    readln();
}
