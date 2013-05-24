module pointset;

debug import std.stdio;

import std.exception : enforce;

import std.random : uniform;
public import sobol : defaultSobols;
import graycode;

/** Generate a point set of dimension, precision and lg(length) specified by choosing its basis randomly.
*/
auto randomPoints(immutable size_t dimension, immutable size_t precision, immutable size_t lg_length)
{
    return BasisPoints(dimension.random_basis(precision, lg_length), precision);
}

debug (none) unittest
{
    foreach (x; randomPoints(3, 10, 4))
    {
        x.writeln();
    }
    "Press Enter to continue:".writeln();
    readln();
}

/** Input Range for point set.

Generate linear combinations of basis using gray code algorithm.
*/
struct BasisPoints
{
    immutable size_t dimension;
    immutable size_t lg_length;
    immutable ulong length;
    immutable size_t precision;
    alias length opDollar;
    private ulong _position;
    /*private /**/ulong[][] basis; /// baisis[i][j] = i-th component of v[j]
    private ulong[] current;
    @property ulong position()
    {
        return _position;
    }
    this (ulong[][] basis, size_t precision) // TODO: basis[i][j] -> basis[j][i] is better?
    {
        this.dimension = basis.length;
        assert (0 < this.dimension);
        this.precision = precision;
        this.lg_length = basis[0].length;
        this.length = 1UL << this.lg_length;
        this._position = 0;
        this.current.length = this.dimension;
        this.basis = basis;
    }
    this (BasisPoints other)
    {
        this.dimension = other.dimension;
        this.length = other.length;
        this.precision = other.precision;
        this.lg_length = other.lg_length;
        this._position = other._position;
        this.basis = other.basis;
        this.current.length = other.current.length;
        this.current[] = other.current[];
    }
    @property ulong[] front()
    {
        return this.current;
    }
    void popFront()
    {
        this._position += 1;
        if (this.empty)
        {
            return;
        }
        auto j = this.position.bottom_zeros();
        foreach (i, c; this.basis)
        {
            this.current[i] ^= c[j];
        }
    }
    @property bool empty()
    {
        return this.length <= this.position;
    }
    @property BasisPoints save()
    {
        return BasisPoints(this);
    }
    BasisPoints truncatePrecision()
    {
        ulong[][] new_basis;
        new_basis.length = this.basis.length;
        foreach (i; 0..new_basis.length)
        {
            new_basis[i].length = this.lg_length;
            foreach (j, c; this.basis[i])
            {
                new_basis[i][j] = c >> (this.precision - this.lg_length);
            }
        }
        return BasisPoints(new_basis, this.lg_length);
    }
}

struct ShiftedBasisPoints(T)
if (is (T == ubyte) || is (T == ushort) || is (T == uint) || is (T == ulong))
{
    immutable size_t dimensionF2;
    immutable size_t dimensionR;
    immutable size_t precision;
    immutable ulong length;
    alias length opDollar;

    private ulong position;
    private immutable (T[])[] basis; /// baisis[i][j] = __j-th__ component of __v[i]__
    private immutable T[] shifter;
    private T[] current;

    this (in immutable (T[])[] basis, in size_t precision, in T[] shifter)
    {
        this.dimensionF2 = basis.length;
        this.dimensionR = shifter.length;
        this.precision = precision;
        this.length = 1UL << this.dimensionF2;

        this.position = 0;
        foreach (b; basis)
            enforce(this.dimensionR == b.length);
        this.basis = basis.idup;
        this.shifter = shifter.idup;
        this.current = this.shifter.dup;
    }
    this (in immutable (T[])[] basis, in size_t precision)
    {
        enforce(basis.length);
        this (basis, precision, new T[basis[0].length]);
    }

    @property T[] front()
    {
        return current.dup;
    }
    @property bool empty()
    {
        return position == length;
    }
    void popFront()
    {
        enforce(!empty);
        position += 1;
        if (this.empty)
            return;
        current[] ^= basis[position.bottom_zeros()][];
    }

    @property bool bisectable()
    {
        return 0 < basis.length;
    }
    Tuple!(ShiftedBasisPoints!T, ShiftedBasisPoints!T) half()
    {
        enforce(bisectable);
        return Tuple!(ShiftedBasisPoints!T, ShiftedBasisPoints!T)(
            ShiftedBasisPoints!T(basis[1..$], precision, shifter),
            ShiftedBasisPoints!T(basis[1..$], precision, shifter.XOR(basis[0]))
        );
    }
}

T[] randomVector(T)(size_t precision, size_t dimensionR)
if (is (T == ubyte) || is (T == ushort) || is (T == uint) || is (T == ulong))
{
    T[] ret;
    foreach (i; 0..dimensionR)
    {
        auto x = uniform!("[]", T, T)(T.min, T.max);
        foreach (j; precision..(x.sizeof << 3))
        {
            x >>= 1;
        }
        ret ~= x;
    }
    return ret;
}

immutable (T[])[] randomVectors(T)(size_t precision, size_t dimensionR, size_t count)
if (is (T == ubyte) || is (T == ushort) || is (T == uint) || is (T == ulong))
{
    immutable (T[])[] ret;
    foreach (i; 0..count)
        ret ~= precision.randomVector!T(dimensionR).idup;
    return ret;
}

unittest
{
    auto P = ShiftedBasisPoints!ubyte(randomVectors!ubyte(6, 2, 6), 6);
    foreach (X; P)
        X.writeln();
    " ... SBP output".writeln();
    readln();
}

auto XOR(T)(in T[] x, in T[] y)
if (is (T == ubyte) || is (T == ushort) || is (T == uint) || is (T == ulong))
{
    enforce(x.length == y.length);
    auto ret = new T[x.length];
    ret[] = x[] ^ y[];
    return ret.idup;
}


auto shift(R)(R P, ulong[] x)
{
    struct Result
    {
        R PS;
        alias PS this;
        @property auto front()
        {
            auto ret = PS.front.dup;
            foreach (i; 0..ret.length)
            {
                ret[i] ^= x[i];
            }
            return ret;
        }
    }
    return Result(P.save);
}

debug unittest
{
    auto P = randomPoints(2, 10, 5);
    auto Q = shift(P, random_vector(2, 10));
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

version (none) auto truncate(R)(R P)
{
    struct Result
    {
        R inner = P.save;
        alias inner this;

    }
}


/// BasisPoints.truncatePrecision for general point sets which support lg_length and precision.
auto truncatePrecision(R)(R P)
{
    static if (is (R == BasisPoints))
    {
        return P.truncatePrecision();
    }
    else
    {
        //struct Result
        //{
        //    R P;
        //    alias P this;
        //    auto front()
        //    {
        //        auto ret = P.front().dup();
        //        foreach (i; 0..ret.length)
        //        {
        //            ret[i] >>= P.precision - P.lg_length;
        //        }
        //        return ret;
        //    }
        //}
    }
}

ulong[] random_vector(size_t count, size_t precision)
{
    ulong[] ret;
    foreach (i; 0..count)
        ret ~= precision.random_bits();
    return ret;
}

ulong random_bits(size_t precision)
in
{
    assert (precision <= 64);
}
body
{
    ulong ret;
    ret = uniform(0UL, 1UL << 32UL);
    if (precision == 32)
        return ret;
    if (precision & 31)
    {
        ret &= (1UL << (precision & 31)) - 1;
        if (precision < 32)
            return ret;
    }
    return (ret << 32) ^ uniform(0UL, 1UL << 32UL);
}

ulong[][] random_basis(size_t dimension, size_t precision, size_t lg_length)
in
{
    assert (precision <= 64);
}
body
{
    ulong[][] ret;
    foreach (i; 0..dimension)
        ret ~= random_vector(lg_length, precision);
    return ret;
}


import std.array : split;
import std.conv : to;
import std.typecons : Tuple;
import std.string : strip;
bool lesst(DigitalNet x, DigitalNet y)
{
    return x.t < y.t;
}
bool lessw(DigitalNet x, DigitalNet y)
{
    return x.wafom < y.wafom;
}
struct DigitalNet
{
    BasisPoints ps;
    double wafom;
    ulong t;
}

DigitalNet lineToBP(string line)
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
    return DigitalNet(BasisPoints(basis, basis.guess_precision()), wafom, t);
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

debug (lineToBP) unittest
{
    auto c = 0;
    foreach (x; "5,0.002124192556608,5.236969948020973,,2600265188,692020818,1829963221,894032275,1090497089,651123054,2898340559,1909687544,843513215,1542217271,39519261,3977641622,,2144888475,2941401343,1387697674,1986117176,3702571292,2647056038,3871827325,2263216594,3008901273,4224148358,3048652205,3799831373,,737302895,1233368001,1654098828,2764743171,239054234,249267380,1039474368,3378013260,2468295934,902812364,993745693,2410603677,,3726908047,3018079636,1719761848,2421945980,8259646,1793582138,3611200899,137680621,2493595579,2004711502,1809926346,2378246536\n".
        lineToBP().ps)
    {
        "%s".writefln(x);
        c += 1;
    }
    "unittest passed with %d elements".writefln(c);
}
