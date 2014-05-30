module lib.scramble;

import std.random : uniform;
import std.algorithm : swap;
import std.exception : enforce;
import std.string : format;
import std.array;
import std.conv : text;

auto deepcopy(U)(in U[][] basis)
{
    auto ret = new U[][](basis.length, basis[0].length);
    foreach (i, ref row; ret)
        row[] = basis[i][];
    return ret;
}

/// Given a space by a basis, take another basis of the space so that the first vector of the output is randomly selected from the whole space.
auto changeBasis(U)(U[][] basis)
{
    immutable m = basis.length;
    immutable k = 1.uniform(1 << m);
    bool first = true;
    foreach (i; 0..m)
        if (k >> i & 1)
            if (first)
            {
                swap(basis[0], basis[i]);
                first = false;
            }
            else
                basis[0][] ^= basis[i][];
    return basis;
}

/// Modify the first vector of the basis by a bit specified.
auto modify(U)(U[][] basis, size_t index_s, U b_idx_r)
{
    basis[0][index_s] ^= 1 << b_idx_r;
    return basis;
}

import std.traits : isUnsigned;
auto sumB(U)(U x)
    if (isUnsigned!U)
out (result)
{
    assert (result == 0 || result == 1);
}
body
{
    static if (is (U == ulong))
        {x ^= x >> 32; x &= 0x100000000UL - 1;}
    static if (is (U == ulong) || is (U == uint))
        {x ^= x >> 16; x &= 0x10000U - 1;}
    static if (is (U == ulong) || is (U == uint) || is (U == ushort))
        {x ^= x >> 8; x &= 0x100U - 1;}
    x ^= x >> 4; x &= 0x10U - 1;
    x ^= x >> 2; x &= 4-1;
    x ^= x >> 1; x &= 2-1;
    return x;
}

bool isLinearlyIndependent(U)(U[] a)
{
    import std.algorithm;
    auto x = a.dup;
    while (!x.empty)
    {
        sort(x);
        if (x[$-1] >> (x.length-1) == 0)
            return false;
        foreach (i; 0..x.length-1)
            if (x[i] >> (x.length-1))
                x[i] ^= x[$-1];
        x.length -= 1;
    }
    return true;
}

auto randomMatrix(U)(size_t size)
{
    auto ret = new U[size];
    foreach (ref row; ret)
        row = (cast(U)1).uniform!"[]"(cast(U)(U.max >> ((U.sizeof << 3) - size)));
    return ret;
}

auto randomNonsingularMatrix(U)(size_t size)
{
    while (true)
    {
        auto r = size.randomMatrix!U();
        if (r.isLinearlyIndependent())
            return r;
    }
    assert (false);
}

auto identityMatrix(U)(size_t size)
{
    U r = 1;
    auto ret = new U[size];
    foreach_reverse (i; 0..size)
    {
        ret[i] = r;
        r <<= 1;
    }
    return ret;
}

/// The direct sum of the two matrix given.
auto directSum(U)(U[] a, U[] b)
{
    auto ret = a ~ b;
    foreach (i; 0..a.length)
        ret[i] <<= b.length;
    return ret;
}

/** Multiply a matrix from the left to each coordinate.
*/
auto multiplyMatrices(U)(in U[][] basis, U[][] matrices)
    if (isUnsigned!U)
{
    immutable dimB = enforce(basis.length);
    immutable dimR = enforce(matrices.length);
    immutable precision = matrices[0].length;
    auto ret = new U[][](dimB, dimR);
    foreach (j; 0..dimR)
        foreach (k; 0..precision)
            foreach (i; 0..dimB)
                ret[i][j] ^= (cast(U)(matrices[j][k] & basis[i][j])).sumB() << (precision - k - 1);
    return ret;
}

/** Yoshiki scramble: multiplying the direct sum of a random matrix and the identity matrix from the left.

bugs: no linear independecy check */
auto multiplyRandomMatrices(U)(in U[][] basis, size_t precision, size_t distance)
{
    immutable dimB = enforce(basis.length);
    immutable dimR = enforce(basis[0].length);
    U[][] matrices;
    auto im = (precision - distance).identityMatrix!U();
    foreach (i; 0..dimR)
        matrices ~= distance.randomNonsingularMatrix!U().directSum(im);
    return basis.multiplyMatrices(matrices);
}

alias Precision = size_t;
alias DimensionR = size_t;

struct LargeMatrix(U)
    if (isUnsigned!U)
{
    immutable size_t precision, dimensionR, size;
    enum bpe = U.sizeof << 3;
    U[][] bits;

    this (Precision precision, DimensionR dimensionR)
    {
        this.precision = precision;
        assert (this.precision <= bpe);
        this.dimensionR = enforce(dimensionR);
        this.size = this.precision * this.dimensionR;
        this.bits = new U[][](size, this.dimensionR);
    }
    string toString()
    {
import std.conv : text;
        return text("%(%(%0", this.precision, "b%)\n%)").format(bits);
    }
    const bool opIndex(size_t i, size_t j)
    in
    {
        assert (i < size);
        assert (j < size);
    }
    body
    {
        return bits[i][j/precision] >> (precision - 1 - j % precision) & 1;
    }
    auto setDiagonal()
    {
        foreach (i; 0..this.size)
            this[i, i] = true;
        return this;
    }
    bool opIndexAssign(bool bit, size_t i, size_t j)
    in
    {
        assert (i < size);
        assert (j < size);
    }
    body
    {
        if (bit)
            return !!(bits[i][j/precision] |= cast(U)1 << (precision - 1 - j % precision));
        else
            return !!(bits[i][j/precision] &= ~(cast(U)1 << (precision - 1 - j % precision)));
    }
    U[] opBinary(string op)(U[] vector)
        if (op == "*")
    {
        enforce(vector.length == dimensionR);
        auto ret = new U[dimensionR];
        foreach (r, row; bits)
        {
            immutable size_t s = r / precision, n = precision - 1 - r % precision;
            U tmp = 0;
            foreach (j, e; row)
                tmp ^= e & vector[j];
            ret[s] |= tmp.sumB() << n;
        }
        return ret;
        }
    U[][] opBinary(string op)(U[][] basis)
        if (op == "*")
    {
        U[][] ret;
        foreach (vector; basis)
            ret ~= this * vector;
        return ret;
    }
}

auto getBit(U)(U[] bitseq, size_t index, Precision precision)
{
    //"getBit(%d, %d) = [%d] >> %d & 1".writefln(index, precision, bitseq.length - 1 - index/precision, index % precision);
    return !!(bitseq[bitseq.length - 1 - index/precision] >> (index % precision) & 1);
}

auto isNonsingular(U)(LargeMatrix!U m)
{
    import std.algorithm;
    auto x = m.bits.deepcopy();
    while (!x.empty)
    {
        sort(x);
        //x.outputMatrix(m.precision);
        if (!x[$-1].getBit(x.length-1, m.precision))
        {
            //"returning false".writeln();readln();
            return false;
        }
        auto last_row = x[$-1];
        x.length -= 1;
        foreach (i, ref r; x)
            if (r.getBit(x.length, m.precision))
                foreach (j, ref e; r)
                    e ^= last_row[j];
    }
    //"returning true".writeln();readln();
    return true;
}

auto largeRandomMatrix(U)(Precision precision, DimensionR dimensionR)
{
    auto ret = LargeMatrix!U(precision, dimensionR);
    auto size = precision * dimensionR;
    foreach (i; 0..size)
        foreach (j; 0..size)
            ret[i, j] = uniform(0, 2) == 1;
    return ret;
}

auto largeRandomNonsingularMatrix(U)(Precision precision, DimensionR dimensionR)
{
    while (true)
    {
        auto ret = largeRandomMatrix!U(precision, dimensionR);
        if (ret.isNonsingular())
            return ret;
    }
    assert (false);
}

auto largeRandomIdentity(U)(Precision precision, DimensionR dimensionR, size_t degree_of_scrambling)
{
    auto I = LargeMatrix!U(precision - degree_of_scrambling, dimensionR);
    I.setDiagonal();
    return largeRandomNonsingularMatrix!U(degree_of_scrambling, dimensionR).directSum(I);
}

auto directSum(U)(LargeMatrix!U a, LargeMatrix!U b)
{
    enforce(a.dimensionR == b.dimensionR);
    immutable
        precision = a.precision + b.precision,
        dimR = a.dimensionR;
    auto ret = LargeMatrix!U(a.precision + b.precision, a.dimensionR);
    // a
    foreach (s; 0..dimR)
        foreach (t; 0..a.precision)
            foreach (u; 0..dimR)
                foreach (v; 0..a.precision)
                {
//                    debug "ret[%d, %d] = a[%d, %d]".writefln(s * precision + t, u * precision + v, s * a.precision + t, u * a.precision + v);
                    ret[s * precision + t, u * precision + v]
                    = a[s * a.precision + t, u * a.precision + v];
//                    debug readln();
                }
    // b
    foreach (s; 0..dimR)
        foreach (t; 0..b.precision)
            foreach (u; 0..dimR)
                foreach (v; 0..b.precision)
                {
//                    debug "ret[%d, %d] = b[%d, %d]".writefln(s * precision + a.precision + t, u * precision + a.precision + v, s * b.precision + t, u * b.precision + v);
                    ret[s * precision + a.precision + t, u * precision + a.precision + v]
                    = b[s * b.precision + t, u * b.precision + v];
//                    debug readln();
                }
    return ret;
}

version (old)
auto multiplyRandomMatrices(U)(U[][] basis, size_t precision, size_t distance)
{
    enforce(distance <= precision);
    immutable dimB = enforce(basis.length);
    immutable dimR = enforce(basis[0].length);
    immutable U maskM = cast(U)(U.max >> ((U.sizeof << 3) - distance) << (precision - distance));
    immutable U maskI = U.max >> ((U.sizeof << 3) + distance - precision);

    auto ret = new U[][](dimB, dimR);
    foreach (j; 0..dimR)
    {
        foreach (k; 0..distance)
        {
            auto row = U.min.uniform!"[]"(U.max) & maskM;
            foreach (i; 0..dimB)
                ret[i][j] ^= (row & basis[i][j]).sumB() << (precision - distance + k);
        }
        foreach (i; 0..dimB)
            ret[i][j] ^= basis[i][j] & maskI;
    }
    return ret;
}

version (stand_alone):
import std.stdio;

void main()
{
    testLMall();
}

auto testLM(U)()
{
    enum size_t size = 48, precision = 6 * U.sizeof, dimR = size / precision;
    foreach (d_of_scrambling; [3, 6])
    {
        largeRandomIdentity!U(precision, dimR, d_of_scrambling).toString().writeln();
        writeln();
    }
}

auto oldtestLM(U)()
{
    "in".writeln();
    enum size_t size = 48, precision = 6 * U.sizeof, dimR = size / precision;
    enum U mask = U.max >> ((U.sizeof << 3) - precision);
    enum num1 = 10;
    auto m = LargeMatrix!U(precision, dimR).setDiagonal();
    assert (m.isNonsingular());
    foreach (k; 0..num1)
        m[0.uniform(size), 0.uniform(size)] = true;
    "enter while".writeln();
    while (!m.isNonsingular())
    {
        auto i = 0.uniform(size), j = 0.uniform(size);
        m[i, j] = !m[i, j];
    }
    m.toString().writeln();
    auto v = new U[dimR];
    foreach (i; 0..v.length)
        v[i] = (U.max >> ((U.sizeof << 3) + i - precision));
    text("%(%0", precision, "b %)").writefln(v);
    text("%(%0", precision, "b %)").writefln(m * v);

    enum numtry = 10000;
    double countnonsingular = 0;
    m.bits = LargeMatrix!U(precision, dimR).bits;
    foreach (i; 0..numtry)
    {
        foreach (ref row; m.bits)
            foreach (ref e; row)
                e ^= 0.uniform!("[]", U, U)(mask);
        if (m.isNonsingular())
            countnonsingular += 1;
    }
    writeln("nonsingular probability = ", countnonsingular / numtry);
}

void testLMall()
{
    import std.stdio;
    testLM!ulong();
    testLM!uint();
    testLM!ushort();
    testLM!ubyte();
}

auto output(U)(U[][] basis)
{
    "%(%(%08b %)\n%)\n".writefln(basis);
    return basis;
}

void outputMatrix(U)(U[][] bits, size_t precision)
{
    text("%(%(%0", precision, "b%)\n%)").format(bits).writeln();writeln();
}
