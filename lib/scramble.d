module lib.scramble;

import std.random : uniform;
import std.algorithm : swap;
import std.exception : enforce;

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
        {x ^= x >> 32; x >>= 0x100000000UL - 1;}
    static if (is (U == ulong) || is (U == uint))
        {x ^= x >> 16; x >>= 0x10000U - 1;}
    static if (is (U == ulong) || is (U == uint) || is (U == ushort))
        {x ^= x >> 8; x >>= 0x100U - 1;}
    x ^= x >> 4; x &= 0x10U - 1;
    x ^= x >> 2; x &= 4-1;
    x ^= x >> 1; x &= 2-1;
    return x;
}

bool isLinearlyIndependent(U)(U[] a)
{
    import std.array;
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
        x = x[0..$-1];
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
    ret.writeln();
    return ret;
}

/** Multiply a matrix from the left to each coordinate.
*/
auto multiplyMatrices(U)(U[][] basis, U[][] matrices)
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
auto multiplyRandomMatrices(U)(U[][] basis, size_t precision, size_t distance)
{
    immutable dimB = enforce(basis.length);
    immutable dimR = enforce(basis[0].length);
    U[][] matrices;
    auto im = (precision - distance).identityMatrix!U();
    foreach (i; 0..dimR)
        matrices ~= distance.randomNonsingularMatrix!U().directSum(im);
    return basis.multiplyMatrices(matrices);
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

void notmain()
{
    import std.algorithm;
    import std.array;
    foreach (i; 0..10)
    {
        auto a = 4.randomMatrix!ubyte();
        a.map!(r => [r]).array().output();
        a.isLinearlyIndependent().writeln();
    }
}

void main()
{
    ubyte[][] b = [[128, 255], [0, 128], [128, 127]];
    b.output();
    foreach (i; 0..9)
    {
        writeln(i);
        foreach (j; 0..10)
            b.multiplyRandomMatrices(8, i).output();
    }
}

auto output(U)(U[][] basis)
{
    "%(%(%08b %)\n%)\n".writefln(basis);
    return basis;
}
