module pointset;

debug import std.stdio;

import std.random : uniform;
public import sobol : defaultSobols;
import graycode;

/** Generate a point set of dimension, precision and lg(length) specified by choosing its basis randomly.
*/
auto randomPoints(immutable size_t dimension, immutable size_t precision, immutable size_t lg_length)
{
    return BasisPoints(dimension.random_basis(precision, lg_length), precision);
}

debug unittest
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
    /*private /**/ulong[][] basis;
    private ulong[] current;
    @property ulong position()
    {
        return _position;
    }
    this(ulong[][] basis, size_t precision) // TODO: basis[i][j] -> basis[j][i] is better?
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
    this(BasisPoints other)
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

/// BasisPoints.truncatePrecision for general point sets which support lg_length and precision.
auto truncatePrecision(R)(R P)
{
    static if (is (R == BasisPoints))
    {
        return P.truncatePrecision();
    }
    else
    {
        struct Result
        {
            R P;
            alias P this;
            auto front()
            {
                auto ret = P.front().dup();
                foreach (i; 0..ret.length)
                {
                    ret[i] >>= P.precision - P.lg_length;
                }
                return ret;
            }
        }
    }
}

private ulong[][] random_basis(size_t dimension, size_t precision, size_t lg_length)
in
{
    assert (precision <= 64);
}
body
{
    ulong[][] ret;
    ret.length = dimension;
    foreach (i; 0..dimension)
    {
        ret[i].length = lg_length;
        foreach (j; 0..lg_length)
        {
            ret[i][j] = uniform(0UL, 1UL << 32UL) << 32UL ^ uniform(0UL, 1UL << 32UL);
            if (precision < 64)
            {
                ret[i][j] &= (1UL << precision) - 1;
            }
        }
    }
    return ret;
}
