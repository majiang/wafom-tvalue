module pointset;

debug import std.stdio;

import std.random : uniform;
public import sobol : sobols, direction_numbers;
import graycode;

/** Generate a point set of dimension, precision and lg(length) specified by choosing its basis randomly.
*/
auto randomPoints(size_t dimension, size_t precision, size_t lg_length)
{
    return BasisPoints(dimension.random_basis(precision, lg_length), precision);
}

unittest
{
    debug foreach (x; randomPoints(3, 10, 4))
    {
        x.writeln();
    }
    "Press Enter to continue:".writeln();
    readln();
}

/** Input Range for point set.

Generate linear combinations of basis using gray code algorithm.

TODO: dimension over F2 < precision
*/
struct BasisPoints
{
    immutable size_t dimension;
    immutable ulong length;
    immutable size_t precision;
    alias length opDollar;
    private size_t _position;
    private ulong[][] basis;
    private ulong[] current;
    @property size_t position()
    {
        return _position;
    }
    this(ulong[][] basis, size_t precision)
    {
        this.dimension = basis.length;
        this.precision = precision;
        this.length = 1UL << basis[0].length;
        this._position = 0;
        this.current.length = this.dimension;
        this.basis = basis;
    }
    this(BasisPoints other)
    {
        this.dimension = other.dimension;
        this.length = other.length;
        this.precision = other.precision;
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
}

private ulong[][] random_basis(size_t dimension, size_t precision, size_t lg_length)
{
    ulong[][] ret;
    ret.length = dimension;
    foreach (i; 0..dimension)
    {
        ret[i].length = lg_length;
        foreach (j; 0..lg_length)
        {
            ret[i][j] = uniform(0UL, 1UL << 32UL) << 32UL ^ uniform(0UL, 1UL << 32UL);
            if (precision != 64)
            {
                ret[i][j] &= (1UL << precision) - 1;
            }
        }
    }
    return ret;
}
