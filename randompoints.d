module randompoints;

debug import std.stdio;

import std.random : uniform;

import graycode;

struct RandomPoints
{
    immutable size_t dimension;
    immutable ulong length;
    immutable size_t bits;
    alias length opDollar;
    private size_t _position;
    private ulong[][] basis;
    private ulong[] current;
    @property size_t position()
    {
        return _position;
    }
    this(size_t dimension, size_t bits)
    {
        this.dimension = dimension;
        this.length = 1UL << bits;
        this.bits = bits;
        this._position = 0;
        this.current.length = dimension;
        this.basis = dimension.random_basis(bits);
    }
    this(RandomPoints other)
    {
        this.dimension = other.dimension;
        this.length = other.length;
        this.bits = other.bits;
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
    @property RandomPoints save()
    {
        return RandomPoints(this);
    }
}

ulong[][] random_basis(size_t dimension, size_t bits)
{
    ulong[][] ret;
    ret.length = dimension;
    foreach (i; 0..dimension)
    {
        ret[i].length = bits;
        foreach (j; 0..bits)
        {
            ret[i][j] = uniform(0UL, 1UL << 32UL) << 32UL ^ uniform(0UL, 1UL << 32UL);
            if (bits != 64)
            {
                ret[i][j] &= (1UL << bits) - 1;
            }
        }
    }
    return ret;
}

debug unittest
{
    foreach (x; RandomPoints(3, 4))
    {
        x.writeln();
    }
    readln();
}
