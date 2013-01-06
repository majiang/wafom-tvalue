module sobol;

debug import std.stdio;

struct Sobols
{
    immutable size_t dimension;
    immutable ulong length;
    immutable size_t bits;
    alias length opDollar;
    private ulong[][] direction_numbers;
    private size_t _position;
    @property size_t position()
    {
        return _position;
    }
    private ulong[] current;
    this(ulong[][] direction_numbers)
    {
        this.dimension = direction_numbers.length;
        this.bits = direction_numbers[0].length;
        this.length = 1 << this.bits;
        this.current.length = this.direction_numbers.length = this.dimension;
        foreach (i, c; direction_numbers)
        {
            this.direction_numbers[i] = c.shift();
            assert (1 << c.length == this.length);
        }
        this._position = 0;
    }
    this(Sobols other)
    {
        this.dimension = other.dimension;
        this.length = other.length;
        this.bits = other.bits;
        this.current.length = other.current.length;
        this.current[] = other.current[];
        this.direction_numbers = other.direction_numbers;
        this._position = other._position;
    }
    @property ulong[] front()
    {
        return this.current;
    }
    void popFront()
    {
        this._position += 1;
        if (this.empty())
        {
            return;
        }
        auto j = this.position.bottom_zeros();
        foreach (i, c; this.direction_numbers)
        {
            this.current[i] ^= c[j];
        }
    }
    @property bool empty()
    {
        return this.length <= this.position;
    }
    @property Sobols save()
    {
        return Sobols(this);
    }
}

import std.range;
static assert (isInputRange!Sobols);
static assert (isForwardRange!Sobols);
static assert (!isBidirectionalRange!Sobols);
static assert (!isRandomAccessRange!Sobols);


/** Sobol sequence generator: one-dimensional version.

This is an input range: foreach (ulong x; Sobol(ulong[])) works well.
*/
struct Sobol
{
    ulong[] direction_numbers;
    size_t position;
    ulong current;
    this(ulong[] direction_numbers)
    {
        this.direction_numbers = direction_numbers.shift();
        this.position = 0;
        this.current = 0;
    }
    ulong front()
    {
        return this.current;
    }
    void popFront()
    {
        this.position += 1;
        if (this.empty())
        {
            return;
        }
        this.current ^= direction_numbers[this.position.bottom_zeros()];
    }
    bool empty()
    {
        return 1 << this.direction_numbers.length <= this.position;
    }
}

private ulong[] shift(ulong[] x)
{
    ulong[] ret;
    ret.length = x.length;
    foreach (i, c; x)
    {
        ret[i] = x[i] << (x.length - 1 - i);
    }
    return ret;
}

unittest
{
    assert (shift([1, 1, 1]) == [4, 2, 1]);
    debug "shift: unittest passed!".writeln();
}

/** Compute direction numbers from initial terms and a primitive polynomial which represents the recurrence relation.
*/
ulong[] direction_numbers(ulong[] initial_terms, ulong primitive_polynomial, size_t length)
{
    //debug writeln("direction_numbers(initial_terms=", initial_terms, ", primitive_polynomial=", primitive_polynomial, ", length=", length);
    auto ret = initial_terms;
    ret.length = length;
    auto degree = primitive_polynomial.degree();
    foreach (i; (initial_terms.length)..length)
    {
        foreach (j; 0..degree)
        {
            if (primitive_polynomial >> j & 1)
            {
                //debug writeln("i = ", i, "; j = ", j);
                ret[i] ^= ret[i - degree + j];
            }
            ret[i] <<= 1;
        }
        ret[i] ^= ret[i - degree];
    }
    return ret;
}

unittest
{
    assert (direction_numbers([1, 3, 7], (1 << 3) + (1 << 1) + 1, 6) == [1, 3, 7, 5, 7, 43]);
    debug "direction_numbers: unittest passed!".writeln();
}

private size_t bottom_zeros(ulong x)
{
    assert (x);
    size_t ret;
    while ((x & 1) == 0)
    {
        x >>= 1;
        ret += 1;
    }
    return ret;
}

unittest
{
    auto p = [1, 2, 3, 4, 5, 6, 7, 8, 16, 32];
    auto q = [0, 1, 0, 2, 0, 1, 0, 3, 4, 5];
    foreach (i, x; p)
    {
        assert (x.bottom_zeros() == q[i]);
    }
    debug "bottom_zeros: unittest passed!".writeln();
}

private size_t degree(ulong polynomial)
{
    assert (polynomial != 0);
    size_t ret = 0;
    while (polynomial > 1)
    {
        ret += 1;
        polynomial >>= 1;
    }
    return ret;
}

unittest
{
    auto p = [1, 2, 3, 4, 5, 6, 7, 8, 9];
    auto q = [0, 1, 1, 2, 2, 2, 2, 3, 3];
    foreach (i, x; p)
    {
        assert (x.degree() == q[i]);
    }
    debug "degree: unittest passed!".writeln();
}

debug unittest
{
    "Sobol sequence formed from direction numbers [1, 3, 7, 5]: ".write;
    foreach (x; Sobol([1, 3, 7, 5]))
    {
        x.write(", ");
    }
    writeln();
}
