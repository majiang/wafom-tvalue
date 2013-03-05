module sobol;

debug import std.stdio;
//import graycode : bottom_zeros;
import pointset : BasisPoints;

auto defaultSobols(immutable size_t dimension, immutable size_t precision, immutable size_t lg_length)
{
    assert (precision <= lg_length);
    ulong[][] _direction_numbers;
    auto trailing_zeros = precision - lg_length;
    foreach (i; 0..dimension)
    {
        _direction_numbers ~= defaultDirectionNumbers(i).generate(lg_length).all_left_shift(trailing_zeros);
    }
    return sobols(_direction_numbers);
}

debug unittest
{
    "default sobol sequence of dimension 3 and precision 4:".writeln();
    foreach (x; defaultSobols(3, 4, 4))
    {
        x.writeln();
    }
}


struct DirectionNumbersGenerator
{
    ulong[] initialTerms;
    ulong primitivePolynomial;
    /** Compute direction numbers from initial terms and a primitive polynomial which represents the recurrence relation.
    */
    ulong[] generate(immutable size_t length)
    {
        auto initial_terms = this.initialTerms;
        auto primitive_polynomial = this.primitivePolynomial;
        auto ret = initial_terms;
        ret.length = length;
        auto degree = primitive_polynomial.degree();
        foreach (i; (initial_terms.length)..length)
        {
            foreach (j; 0..degree)
            {
                if (primitive_polynomial >> j & 1)
                {
                    ret[i] ^= ret[i - degree + j];
                }
                ret[i] <<= 1;
            }
            ret[i] ^= ret[i - degree];
        }
        return ret;
    }
}


ulong[] all_left_shift(ulong[] x, immutable size_t trailing_zeros)
{
    auto ret = new ulong[x.length];
    foreach (i, c; x)
    {
        ret[i] = c << trailing_zeros;
    }
    return ret;
}

import std.array : split;
//import std.algorithm : map;
import std.conv : to;

private auto defaultDirectionNumbers(immutable size_t dimension)
in
{
    assert (dimension <= 10);
}
body
{
    ulong[] initial_terms;
    foreach (x; "1,,1,1,,1,3,7,,1,1,5,,1,3,1,1,,1,1,3,7,,1,3,3,9,9,,1,3,7,13,3,,1,1,5,11,27,,1,3,5,1,15"
        .split(",,")[dimension]
        .split(","))
    {
        initial_terms ~= x.to!ulong();
    }
    auto primitive_polynomial = "3,7,11,13,19,25,37,59,47,61"
        .split(",")[dimension]
        .to!ulong();
    return DirectionNumbersGenerator(initial_terms, primitive_polynomial);
}

/// Generate point set for the direction numbers specified.
auto sobols(ulong[][] direction_numbers)
{
    ulong[][] _direction_numbers;
    _direction_numbers.length = direction_numbers.length;
    foreach (i, c; direction_numbers)
    {
        _direction_numbers[i] = c.shift();
    }
    return BasisPoints(_direction_numbers, _direction_numbers[0].length);
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
    assert ([1, 1, 1].shift() == [4, 2, 1]);
    debug "shift: unittest passed!".writeln();
}


unittest
{
    assert (DirectionNumbersGenerator([1UL, 3, 7], (1 << 3) + (1 << 1) + 1).generate(6) == [1, 3, 7, 5, 7, 43]);
    assert (DirectionNumbersGenerator([1UL, 3, 7], (1 << 3) + (1 << 1) + 1).generate(2) == [1, 3]);
    debug "direction_numbers: unittest passed!".writeln();
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

version (old){
    /** Sobol Sequence generator: Multidimensional version.

    Examples:
    ----------------
    double integral = 0;
    auto sobol = Sobol([direction_numbers(...), ...]);
    foreach (x; sobol) integral += x.f();
    ----------------

    Remarks:
    t <= (primitive).(degree-1).sum
    */
    struct Sobols
    {
        immutable size_t dimension;
        immutable ulong length;
        immutable size_t precision;
        alias length opDollar;
        private ulong[][] direction_numbers;
        private size_t _position;
        @property size_t position()
        {
            return _position;
        }
        private ulong[] current;
        this(ulong[][] direction_numbers, size_t precision)
        {
            this.dimension = direction_numbers.length;
            this.precision = precision;
            this.length = 1UL << direction_numbers[0].length;
            this.current.length = this.direction_numbers.length = this.dimension;
            foreach (i, c; direction_numbers)
            {
                this.direction_numbers[i] = c.shift();
                assert (1UL << c.length == this.length);
            }
            this._position = 0;
        }
        this(Sobols other)
        {
            this.dimension = other.dimension;
            this.length = other.length;
            this.precision = other.precision;
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
    debug unittest
    {
        "Sobol sequence formed from direction numbers [1, 3, 7, 5]: ".write;
        foreach (x; Sobol([1, 3, 7, 5]))
        {
            x.write(", ");
        }
        writeln();
    }
}
