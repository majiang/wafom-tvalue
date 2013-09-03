module sobol;

debug import std.stdio;
import pointset : transposedBasisPoints;

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

debug = verbose;
unittest
{
    debug (verbose) "default sobol sequence of dimension 3 and precision 4:".writeln();
    foreach (x; defaultSobols(3, 4, 4))
    {
        debug (verbose) x.writeln();
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

unittest
{
    assert (DirectionNumbersGenerator([1UL, 3, 7], (1 << 3) + (1 << 1) + 1).generate(6) == [1, 3, 7, 5, 7, 43]);
    assert (DirectionNumbersGenerator([1UL, 3, 7], (1 << 3) + (1 << 1) + 1).generate(2) == [1, 3]);
    debug (verbose) "direction_numbers: unittest passed!".writeln();
}


import std.array : split;
import std.conv : to;

/// Provide direction numbers from kuo: http://web.maths.unsw.edu.au/~fkuo/sobol/new-joe-kuo-6.21201
auto defaultDirectionNumbers(immutable size_t dimension)
{
    static dn = import("sobol.csv").split();
    auto buf = dn[dimension].split(",");
    ulong primitive_polynomial = buf[0].to!ulong();
    ulong[] initial_terms;
    foreach (w; buf[1..$])
        initial_terms ~= w.to!ulong();
    return DirectionNumbersGenerator(initial_terms, primitive_polynomial);
}

version (old){
/// Provide direction numbers from kuo: http://web.maths.unsw.edu.au/~fkuo/sobol/joe-kuo-old.1111
auto defaultDirectionNumbers(immutable size_t dimension)
in
{
    assert (dimension <= 10);
}
body
{
    ulong[] initial_terms;
    foreach (x; "1,,1,1,,1,3,7,,1,1,5,,1,3,1,1,,1,1,3,7,,1,3,3,9,9,,1,3,7,13,3,,1,1,5,11,27,,1,3,5,1,15".split(",,")[dimension].split(","))
    {
        initial_terms ~= x.to!ulong();
    }
    auto primitive_polynomial = "3,7,11,13,19,25,37,59,47,61".split(",")[dimension].to!ulong();
    return DirectionNumbersGenerator(initial_terms, primitive_polynomial);
}}


/// Generate point set for the direction numbers specified.
auto sobols(ulong[][] direction_numbers)
{
    ulong[][] _direction_numbers;
    _direction_numbers.length = direction_numbers.length;
    foreach (i, c; direction_numbers)
    {
        _direction_numbers[i] = c.shift();
    }
    return transposedBasisPoints(_direction_numbers, _direction_numbers[0].length);
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
    debug (verbose) "shift: unittest passed!".writeln();
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
    debug (verbose) "degree: unittest passed!".writeln();
}
