module tvalue;

import sobol : Sobols;
debug import std.stdio;
import std.bigint;

ulong tvalue1(Sobols P)//, immutable size_t bits, immutable size_t dimension)
{
    BigInt[] total;
    total.length = P.bits + 1;
    foreach (x; P)
    {
        auto current = P.bits.empty_product();
        foreach (e; x)
        {
            auto nu = e.reciprocal(P.bits);
            foreach_reverse (i; 0..(P.bits+1))
            {
                if (P.bits < i + nu)
                {
                    continue;
                }
                current[i+nu] -= current[i] << nu;
            }
        }
        total = total.add(current);
    }
    BigInt[] Q;
    Q.length = P.bits + 1;
    Q[0] = 1;
    foreach (i; 0..P.bits)
    {
        Q[i+1] = 1;
        Q[i+1] <<= i;
    }
    Q = Q.polynomial_power(P.dimension);
    foreach (i, x; Q.polynomial_product_polynomial(total))
    {
        if (i == 0)
            continue;
        if (x.toLong)
            return P.bits + 1 - i;
            //return bits - i;
    }
    return -1;
}

/** Compute t-value of a digital net formed from Sobol sequences.

P = a subset of  [0..(1<<m))^s
*/
ulong tvalue(Sobols P)//, immutable size_t bits, immutable size_t dimension)
{
    BigInt[] total;
    total.length = P.bits + 1;
    int j = 0;
    foreach (x; P)
    {
        auto current = P.bits.empty_product();
        foreach (e; x)
        {
            assert (current[0].toLong);
            current = current.polynomial_product_sparse(e.reciprocal(P.bits));
        }
        assert (current[0].toLong);
        total = total.add(current);
        //debug j.write(": ");
        assert (total[0].toLong);
        //debug "OK".writeln();
        j++;
    }
    foreach (i; 0..(P.bits + 1)) total[i] <<= ((P.dimension - 1) * P.bits);
    //debug "total = ".writeln(total);
    auto rest = P.bits.power_part().polynomial_power(P.dimension);
    foreach (i; 0..(P.bits + 1)) rest[i] = -rest[i];
    //debug "rest = ".writeln(total);
    foreach (i, c; total.add(rest))
    {
        //debug c.write(", ");
        if (c.toLong)
        {
            assert (i);
            return P.bits + 1 - i;
        }
    }
    return 0;
}


BigInt[] power_part(immutable size_t m)
{
    BigInt[] ret;
    ret.length = m + 1;
    ret[0] = -1;
    ret[0] <<= m;
    foreach (i; 1..m+1)
    {
        ret[i] = 1;
        ret[i] <<= m - i;
    }
//    debug writeln(ret);
    return ret;
}

/** Construct reciprocal polynomial.
*
* Params:
*     m = lg|P|
*     x = a component of a vector in P; 0 <= x < 2^m
* Returns:
*     m + 1 if x = 0
*     m - (x.lg.floor) otherwise
*/
size_t reciprocal(ulong x, immutable size_t m)
in
{
    assert (0 <= x);
    if (!(x < 1UL << m))
    {
        x.writeln(".reciprocal illegally called where m = ", m);
        assert (false);
    }
}
out (result)
{
    assert (0 < result);
    assert (result <= m+1);
}
body
{
    size_t ret = m + 1;
    while (x)
    {
        ret -= 1;
        x >>= 1;
    }
    return ret;
}


BigInt[] empty_product(immutable size_t m)
{
    BigInt[] ret;
    ret.length = m + 1;
    ret[0] = 1;
    return ret;
}

BigInt[] polynomial_product_polynomial(BigInt[] f, BigInt[] g)
{
    BigInt[] ret;
    ret.length = f.length;
    foreach (i, c; f)
    {
        foreach (j, d; g)
        {
            if (ret.length <= i + j)
            {
                break;
            }
            ret[i + j] += c * d;
        }
    }
    return ret;
}

BigInt[] polynomial_square(BigInt[] f)
{
    BigInt[] ret;
    ret.length = f.length;
    foreach (i, c; f)
    {
        if (i << 1 < ret.length)
        {
            ret[i << 1] = c * c;
        }
        foreach (j; 0..i)
        {
            if (i + j < ret.length)
            {
                ret[i + j] += (c * f[j]) << 1;
            }
            else
            {
                break;
            }
        }
    }
    return ret;
}

BigInt[] polynomial_power(BigInt[] f, ulong s)
{
    auto ret = (f.length - 1).empty_product();
    BigInt[] square = f.dup;
    while (s)
    {
        if (s & 1)
        {
            ret = ret.polynomial_product_polynomial(square);
        }
        square = square.polynomial_square();
        s >>= 1;
    }
    return ret;
}

BigInt[] polynomial_product_sparse(BigInt[] f, size_t g)
in
{
    assert (f[0].toLong);
    assert (g);
}
out (result)
{
    assert (result[0].toLong);
}
body
{
    //debug f.writeln(" times -1+z^", g); 
    BigInt[] ret;
    ret.length = f.length;
    foreach (i, c; f)
    {
        ret[i] = -c;
    }
    //debug ret.writeln();
    foreach (i, c; f)
    {
        if (ret.length <= i + g)
        {
            break;
        }
        ret[i + g] += c;
    }
    //debug ret.writeln();
    return ret;
}

BigInt[] add(BigInt[] f, BigInt[] g)
{
    BigInt[] ret;
    ret.length = f.length;
    foreach (i, c; f)
    {
        ret[i] = c + g[i];
    }
    //ret[] = f[] + g[];
    return ret;
}

version (digitalnet){
    /** Compute t-value of digital net.

    Algorithm 2 of MacWilliams12.pdf

    Params:
    b = |G|
    P = digital net, subset of [0..b^m) ^ s

    Returns:
    t-value
    */
    ulong tvalue(ulong b, ulong m, ulong s)(ulong[s][b ^^ m] P)
    {
        ulong[m+1] total;
        ulong[m+1] current;
        foreach (x; P)
        {
            current = empty_product!m();
            foreach (e; x)
            {
                current = current.polynomial_product_sparse(reciprocal!(cast(size_t)m, b)(e));
            }
            total = total.add!m(current);
        }
        total[] *= b ^^ ((s - 1) * m);
        auto rest = power_part!(b, m)().polynomial_power!m(s);
        rest[] *= -1;
        foreach (i, c; total.add!m(rest))
        {
            if (c)
            {
                return m + 1 - i;
            }
        }
        return 0;
    }}
