module tvalue;

import std.functional : memoize;

debug import std.stdio;
import std.bigint;

version = 1;
version (1) {alias tvalue1 tvalue;}
version (2) {alias tvalue2 tvalue;}


/** Compute t-value of a digital net formed from Sobol sequences.

Algorithm:
Algorithm 1 of MacWilliams.
* 1: input P = (x[n] : 0 <= n < b^m) in [0..1)^s.
** b = 2
** m = P.precision
** s = P.dimension
*
* 2: compute the coefficients N[a] of z^a for 0 <= a <= m of the polynomial
*     Q[m](z)(b^-m)sum[n in b^m]prod[i in s](1-(bz)^nu*(x[n,i]))
* modulo z^(m+1), where
*     Q[m](z) := (1+(b-1)z+(b^2-b)z^2 + ... + (b^m-b^{m-1})z^m)^s
* and
*     nu*(x) = ceil(-lg x)
*
* 3: t = m + 1 - min(a : 0 < a <= m; N[a] = 0)
*
* 4: return t.

Params:
P = a subset of [0..(2^m))^s whose cardinality is 2^m

Remarks:
Remark 2 of the paper states that:<blockquote>
<p>If one uses this function for a general point set, then the returned 
value of t is a lower bound on the t-value of the point set, 
i.e., it implies that P is not a (t-1, m, s)-net in base two.</p></blockquote>
but this function may have a bug for general point set.

TODO (for developper):
Check whether lower t-value bound of a point set is properly calculated
by this function. P.precision is used in some code, but P.length.lg may
be suitable.
*/
ulong tvalue1(R)(R P)
{
    auto total = new BigInt[P.precision + 1] ; // sum
    foreach (x; P)
    {
        auto current = P.precision.empty_product(); // prod
        foreach (e; x)
        {
            auto nu = e.reciprocal(P.precision);
            foreach_reverse (i; 0..(P.precision+1-nu))
            {
                current[i+nu] -= current[i] << nu;
            }
        }
        total = total.plus(current);
    }
    // construct Q.
    // TODO: compare with power_part.
    BigInt[] Q;
    Q.length = P.precision + 1;
    Q[0] = 1;
    foreach (i; 0..P.precision)
    {
        Q[i+1] = 1;
        Q[i+1] <<= i;
    }
    foreach (i, x; Q.power(P.dimension).times(total)) // (b^-m) is unnecessary
    {
        if (i == 0)
        {
            continue;
        }
        if (x.toLong)
        {
            return P.precision + 1 - i;
        }
    }
    return 0; // all zero, t = m + 1 - (m+1) according to Algorithm 1.
}

/** Compute t-value of a digital net formed from Sobol sequences.

Algorithm:
Algorithm 2 of MacWilliams.
* 1: input P = (x[n] : 0 <= n < b^m) in [0..1)^s
** b = 2
** m = P.precision
** s = P.dimension
*
* 2: compute the coefficients of z^a for (s-1)(m+1) <= a < s(m+1) of the polynomial
*     Q(z) = b^{(s-1)m}sum[n in b^m]prod[i in s](z^mu(x[n,i])-z^(m+1))
*                 - (1+(b-1)z+...+(b^m-b^{m-1})-b^mz^{m+1})^s
* where mu(x) = 1+floor(lg(b^m x)) = m+1 + floor(lg x) = m+1 - ceil(-lg x)
*
* 3: t = (1 - s)(m + 1) + Q.degree
*
* 4: return t.

This function computes only the m+1 most significant coefficients
by computing the m+1 lowest coefficients of the product of the reciprocal
prod[i in s](-1 + y^(m+1-mu(x))) = prod[i in s](-1 + y^nu*(x))

Params:
P = a subset of  [0..2^m)^s
*/
ulong tvalue2(R)(R P)
{
    auto total = new BigInt[P.precision + 1];
    foreach (x; P)
    {
        auto current = P.precision.empty_product();
        foreach (e; x)
        {
            auto nu = e.reciprocal(P.precision);
            foreach_reverse (i; nu..total.length)
            {
                current[i] = current[i - nu] - current[i];
            }
            foreach_reverse (i; 0..nu)
            {
                current[i] = -current[i];
            }
            assert (current[0].toLong); // 最上位の係数は0であってはならない。
        }
        total = total.plus(current);
        assert (total[0].toLong); // ditto
    }
    auto shift = (P.dimension - 1) * P.precision;
    // BigInt mask = (BigInt(1) << shift) - 1;
    foreach (i, c; P.precision.power_part().power(P.dimension))
    {
        if (total[i] << shift != c)
        //if ((c & mask) || total[i] != c >> shift) と書きたいが、BigInt が cast(bool) や &, |, ^ を持っていないのでできない。
        {
            assert (i);
            return P.precision + 1 - i;
        }
    }
    //assert (!(total[0].toLong)); // 最上位の係数は (ここでは) 0でなければならない。
    //foreach (i, c; total)
    //{
        //if (c.toLong)
        //{
            //return P.precision + 1 - i;
        //}
    //}
    return 0;
}


private BigInt[] power_part(immutable size_t m)
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
private size_t reciprocal(ulong x, immutable size_t m)
in
{
    assert (0 <= x);
    if (m != 64 && !(x < 1UL << m))
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

private size_t reciprocal_binary_search(ulong x, immutable size_t m)
{
    if (x == 0) return m + 1;
    size_t ret = m + 1;
    size_t current = 32;
    while (current)
    {
        if (x >> current)
        {
            ret -= current;
            x >>= current;
        }
        current >>= 1;
    }
    return ret - 1;
}

unittest
{
    "testing reciprocal: binary-search version".writeln();
    foreach (x; 0..10000)
    {
        assert (reciprocal(x, 14) == reciprocal_binary_search(x, 14));
        assert (reciprocal(x, 15) == reciprocal_binary_search(x, 15));
        assert (reciprocal(x, 16) == reciprocal_binary_search(x, 16));
        assert (reciprocal(x, 17) == reciprocal_binary_search(x, 17));
        assert (reciprocal(x, 32) == reciprocal_binary_search(x, 32));
        assert (reciprocal(x, 33) == reciprocal_binary_search(x, 33));
        assert (reciprocal(x, 60) == reciprocal_binary_search(x, 60));
        assert (reciprocal(x, 63) == reciprocal_binary_search(x, 63));
        assert (reciprocal(x, 64) == reciprocal_binary_search(x, 64));
    }
    foreach (m; 1..65)
    {
        foreach (x; 0..m)
        {
            assert (reciprocal(1UL << x, m) == reciprocal_binary_search(1UL << x, m));
            assert (reciprocal((1UL << x) - 1, m) == reciprocal_binary_search((1UL << x) - 1, m));
        }
    }
    "...OK.".writeln();
}

private BigInt[] empty_product(immutable size_t m)
{
    return [BigInt(1)] ~ new BigInt[m];
}

private BigInt[] times(BigInt[] f, BigInt[] g)
in
{
    assert (f.length == g.length);
}
body
{
    auto ret = new BigInt[f.length];
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

private BigInt[] square(BigInt[] f)
{
    auto ret = new BigInt[f.length];
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

private BigInt[] power(BigInt[] f, ulong s)
{
    auto ret = (f.length - 1).empty_product();
    BigInt[] sq = f.dup;
    while (s)
    {
        if (s & 1)
        {
            ret = ret.times(sq);
        }
        sq = sq.square();
        s >>= 1;
    }
    return ret;
}

version (old) {private BigInt[] times_sparse(BigInt[] f, size_t g) // in-place なものに置き換えた
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
    BigInt[] ret;
    ret.length = f.length;
    foreach (i, c; f)
    {
        ret[i] = -c;
    }
    foreach (i, c; f)
    {
        if (ret.length <= i + g)
        {
            break;
        }
        ret[i + g] += c;
    }
    return ret;
}}

private BigInt[] plus(BigInt[] f, BigInt[] g)
in
{
    assert (f.length == g.length);
}
body
{
    auto ret = new BigInt[f.length];
    foreach (i, c; f)
    {
        ret[i] = c + g[i];
    }
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
                current = current.times_sparse(reciprocal!(cast(size_t)m, b)(e));
            }
            total = total.plus!m(current);
        }
        total[] *= b ^^ ((s - 1) * m);
        auto rest = power_part!(b, m)().power!m(s);
        rest[] *= -1;
        foreach (i, c; total.plus!m(rest))
        {
            if (c)
            {
                return m + 1 - i;
            }
        }
        return 0;
    }}
