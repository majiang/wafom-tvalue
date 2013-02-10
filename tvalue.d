module tvalue;

import std.functional : memoize;

debug import std.stdio;
import std.bigint;

version = 1;
version (1) {alias tvalue1 tvalue;}
version (2) {alias tvalue2 tvalue;}


/** Compute t-value of a digital net.

* Algorithm:
* Algorithm 1 of MacWilliams.<ol>
*   <li>input P = (x[n] : 0 <= n < b<sup>m</sup>), a subset of [0..1)<sup>s</sup>. This function requires P to have three properties: lg_length, precision and dimension. Here,<ul>
*       <li>b = 2</li>
*       <li>m = P.lg_length</li>
*       <li>s = P.dimension</li>
*   </ul></li>
*   <li>compute the coefficients N[a] of z<sup>a</sup> for 0 <= a <= m of the polynomial<ul>
*       <li>b<sup>-m</sup>Q[m](z)sum[x in P]prod[i in s](1-(bz)<sup>nu(x[i])</sup>)</li>
*   </ul>modulo z<sup>(m+1)</sup>, where<ul>
*       <li>Q[m](z) := (1 + (b-1)z + (b<sup>2</sup>-b)z<sup>2</sup> + ... + (b<sup>m</sup>-b<sup>m-1</sup>)z<sup>m</sup>)<sup>s</sup> = (1 + z + 2z<sup>2</sup> + ... + b<sup>m-1</sup>z<sup>m</sup>)<sup>s</sup></li>
*   </ul>and<ul>
*       <li>nu(x) = ceil(-lg x)</li>
*   </ul></li>
*   <li>return t := m + 1 - min{a : 0 < a <= m; N[a] = 0}</li>
* </ol>

* Params:
* P = an InputRange which yields ulong[]s, and have properties length, lg_length, precision, dimension.

* Conditions:
* P must satisfy three conditions:<ul>
*   <li>1UL << P.lg_length == P.length</li>
*   <li>P.front.length == P.dimension</li>
*   <li>0 <= P.front < (1 << P.precision)</li>
* </ul>

* Remarks:
* Remark 2 of the paper is the following:<blockquote>
* <p>If one uses this function for a general point set, then the returned value of t is a lower bound on the t-value of the point set, i.e., it implies that P is not a (t-1, m, s)-net in base two.</p></blockquote>
*/
ulong tvalue1(R)(R P)
{
    if (P.lg_length < P.precision)
    {
        import pointset : truncatePrecision;
        return truncatePrecision(P).tvalue1();
    }
    auto total = new BigInt[P.lg_length + 1] ; // sum
    foreach (x; P)
    {
        auto current = P.lg_length.empty_product(); // prod
        foreach (e; x)
        {
            immutable nu = e.nu_star(P.lg_length);
            foreach_reverse (i; nu..current.length)
            {
                current[i] -= current[i - nu] << nu;
            }
        }
        total = total.plus(current);
    }
    foreach (i, x; P.lg_length.Q1(P.dimension).times(total)) // (b^-m) is unnecessary
    {
        if (i && x.toLong)
        {
            return P.lg_length + 1 - i;
        }
    }
    return 0; // all zero, t = m + 1 - (m+1) according to Algorithm 1.
}

/** Compute Q[m](z) for Algorithm 1.

Q[m](z) = (1 + (b-1)z + (b<sup>2</sup>-b)z<sup>2</sup> + ... + (b<sup>m</sup>-b<sup>m-1</sup>)z<sup>m</sup>)<sup>s</sup> mod z<sup>m+1</sup> = [1, 2<sup>0</sup>, 2<sup>1</sup>, ..., 2<sup>m</sup>]<sup>s</sup>
*/
BigInt[] Q1(immutable size_t lg_length, immutable size_t dimension)
{
    auto ret = new BigInt[lg_length + 1];
    ret[0] = 1;
    ret[1] = 1;
    foreach (i; 1..lg_length)
    {
        ret[i + 1] = ret[i] << 1;
    }
    return ret.power(dimension);
}

/** Compute t-value of a digital net.

* Algorithm:
* Algorithm 2 of MacWilliams.<ol>
*   <li>input P = (x[n] : 0 <= n < b<sup>m</sup>), a subset of [0..1)<sup>s</sup>. Remarks are the same as Algorithm 1.</li>
*   <li>compute the coefficients of z<sup>a</sup> for (s-1)(m+1) <= a < s(m+1) of the polynomial<ul>
*       <li>Q(z) = b<sup>(s-1)m</sup>sum[x in P]prod[i in s](z<sup>mu(x[i])</sup>-z<sup>m+1</sup>) - (1+(b-1)z+...+(b<sup>m</sup>-b<sup>m-1</sup>)-b<sup>m</sup>z<sup>m+1</sup>)<sup>s</sup></li>
*   </ul>where<ul>
*       <li>mu(x) = 1+floor(lg(b<sup>m</sup>x)) = m+1 + floor(lg x) = m+1 - ceil(-lg x) = m+1 - nu(x)</li>
*   </ul></li>
*   <li>return t = (1 - s)(m + 1) + Q.degree</li>
* </ol>

* Computation:
* This function computes the reciprocal, namely,<ul>
*   <li>sum[x in P]prod[i in s](1-y<sup>nu(x[i])</sup>)</li>
* </ul>modulo y<sup>m+1</sup>, instead of<ul>
*   <li>sum[x in P]prod[i in s](z<sup>mu(x[i])</sup>-z<sup>m+1</sup>)</li>
* </ul> and compare the result with Q2.

* Params:
* same as Algorithm 1.

* Conditions:
* same as Algorithm 1.
*/
ulong tvalue2(R)(R P)
{
    if (P.lg_length < P.precision)
    {
        import pointset : truncatePrecision;
        return truncatePrecision(P).tvalue2();
    }
    auto total = new BigInt[P.lg_length + 1];
    foreach (x; P)
    {
        auto current = P.lg_length.empty_product();
        foreach (e; x)
        {
            immutable nu = e.nu_star(P.lg_length);
            foreach_reverse (i; nu..current.length)
            {
                current[i] -= current[i - nu];
            }
        }
        total = total.plus(current);
    }
    auto shift = (P.dimension - 1) * P.lg_length;
    // BigInt mask = (BigInt(1) << shift) - 1; // としたい
    foreach (i, c; P.lg_length.Q2(P.dimension))
    {
        if (total[i] << shift != c)
        //if ((c & mask) || total[i] != c >> shift) // と書きたいが、BigInt が cast(bool) や &, |, ^ を持っていないのでできない。まあ大したパフォーマンスダウンにはならないだろう。
        {
            if (i == 0)
            {
                total[i].writeln(" != ", c);
            }
            assert (i); return P.lg_length + 1 - i;
        }
    }
    return 0;
}

/** Compute Q(z) for algorithm 2.

Definition:
Q(z) = (1 + (b-1)z + (b<sup>2</sup>-b)z<sup>2</sup> + ... + (b<sup>m</sup>-b<sup>m-1</sup>)z<sup>m</sup> - b<sup>m</sup>z<sup>m+1</sup>)<sup>s</sup> = (1 + z + 2z<sup>2</sup> + ... + 2<sup>m-1</sup>z<sup>m</sup> - 2<sup>m</sup>z<sup>m+1</sup>)<sup>s</sup>

Computation:
This function computes the reciprocal of Q, namely, R(y) := (-y<sup>m+1</sup>)<sup>s</sup>Q(1/y), modulo y<sup>m+1</sup>.

Returns:
[b<sup>m</sup>, -b<sup>m-1</sup>, -b<sup>m-2</sup>, ..., -b<sup>0</sup><del>, -1</del>]<sup>s</sup>
*/
BigInt[] Q2(immutable size_t lg_length, immutable size_t dimension)
{
    auto ret = new BigInt[lg_length + 1];
    ret[$ - 1] = -1;
    foreach_reverse (i; 1..lg_length)
    {
        ret[i] = ret[i + 1] << 1;
    }
    ret[0] = -(ret[1] << 1);
    return ret.power(dimension);
}

/** Compute ceil(-lg(x >> m))
*
* Params:
*     m = lg|P|
*     x = a component of a vector in P; 0 <= x < 2^m
* Returns:
*     m + 1 if x = 0
*     m - (x.lg.floor) otherwise
*/
size_t nu_star(ulong x, immutable size_t m)
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

/// ditto
/// binary search version. maybe faster.
size_t nu_star_bs(ulong x, immutable size_t m)
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
        assert (nu_star(x, 14) == nu_star_bs(x, 14));
        assert (nu_star(x, 15) == nu_star_bs(x, 15));
        assert (nu_star(x, 16) == nu_star_bs(x, 16));
        assert (nu_star(x, 17) == nu_star_bs(x, 17));
        assert (nu_star(x, 32) == nu_star_bs(x, 32));
        assert (nu_star(x, 33) == nu_star_bs(x, 33));
        assert (nu_star(x, 60) == nu_star_bs(x, 60));
        assert (nu_star(x, 63) == nu_star_bs(x, 63));
        assert (nu_star(x, 64) == nu_star_bs(x, 64));
    }
    foreach (m; 1..65)
    {
        foreach (x; 0..m)
        {
            assert (nu_star(1UL << x, m) == nu_star_bs(1UL << x, m));
            assert (nu_star((1UL << x) - 1, m) == nu_star_bs((1UL << x) - 1, m));
        }
    }
    "...OK.".writeln();
}

/// '1' as polynomial
BigInt[] empty_product(immutable size_t m)
{
    return [BigInt(1)] ~ new BigInt[m];
}

/// polynomial multiplication
BigInt[] times(BigInt[] f, BigInt[] g)
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

import std.algorithm : min;

/// polynomial squaring
BigInt[] square(BigInt[] f)
{
    auto ret = new BigInt[f.length];
    foreach (i, c; f)
    {
        if (i << 1 < ret.length)
        {
            ret[i << 1] = c * c;
        }
        foreach (j; 0..i.min(ret.length-i))
        {
            ret[i + j] += (c * f[j]) << 1;
        }
    }
    return ret;
}

/// polynomial power
BigInt[] power(BigInt[] f, ulong s)
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

/// polynomial addition
BigInt[] plus(BigInt[] f, BigInt[] g)
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
