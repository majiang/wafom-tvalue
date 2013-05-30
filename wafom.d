module wafom;

import std.functional : memoize;

debug import std.stdio;


debug = speedup;
/** Compute wafom of a general quasi-Monte Carlo point set.

* Algorithm:
* Equation (4.2) of wafom-arxiv.<ul>
*   <li>input P = (x[n] : 0 <= n < b<sup>m</sup>) in [0..1)<sup>s</sup>.</li>
*   <li>return sum[x in P](prod[i in s][j in m](1+(-1)<sup>x[i,j]</sup>2<sup>-j-1</sup>)-1) / |P|.</li>
* </ul>

* Params:
* P = an array-of-integral-type range which has attributes precision and dimensionF2. If v is an output vector, v*2<sup>-precision</sup> is an element of [0..2<sup>m</sup>)<sup>s</sup>.

Remarks:
Using double, precision > 54 is no better than precision = 54.
*/
double wafom(R)(R P)
{
    assert (P.precision);
    double ret = 0;
    debug (speedup) auto f = factors(P.precision, 2);
    foreach (B; P)
    {
        double cur = 1;
        debug (speedup) double cur_backup = 1;
        foreach (l; B)
        {
            cur *= l.wafom_factor(P.precision);
            debug (speedup) foreach (j, c; f) cur_backup *= c[(l >> j) & 1];
        }
        debug (speedup) auto diff = cur - cur_backup;
        debug (speedup) assert (diff * diff < 1e-10);
        ret += cur;
    }
    foreach (i; 0..P.dimensionF2)
        ret *= 0.5;
    return ret - 1;
}

double mswafom(R)(R P)
{
    import std.math : sqrt;
    double ret = 0;
    auto f = factors(P.precision, 4);
    foreach (B; P)
    {
        double cur = 1;
        foreach (l; B)
        {
            foreach (j, c; f)
            {
                cur *= c[(l >> j) & 1];
            }
        }
        ret += cur;
    }
    return sqrt((ret / P.length) - 1);
}

unittest
{
    import pointset : randomPoints;
    foreach (i; 0..10)
    {
        auto P = randomPoints!uint(4, 32, 10);
        auto w = P.wafom();
        auto m = P.mswafom();
        debug (verbose) "wafom = %.15f; mswfm = %.15f".writefln(w, m);
    }
}


double wafom_factor(ulong x, ptrdiff_t precision)
{
    debug auto memo = memoize!get_memo(); // just memoize for debug, to compile faster.
    else static memo = get_memo(); // CTFE for release, to execute faster.
    double ret = 1;
    while (0 < precision)
    {
        ret *= memo[precision - 1][x & 255];
        precision -= 8;
        x >>= 8;
    }
    return ret;
}


double[256][64] get_memo()
{
    import std.algorithm : min, max;
    double[256][64] ret;
    auto f = _factors(64);
    foreach (i; 0..64)
    {
        foreach (j; 0..2 << min(i, 7))
        {
            ret[i][j] = 1;
            foreach (k, c; f[$-(i+1)..$-max(0, i-7)])
            {
                ret[i][j] *= c[(j >> k) & 1];
            }
        }
    }
    return ret;
}


double[2][] _factors(size_t precision, size_t base = 2)
{
    assert (precision);
    auto ret = new double[2][precision];
    double recip = 1.0 / base;
    foreach (i; 0..2)
    {
        ret[$-1][i] = recip; // Mr. Yoshiki suggests recip * recip
        foreach_reverse (j; 1..precision)
        {
            ret[j-1][i] = ret[j][i] * recip;
        }
    }
    foreach (j; 0..precision)
    {
        ret[j][0] += 1;
        ret[j][1] *= -1;
        ret[j][1] += 1;
    }
    return ret;
}
alias memoize!_factors factors;

unittest
{
    debug (verbose) "factors(64) = ".writeln();
    foreach (x; factors(64, 2))
        debug (verbose)
            if (!(x[0] == 1.0 && x[1] == 1.0))
                x.writeln();
}
