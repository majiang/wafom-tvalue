module wafom;

import std.functional : memoize;

debug import std.stdio;

double[2][] _factors(size_t precision)
{
    auto ret = new double[2][precision];
    foreach (i; 0..2)
    {
        ret[$-1][i] = 0.5;
        foreach_reverse (j; 1..precision)
        {
            ret[j-1][i] = ret[j][i] * 0.5;
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
    factors(10).writeln();
}

/** Compute wafom of a general quasi-Monte Carlo point set.

Algorithm:
Equation (4.2) of wafom-arxiv.
* 1: input P = (x[n] : 0 <= n < b^m) in [0..1)^s.
** b = 2
** m = P.precision
** s = P.dimenson
*
* 2: return sum[n in b^m](prod[t in s][j in m](1+(-1)^x[n,t,j]2^(-j-1))-1) / b^m.

Params:
P = a subset of  [0..2^m)^s
*/
double wafom(R)(R P)
{
    double ret = 0;
    auto f = factors(P.precision);
    //auto power_of_half = new double[P.precision];
    //power_of_half[$ - 1] = 0.5;
    //foreach_reverse (i; 1..P.precision)
    //{
    //    power_of_half[i - 1] = power_of_half[i] * 0.5;
    //}
    //power_of_half.writeln();
    foreach (B; P)
    {
        double cur = 1;
        foreach (l; B)
        {
            foreach (j, c; f)
            {
                cur *= c[(l >> j) & 1];
            }
            //foreach (j, c; power_of_half)
            //{
            //    if (l >> j & 1)
            //    {
            //        cur *= 1 - c;
            //    }
            //    else
            //    {
            //        cur *= 1 + c;
            //    }
            //}
        }
        ret += cur - 1;
    }
//    debug "wafom is returning %f.".writefln(ret / P.length);
    return ret / P.length;
}

version (digitalnet){
double wafom(ulong n, ulong b, ulong m, ulong s)(ulong[s][b ^^ m] P)
{
    double ret = 0;
    foreach (B; P)
    {
        double cur = 1;
        foreach (l; B)
        {
            foreach (j; 0..n)
            {
                if (l >> j & 1)
                {
                    cur *= 1 - 0.5 ^^ (j + 1);
                }
                else
                {
                    cur *= 1 + 0.5 ^^ (j + 1);
                }
            }
        }
        ret += cur - 1;
    }
    return ret / N;
}

double wafom(ulong n, ulong s, ulong N)(bool[n][s][N] P)
{
    double ret = 0;
    foreach (B; P)
    {
        double cur = 1;
        foreach (l; B)
        {
            foreach (j, e; l)
            {
                if (e)
                {
                    cur *= 1 - 0.5 ^^ (j + 1);
                }
                else
                {
                    cur *= 1 + 0.5 ^^ (j + 1);
                }
            }
        }
        ret += cur - 1;
    }
    return ret / N;
}
}