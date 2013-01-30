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
    "factors(64) = ".writeln();
    foreach (x; factors(64))
    {
        if (x[0] == 1.0 && x[1] == 1.0)
        {
            // pass
        }
        else
        {
            x.writeln();
        }
    }
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

Remarks:
Using double, precision > 54 means factor = 1.
*/
double wafom(R)(R P)
{
    double ret = 0;
    //auto f = factors(P.precision);
    foreach (B; P)
    {
        double cur = 1;
        debug (speedup) double cur_backup = 1;
        foreach (l; B)
        {
            cur *= l.wafom_factor(P.precision);
            debug (speedup) foreach (j, c; f)
            {
                cur_backup *= c[(l >> j) & 1];
            }
        }
        debug (speedup) auto diff = cur - cur_backup;
        debug (speedup) assert (diff * diff < 1e-10);
        ret += cur - 1;
    }
//    debug "wafom is returning %f.".writefln(ret / P.length);
    return ret / P.length;
}

double wafom_factor(ulong x, ptrdiff_t precision)
{
    debug {
        auto memo = memoize!get_memo();
    } else {
        static memo = get_memo();
    }
    double ret = 1;
    while (0 < precision)
    {
        ret *= memo[precision - 1][x & 255];
        precision -= 8;
        x >>= 8;
    }
    return ret;
}


unittest
{
    /// speedup for foreach (j, c; f) cur *= c[(l >> j) & 1].
    void test_speedup(int t_max)
    {
        import std.datetime : StopWatch;
        double[256] memo_first; double[256] memo_second;
        auto f = factors(16);
        foreach (i; 0..256)
        {
            memo_first[i] = 1;
            foreach (j, c; f[8..$])
            {
                memo_first[i] *= c[(i >> j) & 1];
            }
            memo_second[i] = 1;
            foreach (j, c; f[0..8])
            {
                memo_second[i] *= c[(i >> j) & 1];
            }
        }
        "memo_first = ".writeln;
        memo_first.writeln();
        "memo_second = ".writeln;
        memo_second.writeln();
        foreach (i; 0..65536)
        {
            double cur = 1;
            foreach (j, c; f)
            {
                cur *= c[(i >> j) & 1];
            }
            auto diff = (memo_first[i >> 8] * memo_second[i & 255] - cur);
            assert (diff * diff < 1e-15);
        }
        StopWatch sw;
        sw.start();
        double total = 0;
        foreach (t; 0..t_max)
            foreach (i; 0..65536)
            {
                auto cur = 1;
                foreach (j, c; f)
                {
                    cur *= c[(i >> j) & 1];
                }
                //total += cur;
            }
        sw.stop();
        sw.peek().msecs.writeln(" ms by naive");
        total = 0;
        StopWatch swn;
        swn.start();
        foreach (t; 0..t_max)
        {
            foreach (i; 0..65536)
            {
                auto cur = memo_first[i >> 8] * memo_second[i & 255];
                //total += cur;
            }
        }
        swn.stop();
        swn.peek().msecs.writeln(" ms by memo-256");
    }
    version (speedup)
    {
        "test_speedup returns:".writeln();
        test_speedup(1000);
    }
    else test_speedup(0);
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

//debug = debug_memo;
debug (debug_memo) unittest
{
    "get_memo returns:".writeln();
    foreach (line; get_memo())
    {
        foreach (x; line)
        {
            x.write(", ");
        }
        writeln();
    }
    ".".writeln();
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