module wafom;

import sobol : Sobols;


/** Compute wafom of a digital net formed from Sobol sequences.

R must support iteration and have properties bits, length, dimension.

Params:
P = a subset of  [0..(1<<m))^s
*/
double wafom(R)(R P)//, size_t bits, size_t dimension)
{
    double ret = 0;
    double[] power_of_half;
    power_of_half.length = P.bits;
    power_of_half[0] = 0.5;
    foreach (i; 1..P.bits)
    {
        power_of_half[i] = power_of_half[i - 1] * 0.5;
    }
    foreach (B; P)
    {
        double cur = 1;
        foreach (l; B)
        {
            foreach (j, c; power_of_half)
            {
                if (l >> j & 1)
                {
                    cur *= 1 - c;
                }
                else
                {
                    cur *= 1 + c;
                }
            }
        }
        ret += cur - 1;
    }
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