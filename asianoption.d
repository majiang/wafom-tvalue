module func.asianoption;

import std.math;
import std.typecons : Tuple;
alias Tuple!(double, "x", double, "y") Double;

/** risk-neutral rate = 3%/year
maturity = 1/3 year = 4 month
average = 4 times discrete monitoring (monthly)
initial price = $40
voratility = 20%
strike price = 35
*/
alias ctintegrand!(0.0295588022, 1.0/3.0, 4, 40, 0.2, 35) default_integrand;

/** Asian option integrand: ordinary function version.

sum[1 <= i <= S] exp (ai+by[i])
= exp(a+by[1]) (1+exp(a+by[2]) (1+exp() ... ))

to calculate this efficiently we can write
v <- exp(a+by[S])
foreach i in [1..S-1]:
    v <- v+1
    v <- exp(a+by[i])
*/
double integrand(in double r, in double T, in int S, in double P0, in double sigma, in double K, in double[] x)
in
{
    assert (S == x.length);
}
body
{
    immutable t = T / S;
    immutable stdv = sigma * sqrt(t);
    const y = gaussinv(S >> 1, x);
    immutable A = exp((r - sigma * sigma / 2) * t);
    double v = A * exp(stdv * y[0]);
    foreach_reverse (i; 1..S)
    {
        v += 1;
        v *= A * exp(stdv * y[i]);
    }
    v *= P0 / S;
    v -= K;
    if (v > 0) return v * exp(-r * T);
    return 0;
}


// compile-time exponential
real ctexp(real x)
{
    real ret = 0.0;
    real current = 1.0;
    int i = 1;
    while (true)
    {
        auto oret = ret;
        ret += current;
        if (ret == oret) return ret;
        current *= x / i;
        i += 1;
    }
}

version = speedup;
/** Asian option integrand: template version.

Remarks:
Not so fast though many constants are CTFEed.
*/
template ctintegrand(double r, double T, int S, double P0, double sigma, double K)
{
    static immutable t = T / S;
    static immutable stdv = sigma * sqrt(t);
    static immutable halfS = S >> 1;
    static immutable A = ctexp((r - sigma * sigma / 2) * t);
    static immutable factor = P0 / S;
    static immutable present_value = ctexp(-r * T);
    version (speedup) double ctintegrand(double[] x)
    {
        const y = gaussinv(halfS, x);
        double v = A * exp(stdv * y[0]);
        static if (S == 4) // loop unrolling
        {
            v += 1;
            v *= A * exp(stdv * y[1]);
            v += 1;
            v *= A * exp(stdv * y[2]);
            v += 1;
            v *= A * exp(stdv * y[3]);
        }
        else foreach_reverse (i; 1..S)
        {
            v += 1;
            v *= A * exp(stdv * y[i]);
        }
        v *= factor;
        v -= K;
        if (v > 0) return v * present_value;
        return 0;
    }
    else double ctintegrand(double[] x)
    {
        return integrand(r, T, S, P0, sigma, K, x);
    }
}

double[] gaussinv(in int n, in double[] x)
in
{
    assert(n << 1 == x.length);
}
body
{
    double[] ret;
    ret.length = n << 1;
    foreach (i; 0..n)
    {
        auto y = gaussinv_(Double(x[i << 1], x[i << 1 | 1]));
        ret[i << 1] = y.x;
        ret[i << 1 | 1] = y.y;
    }
    return ret;
}

double[n << 1] ctgaussinv(int n)(double[] x)
in
{
    assert(n << 1 == x.length);
}
body
{
    double[n << 1] ret;
    foreach (i; 0..n)
    {
        auto y = gaussinv_(Double(x[i << 1], x[i << 1 | 1]));
        ret[i << 1] = y.x;
        ret[i << 1 | 1] = y.y;
    }
    return ret;
}

Double gaussinv_(Double x)// [0..1)^2 -> gauss
{
    auto r = sqrt(-2 * log(1.0 - x.x));
    auto theta = 2 * PI * x.y;
    return Double(r * cos(theta), r * sin(theta));
}

version (none) unittest // not appropriate to place here
{
    import std.stdio;
    import lib.pointset : randomPoints;
    import lib.integral : bintegral;
    import std.traits : ReturnType;
    auto price = bintegral!(default_integrand, ReturnType!(randomPoints!uint))(randomPoints!uint(4, 20, 20));
    debug "asian option price = ".writeln(price);
}
