module genz;

import std.stdio;
import std.math;

/// simple exponential function: integral is (e-1)^s.
template exponential(size_t s)
{
    double f(double[] x)
    in
    {
        assert (x.length == s);
    }
    body
    {
        auto e = 0.0;
        foreach (c; x)
            e += c;
        return exp(e);
    }
    enum I = (E - 1) ^^ s;
}

template negexp(size_t s)
{
    double f(double[] x)
    in
    {
        assert (x.length == s);
    }
    body
    {
        auto e = 0.0;
        foreach (c; x)
            e += c;
        return exp(s - e);
    }
    enum I = (E - 1) ^^ s;
}


version (none)
{
double oscillatory(double[] x, double u, double[] a)
{
    auto t = 2 * PI * u;
    foreach (i, c; x)
        t += a[i] * c;
    return t.cos();
}

double prodpeak(double[] x, double[] u, double[] a)
{
    auto ret = 1.0;
    foreach (i, c; x)
        ret /= 1 / (a[i] * a[i]) + (c - u[i]) * (c - u[i]);
    return ret;
}

double cornpeak(double[] x, double[] a)
{
    auto ret = 1.0;
    foreach (i, c; x)
        ret += c + a[i];
    return ret ^^ -(a.length + 1);
}

double gaussian(double[] x, double[] a, double[] u)
{
    auto exponent = 0.0;
    foreach (i, c; x)
        exponent += a[i] * a[i] * (c - u[i]) * (c - u[i]);
    return exp(-exponent);
}

double continuous(double[] x, double[] a, double[] u)
{
    auto exponent = 0.0;
    foreach (i, c; x)
        exponent += a[i] * abs(c - u[i]);
    return exp(-exponent);
}

double disconti(double[] x, double[] a, double u, double v)
{
    if (x[0] > u || x[1] > v)
        return 0;
    auto exponent = 0.0;
    foreach (i, c; x)
        exponent += a[i] * c;
    return exp(exponent);
}

}