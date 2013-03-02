import std.stdio;
import pointset : BasisPoints;

/** Perform integration of a function f: [0..1)<sup>s</sup> -> R by the point set P.

* s is implicitly given as P.dimension. f must support double opCall(double[] x)).

Usage:
----------------
immutable size_t dimension;
double f(double[] x)
{
    assert (x.length == dimension);
    //...
    return something;
}
double result = integral!f(randomPoints(dimension, precision, lg_length));
----------------
*/
double integral(alias f)(BasisPoints P)
{
    double factor = 1.0 / (1UL << P.lg_length);
    double shift = 0.5 / (1UL << P.precision);
    double result = 0;
    foreach (x; P)
    {
        result += f(x.affine(factor, shift));
    }
    return result / P.length;
}

double integral_ao(double r, double T, int S, double P0, double sigma, double K, BasisPoints P)
{
    import asianoption : integrand;
    double factor = 1.0 / (1UL << P.lg_length);
    double shift = 0.5 / (1UL << P.precision);
    double result = 0;
    foreach (x; P)
    {
        result += integrand(r, T, S, P0, sigma, K,
            x.affine(factor, shift));
    }
    return result / P.length;
}

void integral_sobol(BasisPoints P)
{
    auto S = P.dimension;
    foreach (r; [0.02955880224, 0.03, 0.04879016416, 0.05, 0.06765864847, 0.07])
        foreach (T; [1.0/3.0])
//            foreach (S; [])
                foreach (P0; [40.0])
                    foreach (sigma; [0.2, 0.3, 0.4])
                        foreach (K; [35.0, 40.0, 45.0])
                            "Asian Option (r = %f, T = %f, S = %d, P0 = %f, sigma = %f, K = %f) . price = %f".writefln(
                                r, T, S, P0, sigma, K, integral_ao(r, T, S, P0, sigma, K, P.save()));
}

private double[] affine(ulong[] x, double factor, double shift)
{
    auto ret = new double[x.length];
    foreach (i, c; x)
    {
        ret[i] = c * factor + shift;
    }
    return ret;
}
