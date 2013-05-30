module testfunction;
import std.stdio;
//import std;
import std.math : abs, sin, PI, floor;

template Hellekalek(T ...)
{
    double f(double[] x)
    in
    {
        assert (T.length == x.length);
    }
    body
    {
        double ret = 1;
        foreach (i, c; T)
            ret *= x[i] ^^ c - 1 / (1 + c);
        return ret;
    }
    enum I = 0.0;
}

template Sobol2001(T ...)
{
    double f(double[] x)
    in
    {
        assert (T.length == x.length);
    }
    body
    {
        double ret = 1;
        foreach (i, c; T)
            ret *= ((4 * x[i] - 2).abs() + c) / (1 + c);
        return ret;
    }
    enum I = 1.0;
}

template Sobol1994(int s)
{
    double f(double[] x)
    in
    {
        assert (x.length == s);
    }
    body
    {
        double ret = 1;
        foreach (i; 0..s)
            ret *= (2 * x[i] + i + 1) / (i + 2);
        return ret;
    }
    enum I = 1.0;
}

template Owen(int s)
{
    double f(double[] x)
    in
    {
        assert (x.length == s);
    }
    body
    {
        double ret = 1;
        foreach (i; 0..s)
            ret *= x[i] - 1/2;
        return ret;
    }
    enum I = 0.0;
}

template RooArnold1(int s)
{
    double f(double[] x)
    in
    {
        assert (x.length == s);
    }
    body
    {
        double ret = 0;
        foreach (i; 0..s)
            ret += (4 * x[i] - 2).abs();
        return ret / s;
    }
    enum I = 1.0;
}

template RooArnold2(int s)
{
    double f(double[] x)
    in
    {
        assert (x.length == s);
    }
    body
    {
        double ret = 1;
        foreach (i; 0..s)
            ret *= (4 * x[i] - 2).abs();
        return ret;
    }
    enum I = 1.0;
}

template RooArnold3(int s)
{
    double f(double[] x)
    in
    {
        assert (x.length == s);
    }
    body
    {
        double ret = 1;
        foreach (i; 0..s)
            ret *= (x[i] * PI).sin() * PI / 2;
        return ret;
    }
    enum I = 1.0;
}

/// test function proposed by Kimikazu Sato on https://twitter.com/hamukazu/status/331605169722781696
/// template argument T ... must be array of _int_s, not _double_s.
template Hamukazu(T ...)
{
    double f(double[] x)
    in
    {
        assert (x.length == T.length);
    }
    body
    {
        double ret = 1;
        foreach (i, c; T)
            ret *= (c * x[i] - floor(c * x[i])) * 2;
        return ret;
    }
    enum I = 1.0;
}

//debug = hamu;
debug (hamu) unittest
{
    foreach (i; 0..100)
    {
        foreach (j; 0..100)
            Hamukazu!(5, 7).f([i * 0.01 + 0.005, j * 0.01 + 0.005]).write(",");
        writeln();
    }
}

version (none) template Genz2(int s)
{
    template G(T ...)
    {
        template H(U ...)
        {
            string gen_prod()
            {
                string ret;
                foreach (i; 0..s)
                {
                    ret ~= "ret *= " ~ (T[i] ^^ -2).to!string() ~ " + (x[" ~ i.to!string() ~ "] - " ~ U[i].to!string() ~  ") ^^ 2;";
                }
                return ret;
            }
            static assert (T.length == U.length);
            double f(double[] x)
            in
            {
                assert (x.length == T.length);
            }
            body
            {
                double ret = 1;
                mixin (gen_prod());
                return ret;
            }
        }
    }
}
// usage: "Genz2: ".writeln(Genz2!4.G!(1, 1, 1, 1).H!(1,1,1,1).f([1.0, 1.0, 1.0, 1.0]));
// which is lexically illegal.

debug (verbose) unittest
{
    "test function values at (1.0, 1.0, 1.0, 1.0):".writeln();
    "Hellekalek: ".writeln(Hellekalek!(1, 1, 1, 1).f([1.0, 1.0, 1.0, 1.0]));
    "Sobol2001: ".writeln(Sobol2001!(1, 1, 1, 1).f([1.0, 1.0, 1.0, 1.0]));
    "Sobol1994: ".writeln(Sobol1994!4.f([1.0, 1.0, 1.0, 1.0]));
    "Owen: ".writeln(Owen!4.f([1.0, 1.0, 1.0, 1.0]));
    "RooArnold1: ".writeln(RooArnold1!4.f([1.0, 1.0, 1.0, 1.0]));
    "RooArnold2: ".writeln(RooArnold2!4.f([1.0, 1.0, 1.0, 1.0]));
    "RooArnold3: ".writeln(RooArnold3!4.f([1.0, 1.0, 1.0, 1.0]));
    "Hamukazu: ".writeln(Hamukazu!(2, 3, 4, 5).f([1.0, 1.0, 1.0, 1.0]));
    "testfunction: unittest passed!".writeln();
}
