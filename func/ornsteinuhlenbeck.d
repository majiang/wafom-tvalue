module func.ornsteinuhlenbeck;

alias oup!(0.03, 5, 100, 100) integrandLowShort;
alias oup!(0.03, 5, 500, 100) integrandHighShort;
alias oup!(0.03, 5, 100, 2000) integrandLowLong;
alias oup!(0.03, 5, 500, 2000) integrandHighLong;

unittest
{
    import std.stdio;
    "ornstein uhlenbeck:".writeln();
    import lib.pointset : defaultSobols;
    import lib.integral : bintegral;
    foreach (m; [8, 10, 12, 14, 16, 18])
    //foreach (m; [4, 6, 8])
        foreach (s; [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384])
        //foreach (s; [2, 4, 8])
        {
            writefln("low-short,s=,%d,m=,%d,%.15f", s, m, bintegral!integrandLowShort(defaultSobols(s, 32, m)));
            writefln("low--long,s=,%d,m=,%d,%.15f", s, m, bintegral!integrandLowLong(defaultSobols(s, 32, m)));
            writefln("highshort,s=,%d,m=,%d,%.15f", s, m, bintegral!integrandHighShort(defaultSobols(s, 32, m)));
            writefln("high-long,s=,%d,m=,%d,%.15f", s, m, bintegral!integrandHighLong(defaultSobols(s, 32, m)));
        }
}

template oup(double theta, double sigma, double R0, double T)
{
    import std.math : sqrt;
    immutable Theta = T * theta;
    immutable s0 = R0 / sigma; // these values are evaluated on compile time
    double oup(in double[] x)
    {
        immutable n = x.length;
        immutable r = 1 - Theta / n;
        const w = x.boxmuller();
        immutable sw = (T / n).sqrt();
        double s = s0;
        foreach (c; w)
        {
            s *= r;
            s += c * sw;
        }
        return sigma * s;
    }
}

double[] boxmuller(in double[] x)
{
    import std.math;
    if (x.length & 1) assert (false);
    auto w = new double[x.length];
    foreach (p; 0..(x.length >> 1))
    {
        immutable i = p << 1;
        immutable j = i | 1;
        immutable r = sqrt(-2 * log(1.0 - x[i]));
        if (r.isNaN())
        {
            import std.stdio;
            "sqrt(-2 * log(1.0-%f)) = sqrt(-2 * %f) = sqrt(%f)".writefln(x[i], log(1.0-x[i]), -2*log(1.0-x[i]));
        }
        immutable theta = 2 * PI * x[j];
        w[i] = theta.cos() * r;
        w[j] = theta.sin() * r;
        assert (!w[i].isNaN());
        assert (!w[j].isNaN());
    }
    return w;
}

