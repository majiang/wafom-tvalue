import std.algorithm;
import std.math;
import std.mathspecial : erf;
import std.range : iota;
alias reduce!((a, b) => a + b) sumation;
alias reduce!((a, b) => a * b) product;

enum size_t S = 4;

//alias monomial!(1, S) x1; alias monomial!(2, S) x2; alias monomial!(3, S) x3; alias monomial!(4, S) x4;
alias monomial!(5, S) x5; alias monomial!(6, S) x6; alias monomial!(7, S) x7; alias monomial!(8, S) x8;
alias exponential!(.25, S) exp025; alias exponential!(0.5, S) exp050; alias exponential!(1.0, S) exp100; alias exponential!(2.0, S) exp200; alias exponential!(4.0, S) exp400;
//alias cosineProd!(.25, S) cosp025; alias cosineProd!(0.5, S) cosp050; alias cosineProd!(1.0, S) cosp100; alias cosineProd!(2.0, S) cosp200; alias cosineProd!(4.0, S) cosp400;
alias cosineOfSum!(.25, S) coss025; alias cosineOfSum!(0.5, S) coss050; alias cosineOfSum!(1.0, S) coss100; alias cosineOfSum!(2.0, S) coss200; alias cosineOfSum!(4.0, S) coss400;
alias gaussian!(.25, S) gauss025; alias gaussian!(0.5, S) gauss050; alias gaussian!(1.0, S) gauss100; alias gaussian!(2.0, S) gauss200; alias gaussian!(4.0, S) gauss400;
//alias characteristic!(1.0 / 2.0, S) char1o2; alias characteristic!(1.0 / 3.0, S) char1o3; alias characteristic!(.25, S) char1o4; alias characteristic!(.20, S) char1o5; alias characteristic!(1.0 / 6.0, S) char1o6;
//alias sharppp!(1.0 / 2.0, S) sppp1o2; alias sharppp!(1.0 / 3.0, S) sppp1o3; alias sharppp!(.25, S) sppp1o4; alias sharppp!(.20, S) sppp1o5; alias sharppp!(1.0 / 6.0, S) sppp1o6;
alias productPeak!(.25, S) pp025; alias productPeak!(0.5, S) pp050; alias productPeak!(1.0, S) pp100; alias productPeak!(2.0, S) pp200; alias productPeak!(4.0, S) pp400;
alias disconti!S disco;
alias continuous!S conti;

real monomial_int(in size_t u, in size_t s)
{
    assert (s);
    if (s == 1)
        return 1.0 / (u + 1);
    return (u + 1).iota().map!(i => monomial_int(u - i, s - 1) * u.combin(i) / (i + 1))().sumation();
}

unittest
{
    import std.stdio;
    foreach (u; 1..13)
        foreach (s; 1..13)
            "%d %d -> %.15e".writefln(u, s, u.monomial_int(s));
}

template monomial(size_t u, size_t s)
{
    real f(in real[] xs)
    in
    {
        assert (xs.length == s);
    }
    body
    {
        return xs.sumation() ^^ u;
    }
    enum I = monomial_int(u, s);
}

// a = (large weight; good function) 0.25, 0.5, 1.0, 2.0, 4.0 (small weight; bad function)
template exponential(real a, size_t s)
{
    real f(in real[] xs)
    in
    {
        assert (xs.length == s);
    }
    body
    {
        return (xs.sumation() * a).exp();
    }
    auto I()@property{return (expm1(a) / a) ^^ s;}
}

template cosineProd(real a, size_t s)
{
    real f(in real[] xs)
    in
    {
        assert (xs.length == s);
    }
    body
    {
        return xs.map!(x => (a * x).cos())().product();
    }
    auto I()@property{return (a.sin() / a) ^^ s;}
}

private real combin(size_t n, size_t r)
in
{
    assert (r <= n);
}
body
{
    if (r * 2 > n)
        return n.combin(n - r);
    if (r == 0)
        return 1;
    return r.iota().map!(i => (n - i) / (i + 1.0L)).product();
}

template cosineOfSum(real a, size_t s)
{
    real f(in real[] xs)
    in
    {
        assert (xs.length == s);
    }
    body
    {
        return (xs.sumation() * a).cos();
    }
    auto I()@property
    {
        return (s + 1).iota().map!(i => ((s & 3) * -PI_2 + a * i).cos() * (-1.0) ^^ (s - i) * combin(s, i))().sumation() / a ^^ s;
    }
}

template gaussian(real a, size_t s)
{
    real f(in real[] xs)
    in
    {
        assert (xs.length == s);
    }
    body
    {
        return (xs.map!(x => -x * x).sumation() * a).exp();
    }
    auto I()@property{return (a.sqrt().erf() / (a.sqrt() * M_2_SQRTPI)) ^^ s;}
}

template continuous(size_t s)
{
    real f(in real[] xs)
    in
    {
        assert (xs.length == s);
    }
    body
    {
        return xs.map!(x => (3 * x).min((3 * x - 2).abs()))().product();
    }
    auto I()@property{return 0.5 ^^ s;}
}

template disconti(size_t s)
{
    real f(in real[] xs)
    in
    {
        assert (xs.length == s);
    }
    body
    {
        return xs.map!((real x) => (((1 <= 3 * x) && (3 * x < 2)) ? -1.0 : 1.0))().product();
    }
    auto I()@property{return (1.0/3.0) ^^ s;}
}

// a = 1/3 and binary rational
template characteristic(real a, size_t s)
{
    real f(in real[] xs)
    in
    {
        assert (xs.length == s);
    }
    body
    {
        return xs.map!(x => x < a ? -1 : 1).product();
    }
    auto I()@property{return (1 - 2 * a) ^^ s;}
}

// a = 1/3 and binary rational
template sharppp(real a, size_t s)
{
    real f(in real[] xs)
    {
        return xs.map!(x => (x - a).abs()).product();
    }
    auto I()@property{return (0.5 + a * a - a) ^^ s;}
}

// genz' product peak with u = 0 and a
template productPeak(real a, size_t s)
{
    real f(in real[] xs)
    {
        return 1 / xs.map!(x => a * x * x + 1).product();
    }
    auto I()@property{return (a.sqrt().atan() / a.sqrt()) ^^ s;}
}
