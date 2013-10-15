module exp;

import lib.wafom;

// standard exponential for small x.
real exp_custom(real x, size_t L)
in
{
    assert (0 < x);
    assert (x <= 1);
}
body
{
    real ret = 1;
    foreach (i; 0..L)
    {
        ret *= x;
        ret /= L - i;
        ret += 1;
    }
    return ret;
}

real expm1_custom(real x, size_t L)
in
{
    assert (0 < x);
    assert (x <= 1);
}
body
{
    real ret = 0;
    foreach (i; 0..L)
    {
        ret += 1;
        ret *= x;
        ret /= L - i;
    }
    return ret;
}

auto expm1div_custom(real x, size_t L)
in
{
    assert (0 < x);
    assert (x <= 1);
}
body
{
    real ret = 1;
    foreach (i; 0..L)
    {
        ret *= x;
        ret /= L+1-i;
        ret += 1;
    }
    return ret;
}

version (none) void main()
{
    import std.stdio;
debug
{
    foreach (i; 2..25)
        "( e )[%d] %.20f %.20a".writefln(i, 1.exp_custom(i), 1.exp_custom(i));
    foreach (i; 2..25)
        "(e-1)[%d] %.20f %.20a".writefln(i, 1.expm1_custom(i), 1.expm1_custom(i));
    foreach (i; 2..25)
        "exp' [%d] %.20f %.20a".writefln(i, 1.expm1div_custom(i), 1.expm1div_custom(i));
    auto x = 1.0;
    foreach (i; 1..100)
    {
        x *= 0.5;
        "%.20a ".writef(x);
        "%.20a".writefln(x.expm1div_custom(21));
    }
}
    //"L = %.20f = %.20a".writefln(get_L(), get_L());
    //"U = %.20f = %.20a".writefln(get_U(), get_U());
    import lib.pointset, lib.integration_error;
    import std.string;
    import genz;
    alias uint U;
    static immutable
        prec = U.sizeof << 3,
        dimR = 4,
        dimF = 12;
    auto gen_count = readln().strip().to!size_t();
    foreach (i; 0..gen_count)
    {
        //if (i % 100 == 0)
            stderr.writefln("%d / %d done.", i, gen_count);
retry:
        auto P = ShiftedBasisPoints!U(randomVectors!U(prec, dimR, dimF), prec);
        auto bip = P.bipwafom();
        if (2e-5 < bip)
            goto retry;
        auto xw = P.expwafoms();
        "\"%s\",%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f".writefln(
            P.toString(),
            bip,
            P.bipmswafom(),
            P.biwafom(),
            P.bimswafom(),
            P.integrationError!(exponential!dimR)(),
            xw[0],
            xw[1]
        );
    }
}

auto get_L()
{
    real ret = 1, x = 0x1p-64;
    foreach_reverse (i; 0..64)
    {
        ret *= x.expm1div_custom(21);
        x *= 2;
    }
    return ret;
}
auto get_U()
{
    return 1.expm1_custom(21);
}

static immutable U = get_U();
static immutable L = get_L();
double expwafomLeft(R)(R P)
{
    immutable s = P.dimensionR;
    return 0.5 * (
        P.bipwafom() * (L ^^ s - U ^^ s) +
        P.altwafom() * (U ^^ s + L ^^ s));
}
double expwafomRight(R)(R P)
{
    immutable s = P.dimensionR;
    return 0.5 * (
        P.bipwafom() * (U ^^ s - L^^ s) +
        P.altwafom() * (U ^^ s + L ^^ s));
}

import std.typecons : Tuple, tuple;
auto expwafoms(R)(R P)
{
    immutable u = U ^^ P.dimensionR;
    immutable l = L ^^ P.dimensionR;
    auto zero = P.bipwafom() * (u - l) * 0.5;
    auto one = P.altwafom() * (u + l) * 0.5;
    return tuple(-zero + one, zero + one);
}

double altwafom(R)(R P)
{
    return (P + one!(P.ComponentType)(P.precision, P.dimensionR)).bipwafom();
}

private U[] one(U)(size_t precision, size_t dimensionR)
{
    auto ret = new U[dimensionR];
    if (precision == U.sizeof << 3)
        foreach (ref x; ret)
            x = ~x;
    else
    {
        U c = 1;
        c <<= precision;
        c -= 1;
        foreach (ref x; ret)
            x = c;
    }
    return ret;
}
