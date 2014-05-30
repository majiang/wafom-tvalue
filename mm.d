alias uint U;
import std.stdio;
import lib.pointsettype;
import lib.integral;
import std.string;
import std.math;
alias log2 lg;

auto f(real[] x)
{
    import std.algorithm;
    return (x.reduce!((a, b) => a + b)() * -2).exp();
}

auto I(size_t dimension)
{
    immutable I1 = (-2).expm1() / -2;
    real s = dimension;
    return I1 ^^ s;
}

void main()
{
    foreach (line; stdin.byLine())
    {
        auto P = line.fromString!U();
        "%s,%.15f".writefln(line.strip(), (P.bintegral!f() - P.dimensionR.I()).abs().lg());
    }
}
