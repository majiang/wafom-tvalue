import lib.sobol;
import lib.pointsettype;
import lib.wafom;

void main(string[] args)
{
    import std.stdio;
    import std.array;

    args.popFront();
    foreach (dimensionF2; args[1].to!size_t() .. args[2].to!size_t())
    {
        auto P = defaultSobols(DimensionR(args[0].to!size_t()), Precision(32), DimensionF2(dimensionF2));
        foreach (i; 0..args[3].to!size_t())
        {
            auto Q = P.scrambleRandomly();
            "%s,%.15f,%.15f".writefln(Q, Q.bimswafom().lg(), Q.bimsnrtwafom().lg());
        }
    }
}
