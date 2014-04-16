import lib.newsobol;
import lib.pointsettype;
import lib.wafom;
import std.math : lg  = log2;

void main(string[] args)
{
    import std.stdio;
    import std.array;

    args.popFront();
    if (args.length < 4)
        return "usage: w2_infty_one dimR dimBmin dimBmax count".writeln();
    foreach (dimensionF2; args[1].to!size_t() .. args[2].to!size_t()+1)
    {
        auto P = sobolPointSet!uint(Precision(32), DimensionR(args[0].to!size_t()), DimensionF2(dimensionF2));
        foreach (i; 0..args[3].to!size_t())
        {
            auto sP = P.scrambleRandomly();
            auto Q = sP[0], s = sP[1];
            "%s,%(%b%),%.15f,%.15f".writefln(Q, s, Q.bipmswafom().lg(), Q.bimsnrtwafom().lg());
        }
    }
}
