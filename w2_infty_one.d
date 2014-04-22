import lib.newsobol;
import lib.pointsettype;
import lib.wafom;
import lib.scramble;
import std.math : lg  = log2;

void main(string[] args)
{
    import std.stdio;
    import std.array;
    import std.algorithm;

    args.popFront();
    if (args.length < 6)
        return "usage: w2_infty_one dimR dimBmin dimBmax count r...".writeln();
    foreach (dimensionF2; args[1].to!size_t() .. args[2].to!size_t()+1)
    {
        auto B = sobolBasis!uint(Precision(32), DimensionR(args[0].to!size_t()), DimensionF2(dimensionF2));
        auto P = ShiftedBasisPoints!uint(B, Precision(32));
        "%s,%.15f,%.15f".writefln(P, P.bipmswafom().lg(), P.bimsnrtwafom().lg());
        foreach (r; args[4..$].map!(to!size_t)())
            foreach (i; 0..args[3].to!size_t())
            {
                auto sP = ShiftedBasisPoints!uint(B.multiplyRandomMatrices(32, r), Precision(32));
                "%s,%.15f,%.15f,%d".writefln(sP, sP.bipmswafom().lg(), sP.bimsnrtwafom().lg(), r);
            }
    }
}
