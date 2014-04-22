import lib.newsobol;
import lib.pointsettype;
import lib.wafom;
import lib.scramble;
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
        auto B = randomize(sobolBasis!uint(Precision(32), DimensionR(args[0].to!size_t()), DimensionF2(dimensionF2)));
        auto P = ShiftedBasisPoints!uint(B, Precision(32));
        "%s,%.15f,%.15f".writefln(P, P.bipmswafom().lg(), P.bimsnrtwafom().lg());
        foreach_reverse (n; 0..P.precision)
            foreach (s; 0..P.dimensionR)
            {
                auto sP = ShiftedBasisPoints!uint(modify(B, s, n), Precision(32));
                "%s,%.15f,%.15f".writefln(sP, sP.bipmswafom().lg(), sP.bimsnrtwafom().lg());
                modify(B, s, n);
            }
    }
}
