import lib.smallsearch;
import lib.pointsettype;
import std.stdio : stderr, writeln;
import std.string : strip;
import std.array : split;
import std.conv : to;

void main(string[] args)
{
    args = args[1..$];
    if (args.length != 6)
        return stderr.writeln("smallenumerate s s m m n n");
    immutable
        min_dimR = args[0].to!size_t(),
        max_dimR = args[1].to!size_t(),
        min_dimF = args[2].to!size_t(),
        max_dimF = args[3].to!size_t(),
        min_prec = args[4].to!size_t(),
        max_prec = args[5].to!size_t();
    foreach (precision; min_prec .. max_prec)
        foreach (dimR; min_dimR .. max_dimR)
            foreach (dimF; min_dimF .. max_dimF)
                foreach (basis; ReducedBasis!ubyte(precision, dimR, dimF))
                    ShiftedBasisPoints!ubyte(basis, Precision(precision)).toString().writeln();
}
