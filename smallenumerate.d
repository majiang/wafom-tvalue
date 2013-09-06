import lib.smallsearch;
import lib.pointset;
import std.stdio : writeln, readln;
import std.string : strip;
import std.array : split;

void main()
{
    auto buf = readln().strip().split();
    immutable
        min_dimR = buf[0].to!size_t(),
        max_dimR = buf[1].to!size_t(),
        min_dimF = buf[2].to!size_t(),
        max_dimF = buf[3].to!size_t(),
        min_prec = buf[4].to!size_t(),
        max_prec = buf[5].to!size_t();
    foreach (precision; min_prec .. max_prec)
        foreach (dimR; min_dimR .. max_dimR)
            foreach (dimF; min_dimF .. max_dimF)
                foreach (basis; ReducedBasis!ubyte(precision, dimR, dimF))
                    ShiftedBasisPoints!ubyte(basis, precision).toString().writeln();
}