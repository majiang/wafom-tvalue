module main;

/** Short sample of the library.

<ol>
<li>Randomly generate digital nets P[] with precision, dimensionF2 and dimensionR specified.</li>
<li>Output WAFOM-best point sets.</li>
</ol>
*/

import std.stdio, std.array, std.string, std.conv;
import std.container : heapify;
import lib.pointsettype;
import lib.wafom : criterion = biwafom;
import std.algorithm : min, sort;
alias uint U;
alias ShiftedBasisPoints!U BP;
import std.typecons : Tuple;
import std.math : lg = log2;
alias Tuple!(double, U[][]) DN;

auto ddup(T)(in T[][] arr)
{
    T[][] ret;
    foreach (row; arr)
        ret ~= row.dup;
    return ret;
}

auto ifany(T, T defaultValue)(T[] arr, size_t idx)
{
    if (arr.length <= idx)
        return defaultValue;
    return arr[idx];
}

void main(string[] args)
{
    if (args.length < 5)
        return stderr.writeln("precision(<=32) dimF2 dimR generate choose");
    auto buf = args[1..$];
    immutable precision = buf[0].to!size_t().min(32);
    immutable dimensionF2 = buf[1].to!size_t();
    immutable dimensionR = buf[2].to!size_t();
    immutable count = buf[3].to!size_t();
    immutable best = buf.ifany!(string, "0")(4).to!size_t();


    if (best)
    {
        auto H = (new DN[best]).heapify(0);
        foreach (i; 0..count)
        {
            auto P = randomPointSet!U(Precision(precision), DimensionR(dimensionR), DimensionF2(dimensionF2));
            H.conditionalInsert(DN(P.criterion(), P.basis.ddup()));
        }
        foreach (dn; H.release().sort())
            "%s,%.15f".writefln(BP(dn[1], Precision(precision)).toString(), dn[0].lg());
    }
    else
    {
        import lib.wafom : bipwafom, biwafom, bimswafom, bipmswafom;
        foreach (i; 0..count)
        {
            auto P = randomPointSet!U(Precision(precision), DimensionR(dimensionR), DimensionF2(dimensionF2));
            "%s,%.15f,%.15f,%.15f,%.15f".writefln(P.toString(), P.bipwafom().lg(), P.bipmswafom().lg(), P.biwafom().lg(), P.bimswafom().lg());
        }
    }
}
