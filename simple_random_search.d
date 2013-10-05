module main;

/** Short sample of the library.

<ol>
<li>Randomly generate digital nets P[] with precision, dimensionF2 and dimensionR specified.</li>
<li>Output WAFOM-best point sets.</li>
</ol>
*/

import std.stdio, std.array, std.string, std.conv;
import std.container : heapify;
import lib.pointset : randomVectors, ShiftedBasisPoints;
import lib.wafom : biwafom, bimswafom;
import std.algorithm : min;
import ui.output : writePoints;
alias uint U;
alias ShiftedBasisPoints!U BP;
import std.typecons : Tuple;
alias Tuple!(double, U[][]) DN;

void main()
{
    stderr.writeln("precision(<=32) dimF2 dimR generate choose");
    auto buf = readln().strip().split();
    immutable precision = buf[0].to!size_t().min(32);
    immutable dimensionF2 = buf[1].to!size_t();
    immutable dimensionR = buf[2].to!size_t();
    immutable count = buf[3].to!size_t();
    immutable best = buf[4].to!size_t();

    auto H = (new DN[best]).heapify(0);
    foreach (i; 0..count)
    {
        auto B = randomVectors!uint(precision, dimensionR, dimensionF2);
        H.conditionalInsert(DN(BP(B, precision).biwafom(), B));
    }
    foreach (dn; H.release())
    {
        auto P = BP(dn[1], precision);
        "%s,%.15f,%.15f".writefln(P.toString(), dn[0], P.bimswafom());
    }
}