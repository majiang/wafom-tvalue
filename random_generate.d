module main;
import std.stdio, std.array, std.string, std.conv;
import lib.pointset : randomVector, nonshiftedRandomBasisPoints;
import ui.output : writePoints;

/** Short sample of the library.

<ol>
<li>Randomly generate a digital net P with precision, dimensionF2 and dimensionR specified.</li>
<li>Randomly generate a vector v with the same precision and dimensionR as P.</li>
<li>Output the points of P.</li>
<li>Output the points of P + &lt;v&gt;</li>
</ol>
*/

void main()
{
    auto buf = readln().strip().split();
    immutable precision = buf[0].to!size_t();
    immutable dimensionF2 = buf[1].to!size_t();
    immutable dimensionR = buf[2].to!size_t();
    auto P = nonshiftedRandomBasisPoints!ubyte(precision, dimensionR, dimensionF2);
    auto v = randomVector!ubyte(precision, dimensionR);
    P.writePoints!double();
    P.writePoints!ubyte();
}
