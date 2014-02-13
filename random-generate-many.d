module main;
import std.stdio, std.array, std.string, std.conv;
import lib.pointset : nonshiftedRandomBasisPoints;

/** Short sample of the library.

<ol>
<li>Randomly generate c digital nets P[] with precision, dimensionF2 and dimensionR specified.</li>
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
    immutable count = buf[3].to!size_t();
    foreach (i; 0..count)
        nonshiftedRandomBasisPoints!uint(precision, dimensionR, dimensionF2).writeln();
}
