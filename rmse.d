import lib.wafom;
import lib.pointset : nonshiftedRandomBasisPoints;
import lib.integration_error : shifts, shifteds, integrationErrors, squareRootMeanSquare;
import ui.input : getDigitalNets;
import tf;

alias uint U;

void WRITEERRORS(R, F...)(R P)
{
    import std.stdio;
    "%s,%.15e".writef(P, P.bipmswafom()); // unnecessary recalculate
    foreach (f; F)
        writef(",%.15e", P
            .shifteds(P.precision.shifts!(typeof (P))(P.dimensionR, 10000))
            .integrationErrors!f()
            .squareRootMeanSquare()
            );
    writeln();
}

void main()
{
    foreach (P; getDigitalNets!U())
    {
        P.WRITEERRORS!(typeof (P),
            x3, x4, x5, x6, x7,
            exp100, cosp100, coss100, gauss100, pp100
            );
    }
}
