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
            .shifteds(P.precision.shifts!(typeof (P))(P.dimensionR, 8192))
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
            x6,
            exponential!(2.0/3.0, S),
            exponential!(3.0/2.0, S),
            coss100,
            gauss100,
            pp100,
            conti,
            disco
        );
    }
}
