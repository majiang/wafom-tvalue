import lib.wafom;
import lib.pointset : nonshiftedRandomBasisPoints;
import lib.integration_error : shifts, shifteds, integrationErrors, squareRootMeanSquare;
import ui.input : getDigitalNets;
import tf;

enum precision = 32;
alias uint U;

void WRITEERRORS(R, F...)(R P)
{
    import std.stdio;
    P.toString().write();
    foreach (f; F)
        writef(",%.15f", P
            .shifteds(precision.shifts!(typeof (P))(4, 10000))
            .integrationErrors!f()
            .squareRootMeanSquare()
            );
}

void main()
{
    foreach (P; getDigitalNets!U())
    {
        P.WRITEERRORS!(typeof (P),
        x1, x2, x3, x4, x5, x6, x7, x8,
        exp025, exp050, exp100, exp200, exp400,
        cosp025, cosp050, cosp100, cosp200, cosp400,
        coss025, coss050, coss100, coss200, coss400,
        gauss025, gauss050, gauss100, gauss200, gauss400,
        char1o2, char1o3, char1o4, char1o5, char1o6,
        sppp1o2, sppp1o3, sppp1o4, sppp1o5, sppp1o6,
        pp025, pp050, pp100, pp200, pp400
    );
    import std.stdio;
    writeln();
    }
}
