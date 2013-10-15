module main;

import std.stdio;
import ui.input, ui.output;
import lib.pointset, lib.integration_error;
import genz;
import std.math : log2;

alias uint U;

void main()
{
    immutable numShift = 10000;
    U[][] shift;
    bool first = true;
    try while (true)
    {
        auto P = getDigitalNet!U();
        if (first)
        {
            first = false;
            shift = shifts!(typeof(P))(P.precision, P.dimensionR, numShift);
        }
        "%s,%.15f".writefln(P.toString, P.shifteds!(typeof(P))(shift).integrationErrors!(exponential!4).squareRootMeanSquare().log2());
    }
    catch
    {
    }
}
