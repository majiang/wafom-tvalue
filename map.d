module main;

import lib.wafom;
import ui.input;
import ui.output;
import std.traits : ReturnType;

void main()
{
    foreach (digitalnet; getDigitalNets!ubyte())
    {
        digitalnet.writeWithMany!(ReturnType!(getDigitalNet!ubyte), biwafom, bimswafom)();
    }
}
