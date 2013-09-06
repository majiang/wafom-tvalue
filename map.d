module main;

import lib.wafom;
import ui.input;
import ui.output;
import std.traits : ReturnType;

void main()
{
    foreach (digitalnet; getDigitalNets!uint())
    {
        digitalnet.writeWithMany!(ReturnType!(getDigitalNet!uint), biwafom, bimswafom, binrtwafom, bimsnrtwafom)();
    }
}
