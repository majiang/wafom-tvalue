module main;

import lib.wafom;
import ui.input;
import ui.output;
import std.traits : ReturnType;

void main()
{
    foreach (digitalnet; getDigitalNets!uint())
    {
        digitalnet.writeWith!(biwafom!(ReturnType!(getDigitalNet!uint)))();
        digitalnet.writeWithMany!(ReturnType!(getDigitalNet!uint), 
                                  biwafom!(ReturnType!(getDigitalNet!uint)),
                                  bimswafom!(ReturnType!(getDigitalNet!uint)),)();
    }
}
