import lib.wafom;
import ui.input;
import std.stdio;
import std.math : lg = log2;

void main()
{
    foreach (P; getDigitalNets!uint())
        "%s,%.15f,%.15f%s".writefln(P, P.bipwafom().lg(), P.bipmswafom().lg(), P.yoshiki_wep());
}
