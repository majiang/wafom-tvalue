module ui.input;

import std.stdio : stdin;
import std.array : split;
import std.conv : to;
import std.string : strip;
import std.traits : isUnsigned;

auto getDigitalNet(U)() if (isUnsigned!U)
{
    import lib.pointset : fromString;
    return stdin.readln().fromString!U();
}

auto getDigitalNets(U)() if (isUnsigned!U)
{
    import lib.pointset : fromString;
    import std.traits : ReturnType;
    alias ReturnType!(fromString!U) P;
    struct R
    {
        int opApply(int delegate(P) dg)
        {
            foreach (line; stdin.byLine)
                if (auto result = dg(line.fromString!U()))
                    return result;
            return 0;
        }
    }
    return R();
}
