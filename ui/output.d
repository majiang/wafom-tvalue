module ui.output;

import lib.pointset : ShiftedBasisPoints;

import std.stdio : writefln, writeln, writef, write;
import std.traits : isUnsigned, isFloatingPoint;

/// for type R with toString(), S and function f: R -> S, write x.toString and f(x) in one line.
void writeWith(alias f, R)(R x)
{
    x.toString().write(",");
    auto y = f(x);
    static if (is (typeof(y) : double))
        "%.15f".writefln(f(x));
    else 
        "%s".writefln(f(x));
}

void writeWithMany(R, func...)(R x)
{
    x.toString().write();
    foreach (f; func)
    {
        auto y = f(x);
        static if (is (typeof(y) : double))
            ",%.15f".writef(f(x));
        else
            ",%s".writef(f(x));
    }
    writeln();
}

void writePoints(U, R)(R P) if (isUnsigned!U)
{
    P.reset();
    foreach (X; P)
    {
        foreach (x; X)
            x.write(",");
        writeln();
    }
    writeln();
}

void writePoints(F, R)(R P) if (isFloatingPoint!F)
{
    F scale = 1;
    foreach (i; 0..P.precision)
        scale *= 0.5;
    F centering = scale * 0.5;
    
    P.reset();
    foreach (X; P)
    {
        foreach (x; X)
            "%.15f,".writef(x * scale + centering);
        writeln();
    }
    writeln();
}
