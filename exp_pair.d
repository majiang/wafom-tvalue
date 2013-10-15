module main;

import std.stdio, std.array, std.string, std.conv;
import std.container : heapify;
import lib.pointset : randomVectors, ShiftedBasisPoints;
import lib.wafom : bipwafom, bipmswafom;
import std.algorithm : min;
import std.math : log2;
import exp : expwafoms;

alias uint U;
alias ShiftedBasisPoints!U BP;
import std.typecons : Tuple;
alias Tuple!(double, double, U[][]) DN;
alias Tuple!(double, size_t, size_t) Pair;

void main()
{
    stderr.writeln("precision(<=32) dimF2 dimR generate low-wafom best");
    auto buf = readln().strip().split();
    immutable precision = buf[0].to!size_t().min(32);
    immutable dimensionF2 = buf[1].to!size_t();
    immutable dimensionR = buf[2].to!size_t();
    immutable count = buf[3].to!size_t();
    immutable low = buf[4].to!size_t();
    immutable best = buf[5].to!size_t();

foreach (tmp; 0..best){

    auto PH = (new DN[low]).heapify(0);
    foreach (i; 0..count)
    {
        auto B = randomVectors!uint(precision, dimensionR, dimensionF2);
        auto P = BP(B, precision);
        PH.conditionalInsert(DN(P.bipwafom().log2(), P.bipmswafom().log2(), B));
    }
    auto PS = PH.release();

    auto H = (new Pair[1]).heapify(0);

    foreach (i, P; PS)
        foreach (j, Q; PS)
        {
            if (i == j)
                break;
            auto p = (P[0] - Q[0]) * (P[1] - Q[1]);
            H.conditionalInsert(Pair(p, i, j));
        }
    foreach (pair; H.release())
    {
        auto i = pair[1];
        auto j = pair[2];
        auto P = (PS[i][0] < PS[j][0]) ? PS[i] : PS[j];
        auto Q = (PS[i][0] < PS[j][0]) ? PS[j] : PS[i];
        "%s,%.15f,%.15f".writefln(BP(P[2], precision).toString(), P[0], P[1]);
        "%s,%.15f,%.15f".writefln(BP(Q[2], precision).toString(), Q[0], Q[1]);
    }

}

}
