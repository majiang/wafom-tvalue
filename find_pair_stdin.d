module main;

import std.stdio, std.array, std.string, std.conv;
import std.container : heapify;
import lib.pointset : randomVectors, ShiftedBasisPoints;
import lib.wafom : bipwafom, bipmswafom;
import std.algorithm : min;
import std.math : lg = log2;
import ui.input : getDigitalNets;
alias uint U;
alias getDigitalNets!U getDN;
alias ShiftedBasisPoints!U BP;
import std.typecons : Tuple;
alias Tuple!(real, real, U[][]) DN;
alias Tuple!(real, size_t, size_t) Pair;

U[][] ddup(in U[][] x)
{
    U[][] ret;
    foreach (l; x)
    {
        ret.length += 1;
        foreach (c; l)
            ret[$-1] ~= c;
    }
    return ret;
}

void main()
{
    DN[] PS;
    size_t precision;
    foreach (P; getDN())
    {
        precision = P.precision;
        PS ~= DN(P.bipwafom().lg(), P.bipmswafom().lg(), P.basis.ddup);
    }

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
