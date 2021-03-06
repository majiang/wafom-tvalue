module main;

import std.stdio;
import std.string : strip;
import std.array : split;
import std.conv : to;

alias uint U;

import lib.pointset : nonshiftedRandomBasisPoints;
import lib.wafom : bipwafom, bipmswafom;


version (prep) void main()
{
    stderr.writeln("prec dimR dimF prep rate count(unread here)");
    auto buf = readln().strip().split();
    immutable
        prec = buf[0].to!size_t(),
        dimR = buf[1].to!size_t(),
        dimF = buf[2].to!size_t(),
        prep = buf[3].to!size_t(),
        rate = buf[4].to!size_t(),
        count = buf[5].to!size_t();
    auto threshold = preparation!(U, bipwafom, bipmswafom)(prec, dimR, dimF, prep, rate);
    "%d %d %d %d %.20f %.20f".writefln(prec, dimR, dimF, count, threshold[0], threshold[1]);
}
version (exec) void main()
{
    stderr.writeln("prec dimR dimF count WAFOM RMSW [AND]");
    auto buf = readln().strip().split();
    immutable
        prec = buf[0].to!size_t(),
        dimR = buf[1].to!size_t(),
        dimF = buf[2].to!size_t(),
        count = buf[3].to!size_t(),
        threshold = Tuple!(double, double)(buf[4].to!double(), buf[5].to!double());
    if (6 < buf.length)
        foreach (P; generation!(U, bipwafom, true, bipmswafom)(prec, dimR, dimF, threshold, count))
            P.toString().writeln();
    else
        foreach (P; generation!(U, bipwafom, false, bipmswafom)(prec, dimR, dimF, threshold, count))
            P.toString().writeln();
}

version = silent;

version (verbose)
    auto milestone = [0.00, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.98, 0.99, 1.00];
version (medium)
    auto milestone = [0.00, 0.03, 0.09, 0.27, 0.81];

import std.typecons : Tuple;
auto generation(U, alias criterion0, bool AND, alias criterion1)
(
    in size_t prec, in size_t dimR, in size_t dimF,
    in Tuple!(double, double) threshold, in size_t count
)
{
    static struct R
    {
        size_t length;
        immutable size_t prec;
        immutable size_t dimR;
        immutable size_t dimF;
        immutable Tuple!(double, double) threshold;
        immutable size_t initLength;
        void popFront() {--length; version (silent) {} else{
            import std.stdio;
            import std.datetime : Clock;
            foreach (m; milestone)
                if (length - 1 < initLength * m && initLength * m <= length)
                    stderr.writefln("%.2f left at %s", m, Clock.currTime());
        }}
        @property bool empty() {return !length;}
        auto front()
        {   retry:
            auto P = nonshiftedRandomBasisPoints!U(prec, dimR, dimF);
            static if (AND)
            {
                if (threshold[0] < criterion0(P) || threshold[1] < criterion1(P))
                    goto retry;
            }
            else
            {
                if (threshold[0] < criterion0(P) && threshold[1] < criterion1(P))
                    goto retry;
            }
            return P;
        }
    }
    return R(count, prec, dimR, dimF, threshold, count);
}

auto preparation(U, alias criterion0, alias criterion1)
(
    in size_t prec, in size_t dimR, in size_t dimF,
    in size_t prep, in size_t rate
)
{
    import std.stdio : stderr;
    import std.datetime : Clock;
    stderr.writeln("start preparation at ", Clock.currTime());
    scope (success)
        stderr.writeln("finish preparation at ", Clock.currTime());
    scope (failure)
        stderr.writeln("terminate preparation by error");

    import std.container : heapify;
    auto pq0 = (new double[rate]).heapify(0);
    auto pq1 = (new double[rate]).heapify(0);
    foreach (i; 0..prep)
    {
        auto P = nonshiftedRandomBasisPoints!U(prec, dimR, dimF);
        pq0.conditionalInsert(criterion0(P));
        pq1.conditionalInsert(criterion1(P));
    }
    stderr.writefln("returning %.20e = %.17a", pq0.front, pq0.front);
    stderr.writefln("returning %.20e = %.17a", pq1.front, pq1.front);
    import std.typecons : tuple;
    return tuple(pq0.front, pq1.front);
}
