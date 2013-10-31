module main;

void main()
{
    import std.stdio;
    import std.string : strip;
    import std.array : split;
    import std.conv : to;
    import std.math : isNaN;

    import lib.wafom : criterion = bipwafom;
    alias uint U;

    stderr.writeln("prec dimR dimF prep rate count [abs]");
    auto buf = readln().strip().split();
    immutable
        prec = buf[0].to!size_t(),
        dimR = buf[1].to!size_t(),
        dimF = buf[2].to!size_t(),
        prep = buf[3].to!size_t(),
        rate = buf[4].to!size_t(),
        count = buf[5].to!size_t();
    double abs;
    if (6 < buf.length)
        abs = buf[6].to!double();

    foreach (P; generation!(U, criterion)(prec, dimR, dimF, isNaN(abs) ? 
    preparation!(U, criterion)(prec, dimR, dimF, prep, rate) : abs
    , count))
        P.toString().writeln();
}

import lib.pointset : nonshiftedRandomBasisPoints;

version = verbose;
version (verbose)
    auto milestone = [0.00, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.98, 0.99, 1.00];
version (medium)
    auto milestone = [0.00, 0.03, 0.09, 0.27, 0.81];

auto generation(U, alias criterion)
(
    in size_t prec, in size_t dimR, in size_t dimF,
    in double threshold, in size_t count
)
{
    static struct R
    {
        size_t length;
        immutable size_t prec;
        immutable size_t dimR;
        immutable size_t dimF;
        immutable double threshold;
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
            if (threshold < criterion(P))
                goto retry;
            return P;
        }
    }
    return R(count, prec, dimR, dimF, threshold, count);
}

auto preparation(U, alias criterion)
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
    auto pq = (new double[rate]).heapify(0);
    foreach (i; 0..prep)
        pq.conditionalInsert(criterion(nonshiftedRandomBasisPoints!U(prec, dimR, dimF)));
    stderr.writefln("returning %.20e = %.17a", pq.front, pq.front);
    return pq.front;
}
