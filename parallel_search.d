import std.bigint;

BigInt countSubspaces(size_t n, size_t r)
{
    if (n < r)
        assert (false);
    auto
        ret = BigInt(1),
        u = BigInt(1) << n,
        v = BigInt(1) << r,
        w = BigInt(1);
    foreach (i; 0..r)
    {
        ret *= u - w;
        w <<= 1;
    }
    w = 1;
    foreach (i; 0..r)
    {
        ret /= v - w;
        w <<= 1;
    }
    return ret;
}

unittest
{
    import std.stdio;
    foreach (n; 0..10)
    {
        foreach (k; 0..n)
            n.countSubspaces(k).write(" ");
        writeln(1);
    }
}

import std.typecons;
alias Tuple!(size_t, size_t, size_t) Triad;

Triad[] params(immutable size_t n, immutable size_t maxgen)
{
    Triad[] ret;
    foreach (dimF; 1..n+1)
        foreach (dimR; 1..n+1)
            foreach (prec; 1..n+1)
            {
                if (dimR * prec < dimF || maxgen < (dimR * prec).countSubspaces(dimF))
                    continue; // skip
                ret ~= Triad(prec, dimR, dimF);
            }
    return ret;
}
alias uint U;
alias Tuple!(Triad, U[][], double) DN;

import std.concurrency;

auto generate(Triad param, size_t count, Tid writer)
{
    import lib.smallsearch;
    import std.container;
    import lib.wafom;
    import lib.pointset;
    import std.stdio;
    import std.format;
    import std.array;
    auto heap = (new DN[count]).heapify!((p, q) => p[2] < q[2])(0);
    foreach (basis; ReducedBasis!U(param[0], param[1], param[2]))
        heap.conditionalInsert(DN(param, basis, ShiftedBasisPoints!U(basis, param[0]).biwafom()));
    immutable file_name = text("R", param[1], "-F", param[2], "-P", param[0], ".csv");
    auto content = appender("");
    foreach (dn; heap.release())
        content ~= text(ShiftedBasisPoints!U(dn[1], param[0]).toString(), ",", dn[2], "\n");
    send(writer, Tuple!(string, string)(file_name, content.data()));
}

void main()
{
    import std.conv;
    import std.array;
    import std.string;
    import std.stdio;
    auto buf = readln().strip().split();
    immutable n = buf[0].to!size_t();
    immutable m = buf[1].to!size_t();
    immutable c = buf[2].to!size_t();
    auto threadcount = 0;
    foreach (t; n.params(m))//.parallel())
    {
        stderr.writeln("spawning " , "R", t[1], "-F", t[2], "-P", t[0]);
        spawn(&generate, t, c, thisTid);
        threadcount += 1;
    }
    while (threadcount)
    {
        auto msg = receiveOnly!(string, string);
        stderr.write("writing ", msg[0], " ...");
        File(msg[0], "w").write(msg[1]);
        threadcount -= 1;
        stderr.writeln(" done; ", threadcount, " left");
    }
}
