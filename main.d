module main;

import tvalue : tvalue, nu_star;
import wafom : Bisect, wafom, biwafom, nrtwafom, msnrtwafom, mswafom;
import sobol : defaultSobols;
import pointset : randomPoints, ShiftedBasisPoints, randomVectors;

import integral : bintegral;
import asianoption : default_integrand;
//alias integral!default_integrand tf;
import testfunction;


import std.algorithm : min, max, reduce, map, sort, topN;


import std.stdio;
import std.conv : to, toImpl;
import std.string : strip;
import std.array : split;
import walsh;

//version = nu;
version = calc_all;
void main()
{
    version (sharase)
    {
        immutable precision = 30;
        immutable dimensionF2 = 28;
        immutable dimensionR = 4;
        uint[][] basis;
        basis.length = dimensionF2;
        foreach (ref l; basis)
        {
            l.length = dimensionR;
        }
        auto buf = readln().strip().split();//.map!(s => s.toImpl!(uint, string)(16));
        foreach (i; 0..dimensionF2)
        {
            foreach (j; 0..dimensionR)
            {
                basis[i][j] = buf[i + dimensionF2 * j].toImpl!(uint, string)(16) >> 2;
                assert (basis[i][j] < 1U << 30);
            }
        }
        foreach (dimF2; 1..17)
        {
            stderr.writefln("dimF2 = %d", dimF2);
            auto P = ShiftedBasisPoints!uint(basis[0..dimF2], precision);
            "%d %.15f".writefln(P.dimensionF2, P.biwafom());
        }
    }
    version (nu)
    {
        foreach (ushort i; 0..1024)
        {
            "%d:%d:%d".writefln(i, mu_star!ushort(i, 10), nu_star(i, 10));
        }
    }
    version (walsh)
    {
        stderr.writeln("walsh start.");
        foreach (i; 0..16384)
        {
            "%.15f".writefln(walsh_coefficient_3(i));
        }
        stderr.writeln("walsh end.");
        return;
    }
    version (unittest_only)
    {
        stderr.writeln("unittest passed!");
        return;
    }
    version (random_search)
    {
        foreach (i; 0..100)
        {
            stderr.writefln("%d%% complete", i);
            foreach (j; 0..1000)
                randomPoints(4, 32, 12).write_performance();
        }
    }
    version (large_sobol)
    {
        auto stoptime = readln().strip().to!int();
        foreach (j; 2..33)
        {
            defaultSobols(4, j, j).write_performance();
            if (j == stoptime) 
            {
                readln();
                return;
            }
        }
    }
    version (small_wafom)
    {
        foreach (i; 0..100)
        {
            auto P = randomPoints(2, 6, 6);
            auto w = P.save.wafom();
            writeln(i, ",", w);
            auto f = File("small-" ~ i.to!string() ~ "-" ~ (w * 10000).to!int().to!string() ~ ".csv", "w");
            f.writeln(w);
            foreach (x; P)
            {
                f.writefln("%d,%d", x[0], x[1]);
            }
        }
    }
    version (test_funx)
    {
        DigitalNet!ulong[] x;
        double[] wafoms;
        ulong[] tvalues;
        foreach (line; stdin.byLine())
        {
            auto ps = line.to!string().lineToBP();
            x ~= ps;
            wafoms ~= ps.wafom;
            tvalues ~= ps.t;
        }
        sort(tvalues);
        sort(wafoms);
        ulong current_t = 0;
        ulong[] thresholdt;
        double[] thresholdw;
        foreach (i; 0..x.length)
        {
            if (current_t < tvalues[i])
            {
                current_t = tvalues[i];
                thresholdt ~= tvalues[i];
                thresholdw ~= wafoms[i];
            }
        }
        thresholdt ~= tvalues[$ - 1] + 1;
        thresholdw ~= wafoms[$ - 1] * 2;
        DigitalNet!ulong[][] ts, ws;
        ts.length = thresholdt.length;
        ws.length = thresholdw.length;
        foreach (ps; x)
        {
            foreach (i, t; thresholdt)
            {
                if (ps.t < t)
                {
                    ts[i] ~= ps;
                    break;
                }
            }
            foreach (i, w; thresholdw)
            {
                if (ps.wafom < w)
                {
                    ws[i] ~= ps;
                    break;
                }
            }
        }
        foreach (l; ts)
        {
            stderr.write(l.length, ", ");
        }
        stderr.writeln();
        foreach (l; ws)
        {
            stderr.write(l.length, ", ");
        }
        stderr.writeln();
        foreach (i; 1..thresholdw.length)
        {
            //if (ts[i].length > 100)
            //{
            //    ts[i] = ts[i][0..100];
            //    ws[i] = ws[i][0..100];
            //}
            manytest("t-ordered", ts[i]);
            manytest("wafom-ordered", ws[i]);
            //if (i == 2) break;
        }
    }
    version (calc_all)
    {
        auto i = 0;
        "dimF2,dimR,precision,normal-dick,bisect-dick,Bisect-dick,sqrt-mean-square-dick,nrt,sqrt-mean-square-nrt,t-value,basis".writeln();
        foreach (P; DN!32())
        {
            //if (i % 100 == 0)
                stderr.writefln("processing %d-th point set...", i);
            "%d,%d,%d,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%d".
                writef(
                       P.dimensionF2, P.dimensionR, P.precision,
                       P.wafom(), P.biwafom(), P.Bisect!wafom(),
                       P.mswafom().sqrt(), P.Bisect!nrtwafom(), P.Bisect!msnrtwafom().sqrt(),
                       P.tvalue());
            foreach (l; P.basis) foreach (x; l) ",".
                write(x);
            writeln();
            i += 1;
        }
    }
}
import pointset : lineToBP, DigitalNet;

auto DN(size_t precision)()
{
    ShiftedBasisPoints!uint[] ret;
    foreach (line; stdin.byLine())
    {
        auto ps = line.to!string().lineToBP!uint(precision).ps;
        ret ~= ps;
    }
    return ret;
}


void write_performance(R)(R P)
{
    "%d,%.15f,%.15f%s".writefln(P.tvalue(), P.wafom(), P.tf(), P.basis.tocsv());
}

version (none) auto tocsv(T)(T xss) if (isInputRange!T && isInputRange!(ElementType!T) && !isInputRange!(ElementType!(ElementType!T)))
{
    string ret;
    foreach (xs; xss)
    {
        ret ~= ",";
        foreach (x; xs)
        {
            ret ~= "," ~ x.to!string();
        }
    }
    return ret;
}

import std.range : isInputRange;
auto tocsv(T)(T xs) if (isInputRange!T && !isInputRange!(ElementType!T))
{
    string ret;
    string sep = "";
    foreach (x; xs)
    {
        ret ~= sep ~ x.to!string();
        sep = ",";
    }
    return ret;
}

import std.range : ElementType;

/** Calculate integration error of a function by a point set.

Params:
tf = an alias. the integrand is given as tf.f, with its true integrated value tf.I.
P = a point set.
*/
double integrationError(alias tf, PointSetType)(PointSetType P)
{
    return -tf.I + P.bintegral!(tf.f, ElementType!(ElementType!PointSetType), PointSetType)();
}
/** ditto
Params:
Ps = point sets.
*/
auto integrationErrors(alias tf, PointSetTypeRange)(PointSetTypeRange Ps)
{
    return Ps.map!(integrationError!(tf, ElementType!PointSetTypeRange))();
}
/** Shift a point set by each element of shifts.

Params:
P = point set.
shifts = shifts
*/
auto shifteds(PointSetType, ShifterRange)(PointSetType P, ShifterRange shifts)
{
    return shifts.map!(x => P.shifted(x));
}

import std.math : sqrt;
auto squareRootMeanSquare(NumericRange)(NumericRange r)
{
    ElementType!NumericRange sum = 0;
    ulong count;
    foreach (e; r)
    {
        sum += e * e;
        count += 1;
    }
    return (sum / count).sqrt();
}


version (test_funx) 
{
    import testfunction;

    alias Hellekalek!(1.1, 1.7, 2.3, 2.9) hel;
    alias Sobol2001!(0.0, 0.0, 3.0, 3.0) s01;
    alias Sobol1994!4 s94;
    alias Owen!4 own;
    alias RooArnold1!4 ra1;
    alias RooArnold2!4 ra2;
    alias RooArnold3!4 ra3;
    alias Hamukazu!(5, 7, 11, 13) hmo;
    alias Hamukazu!(8, 8, 8, 8) hme;

    auto test_one(alias tf, alias transformer, PointSetTypeRange)(PointSetTypeRange pss)
    {
        return transformer(pss.integrationErrors!tf());
    }
    void manytest(string header, DigitalNet!ulong[] dns)
    {
        alias tocsv trf;
        //alias squareRootMeanSquare trf;
        auto pss = dns.map!(x => x.ps);
        header.writeln();
        pss.test_one!(hel, trf).writeln();
        pss.test_one!(s01, trf).writeln();
        pss.test_one!(s94, trf).writeln();
        pss.test_one!(own, trf).writeln();
        pss.test_one!(ra1, trf).writeln();
        pss.test_one!(ra2, trf).writeln();
        pss.test_one!(ra3, trf).writeln();
        pss.test_one!(hmo, trf).writeln();
        pss.test_one!(hme, trf).writeln();
    }
    import std.typecons : Tuple;
    void shifttest(DigitalNet!ulong[] pss)
    {
        auto shift = randomVectors!ulong(32, 4, 10);
        stderr.writeln(shift.length);
        assert (pss.length);
        foreach (ps; pss)
        {
            auto P = ps.ps;
            // old: [wafom, squarerootmeansquarewafom, squarerootmeansquareerror, t]
            // new: [wafom, srmsw, t] [srmse]
            // rationale: wafom, srmsw, t is independent of test funx.
            ps.wafom.write(",");
            P.mswafom().write(",");
            ps.t.write(",,");
            P.srmse!hel(shift).write(",");
            P.srmse!s01(shift).write(",");
            P.srmse!s94(shift).write(",");
            P.srmse!own(shift).write(",");
            P.srmse!ra1(shift).write(",");
            P.srmse!ra2(shift).write(",");
            P.srmse!ra3(shift).write(",");
            P.srmse!hmo(shift).write(",");
            P.srmse!hme(shift).writeln();
        }
    }
    import wafom : mswafom;
    import std.conv : text;
    auto srmse(alias tf, PointSetType)(PointSetType P, ulong[][] shifts)
    {
        return P.shifteds(shifts).integrationErrors!tf().squareRootMeanSquare();
    }
}
