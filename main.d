module main;

import tvalue : tvalue, nu_star;
import wafom : biwafom, binrtwafom, bimsnrtwafom, bimswafom, prwafom;
import sobol : defaultSobols;
import pointset : ShiftedBasisPoints, randomVector, nonshiftedRandomBasisPoints;
import randomsearch : minimum;

import integral : bintegral;
import asianoption : default_integrand;
//alias integral!default_integrand tf;
import testfunction;


import std.algorithm : min, max, reduce, map, sort, topN;


import std.stdio;
import std.conv : to, toImpl, text;
import std.string : strip;
import std.array : split;
import walsh;

version = hamukazu;
void main()
{
    version (hamukazu)
    {
        import testfunction : Hamukazu;
        alias Hamukazu!3.f ham;
        immutable precision = 32;
        auto c = 0;
        foreach (line; stdin.byLine())
        {
            version (small) if (10 <= c++) break;
            auto buf = line.strip().split();
            uint[][] basis;
            foreach (i, x; buf)
                if (i)
                    basis ~= [x.to!uint()];
            auto P = ShiftedBasisPoints!uint(basis, precision);
            //"%s,%.15e,%.15e".writefln(buf[0], P.prwafom(), (x => x < 0 ? -x : x)(1.0 - P.bintegral!ham()));
            "%s,%s,%.15e".writefln(buf[0], P.prwafom(), P.biwafom());
        }
    }
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
        foreach (dimF2; 1..28)
        {
            stderr.writefln("dimF2 = %d", dimF2);
            auto P = ShiftedBasisPoints!uint(basis[0..dimF2], precision);
            "%d,".writef(P.dimensionF2);
            P.write_performance!(Q => bintegral!(default_integrand)(Q))();
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
        immutable u = 1024, v = 11;
        foreach (i; 0..u)
        {
            stderr.write(i, ",");
            foreach (j, l; transpose(i.walsh_function(v)))
            {
                if (i && !j) continue;
                "%15f".writef(i.walsh_coefficient_left());
                ",%15f".writef(i.walsh_coefficient_right());
                foreach (x; l)
                    ",".write(x);
                writeln();
            }
        }
        stderr.writeln();
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
        immutable size_t precision = 32, dimensionR = 4, start_dimensionF2 = 8, end_dimensionF2 = 12;
        alias ShiftedBasisPoints!uint PST;

        foreach (i; 0..100)
        {
            stderr.writefln("%d%% complete", i);
            auto P = minimum!
                (biwafom, nonshiftedRandomBasisPoints!(PST.ComponentType), PST)
                (100, precision, dimensionR, start_dimensionF2);
            void rec (PST Q)
            {
                auto R = minimum!
                    (biwafom, (PST S) => S * randomVector!(PST.ComponentType)(S.precision, S.dimensionR), PST)
                    (10, Q);
                R.write_performance!(Q => bintegral!(default_integrand)(Q))();
                if (R.dimensionF2 < end_dimensionF2)
                    rec(R);
            }
            P.write_performance!(Q => Q.bintegral!(default_integrand)())();
            rec(P);
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
        "dimF2,dimR,precision,dick,sqrt-mean-square-dick,nrt,sqrt-mean-square-nrt,t-value,basis".writeln();
        foreach (P; DN!32())
        {
            if (i % 100 == 0)
                stderr.writefln("processing %d-th point set...", i);
            "%d,%d,%d,%.15e,%.15e,%.15e,%.15e,%d".
                writef(
                       P.dimensionF2, P.dimensionR, P.precision,
                       P.biwafom(), P.bimswafom().sqrt(), P.binrtwafom(), P.bimsnrtwafom().sqrt(),
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


void write_performance(alias tf, R)(R P)
{
    "%d,%.15e,%.15f%s".writefln(P.tvalue(), P.prwafom(), tf(P), P.basis.tocsv());
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

auto tocsv(T)(const T[][] xss)
{
    string ret;
    string sep = ",";
    foreach (xs; xss)
    {
        foreach (x; xs)
        {
            ret ~= x.to!string() ~ sep;
        }
        ret ~= sep;
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

double[][m] transpose(size_t m)(double[m][] xs)
{
    double[][m] ret;
    foreach (x; xs)
        foreach (i; 0..m)
            ret[i] ~= x[i];
    return ret;
}
