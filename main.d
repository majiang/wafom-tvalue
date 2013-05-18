module main;

import tvalue : tvalue;
import wafom : wafom;
import sobol : defaultSobols;
import pointset : randomPoints;

import integral : integral;
import asianoption : default_integrand;
//alias integral!default_integrand tf;
import testfunction;

import std.algorithm : min, max, reduce, map, sort, topN;


import std.stdio;
import std.conv : to;
import std.string : strip;

version = test_funx;
void main()
{
    version (unittest_only)
    {
        "unittest passed!".writeln();
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
        DigitalNet[] x;
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
        DigitalNet[][] ts, ws;
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
            if (ts[i].length > 100)
            {
                ts[i] = ts[i][0..100];
                ws[i] = ws[i][0..100];
            }
            writeln("t-ordered");
            shifttest(ts[i]);
            writeln("wafom-ordered");
            shifttest(ws[i]);
            if (i == 2) break;
        }
    }
}
import pointset : lineToBP, DigitalNet, lesst, lessw;


void write_performance(R)(R P)
{
    "%d,%.15f,%.15f%s".writefln(P.save.tvalue(), P.save.wafom(), P.save.tf(), P.basis.tocsv());
}

string tocsv(ulong[][] basis)
{
    string ret;
    foreach (l; basis)
    {
        ret ~= ",";
        foreach (c; l)
        {
            ret ~= "," ~ c.to!string();
        }
    }
    return ret;
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
    import pointset : BasisPoints;

    void test_one(alias tf)(DigitalNet[] pss)
    {
        foreach (ps; pss)
            write(integral!(tf.f, BasisPoints)(ps.ps) - tf.I, ",");
        writeln();
    }
    void manytest(string header, DigitalNet[] pss)
    {
        header.writeln();
        pss.test_one!(hel);
        pss.test_one!(s01);
        pss.test_one!(s94);
        pss.test_one!(own);
        pss.test_one!(ra1);
        pss.test_one!(ra2);
        pss.test_one!(ra3);
        pss.test_one!(hmo);
        pss.test_one!(hme);
    }
    import pointset : shift, random_basis;
    import std.typecons : Tuple;
    void shifttest(DigitalNet[] pss)
    {
        auto shifter = random_basis(1000, 32, 4);
        stderr.writeln(shifter.length);
        assert (pss.length);
        foreach (ps; pss)
        {
            ps.shifted!hel(shifter).toString.write(",,");
            ps.shifted!s01(shifter).toString.write(",,");
            ps.shifted!s94(shifter).toString.write(",,");
            ps.shifted!own(shifter).toString.write(",,");
            ps.shifted!ra1(shifter).toString.write(",,");
            ps.shifted!ra2(shifter).toString.write(",,");
            ps.shifted!ra3(shifter).toString.write(",,");
            ps.shifted!hmo(shifter).toString.write(",,");
            ps.shifted!hme(shifter).toString.writeln();
        }
    }
    import wafom : mswafom;
    struct Stats
    {
        double wafom; double ms_wafom; double ms_error; ulong t;
        string toString()
        {
            return wafom.to!string()
                ~ "," ~ ms_wafom.to!string()
                ~ "," ~ ms_error.to!string()
                ~ "," ~ t.to!string();
        }
    }
    auto shifted(alias tf)(DigitalNet ps, ulong[][] shifter)
    {
        import std.math : sqrt;
        auto bp = ps.ps.save;
        auto mswafom = bp.save.mswafom();
        auto wafom = ps.wafom;
        auto t = ps.t;
        double sse = 0;
        foreach (s; shifter)
        {
            debug stderr.writeln("...");
            sse += (integral!(tf.f)(bp.shift(s)) - tf.I) ^^ 2;
        }
        return Stats(wafom, mswafom, sqrt(sse / shifter.length), t);
    }
}
