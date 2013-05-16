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
        DigitalNet[] x, ts, ws;
        double[] wafoms;
        foreach (line; stdin.byLine())
        {
            auto ps = line.to!string().lineToBP();
            if (ps.t == 3)
                ts ~= ps;
            x ~= ps;
            wafoms ~= ps.wafom;
        }
        topN(wafoms, ts.length + 1);
        //writeln(wafoms.length);
        //writeln(ts.length);
        assert (wafoms[ts.length-1] < wafoms[ts.length]);
        foreach (ps; x)
            if (ps.wafom < wafoms[ts.length])
                ws ~= ps;
        // compare performance of ts and ws
        manytest("t-small", ts);
        manytest("wafom-small", ws);
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
}