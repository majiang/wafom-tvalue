module main;

import tvalue : tvalue;
import wafom : wafom;
import sobol : defaultSobols;
import pointset : randomPoints;

import integral : integral;
import asianoption : default_integrand;
alias integral!default_integrand tf;

import std.algorithm : min, max, reduce;


import std.stdio;

version = large_sobol;
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
        foreach (j; 25..33)
            defaultSobols(4, j, j).write_performance();
    }
}

void write_performance(R)(R P)
{
    "%d,%.15f,%.15f%s".writefln(P.save.tvalue(), P.save.wafom(), P.save.tf(), P.basis.tocsv());
}

import std.conv : to;
string tocsv(ulong[][] basis)
{
    immutable m = basis.length;
    immutable n = basis[0].length;
    string ret;
    foreach (i; 0..m)
    {
        ret ~= ",";
        foreach (j; 0..n)
        {
            ret ~= "," ~ basis[i][j].to!string();
        }
    }
    return ret;
}
