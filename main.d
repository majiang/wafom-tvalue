module main;

import tvalue : tvalue, nu_star;
import wafom : biwafom, binrtwafom, bimsnrtwafom, bimswafom, prwafom, dick_weight_enumerator_polynomial_csv;
import sobol : defaultSobols;
import pointset : ShiftedBasisPoints, randomVector, nonshiftedRandomBasisPoints, fromString;
import randomsearch : minimum;

import integral : bintegral;
import asianoption : default_integrand;
import testfunction;

import std.algorithm : min, max, reduce, map, sort, topN;

import std.stdio;
import std.string : strip;
import std.array : split, replace;
import walsh;
import std.range : isInputRange;

import std.conv : to;
template fromHex(T) if (isUnsigned!T)
{
    T fromHex(S)(S x)
    {
        return x.to!T(16);
    }
}


version = working;
version (working)
{
void write_performance(R)(R P)
{
    import std.datetime;
    "%s,%s,%.15e,%.15e".writefln(P.toString(), Clock.currTime().toSimpleString(), P.biwafom(), P.bimswafom());
}
void main(string[] args)
{
    foreach (i, w; args)
        stderr.writeln(i, ": ", w);
    if (args.length == 1)
        return args[0].split("\\")[$-1].split("/")[$-1].display_usage();
    auto command = args[1];
    args = args[2..$];
    if (command.length >= 3 && command[0..3] == "gen")
        return args.generate_point_sets();
    if (command.length >= 4 && command[0..4] == "filt")
        return args.point_sets_filters();
    if (command.length >= 3 && command[0..3] == "cbc")
        return args.component_by_component();
}
void component_by_component(string[] args)
{
    size_t start_dimensionF2, finish_dimensionF2, dimensionR, precision, initial_count, step_count, total_count;
    foreach (arg; args)
    {
        if (arg.length >= 2)
        {
            auto cmd = arg[0..2], val = arg[2..$].to!size_t();
            if (cmd == "sd") start_dimensionF2 = val;
            if (cmd == "fd") finish_dimensionF2 = val;
            if (cmd == "dr") dimensionR = val;
            if (cmd == "ic") initial_count = val;
            if (cmd == "sc") step_count = val;
            if (cmd == "tc") total_count = val;
        }
        if (arg.length >= 1 && arg[0] == 'p')
            precision = arg[1..$].to!size_t();
    }
    if (!start_dimensionF2)
    {
        stderr.writeln("Start dimension over F is by default 1.  Use sd# to change.");
        start_dimensionF2 = 1;
    }
    if (!finish_dimensionF2)
        stderr.writeln("Use fd# to specify finish dimension over F.");
    if (!dimensionR)
        stderr.writeln("Use dr# to specify dimension over R.");
    if (!step_count)
        stderr.writeln("Use sc# to specify the number of generation per step.");
    if (!total_count)
    {
        stderr.writeln("The number of iteration is by default 1.  Use tc# to change.");
        total_count = 1;
    }
    if (!initial_count)
    {
        stderr.writeln("The number of initial generation is by default the same as generation per step.  Use ic# to change.");
        initial_count = step_count;
    }
    if (!precision)
    {
        stderr.writeln("Precision is by default 32.  To change, use precision# (or prec, p in short).");
        precision = 32;
    }
    if (!(finish_dimensionF2 && dimensionR && step_count))
        return;

    void generation(T)()
    {
        alias ShiftedBasisPoints!T PST;
        auto P = minimum!
            (biwafom, nonshiftedRandomBasisPoints!T, PST)
            (initial_count, precision, dimensionR, start_dimensionF2);
        void rec(PST Q)
        {
            auto R = minimum!
                (biwafom, (PST S) => S * randomVector!T(S.precision, S.dimensionR), PST)
                (step_count, Q);
            R.write_performance();
            if (R.dimensionF2 < finish_dimensionF2)
                rec(R);
        }
        P.write_performance();
        rec(P);
    }
    foreach (i; 0..total_count)
    {
        if (precision <= 16)
            if (precision <= 8)
                generation!ubyte();
            else
                generation!ushort();
        else
            if (precision <= 32)
                generation!uint();
            else
            {
                assert (precision <= 64);
                generation!ulong();
            }
    }
}
void generate_point_sets(string[] args)
// dimR, dimF, prec=32, count=1,
{
    size_t dimensionR, dimensionF2, precision, count;
    foreach (arg; args)
    {
        if (arg.length >= 9 && arg[0..9] == "dimension")
            arg = "d" ~ arg[9..$];
        if (arg.length >= 3 && arg[0..3] == "dim")
            arg = "d" ~ arg[3..$];
        if (arg[0] == 'd')
        {
            if (arg[1] == 'R' || arg[1] == 'r')
                dimensionR = arg[2..$].to!size_t();
            if (arg[1] == 'F' || arg[1] == 'f')
                dimensionF2 = arg[2..$].to!size_t();
        }

        if (arg.length >= 9 && arg[0..9] == "precision")
            arg = "p" ~ arg[9..$];
        if (arg.length >= 4 && arg[0..4] == "prec")
            arg = "p" ~ arg[4..$];
        if (arg[0] == 'p')
            precision = arg[1..$].to!size_t();

        if (arg.length >= 5 && arg[0..5] == "count")
            arg = "c" ~ arg[5..$];
        if (arg.length >= 3 && arg[0..3] == "cnt")
            arg = "c" ~ arg[5..$];
        if (arg.length >= 1 && arg[0] == 'c')
            count = arg[1..$].to!size_t();
    }
    if (!count)
        stderr.writeln("Use count# (or cnt, c in short) to specify the number of point sets to generate.");
    if (!dimensionR)
        stderr.writeln("Use dimensionR# (or dimR, dr in short) to specify dimension of the point sets over [0..1].");
    if (!dimensionF2)
        stderr.writeln("Use dimensionF# (or dimF, df in short) to specify dimension of the point sets over Z/2Z.");
    if (!(count && dimensionR && dimensionF2))
        return;
    if (!precision)
    {
        stderr.writeln("Precision is by default 32.  To change, use precision# (or prec, p in short).");
        precision = 32;
    }

    import pointset : randomBasisPoints;
    import std.typecons : Flag;

    foreach (i; 0..count)
    {
        if (precision <= 8)
            randomBasisPoints!ubyte(precision, dimensionR, dimensionF2, Flag!"shift".no).toString().writeln();
        else if (precision <= 16)
            randomBasisPoints!ushort(precision, dimensionR, dimensionF2, Flag!"shift".no).toString().writeln();
        else if (precision <= 32)
            randomBasisPoints!uint(precision, dimensionR, dimensionF2, Flag!"shift".no).toString().writeln();
        else if (precision <= 64)
            randomBasisPoints!ulong(precision, dimensionR, dimensionF2, Flag!"shift".no).toString().writeln();
        else assert (false);
    }
}


void RapidDick(T)(string arg, T pointset)
{
    if (
        arg == "Dick" || arg == "dick" || arg == "D" || arg == "d" ||
        arg == "RapidDick" || arg == "rDick" || arg == "rd" || arg == "Rd" || arg == "rD" || arg == "RD" ||
        arg == "WAFOM" || arg == "wafom" || arg == "W" || arg == "w")
        "%.15e".writef(pointset.biwafom());
}
void PrecDick(T)(string arg, T pointset)
{
    if (
        arg == "PrecDick" || arg == "pDick" || arg == "pd" || arg == "Pd" || arg == "pD" || arg == "PD" ||
        arg == "PrecWAFOM" || arg == "pWAFOM" || arg == "pW" || "Pw")
        "%.15e".writef(pointset.prwafom());
}
void RMSDick(T)(string arg, T pointset)
{
    if (
        arg == "RMSDick" || arg == "MSDick" || arg == "rmsd" || arg == "msd" || arg == "msw")
        "%.15e".writef(pointset.bimswafom());
}
void DickWEP(T)(string arg, T pointset)
{
    import wafom : dick_weight_enumerator_polynomial;
    if (
        arg == "DickWEP" || arg == "dwep")
        "%s".writef(pointset.dick_weight_enumerator_polynomial());
}
void NRT(T)(string arg, T pointset)
{
    if (
        arg == "NRT" || arg == "nrt")
        "%.15e".writef(pointset.binrtwafom());
}
void RMSNRT(T)(string arg, T pointset)
{
    if (
        arg == "RMSNRT" || arg == "MSNRT" || arg == "rmsnrt" || arg == "msnrt")
        "%.15e".writef(pointset.bimsnrtwafom());
}
void tValue(T)(string arg, T pointset)
{
    if (
        arg == "t" || arg == "tvalue" || arg == "tv")
        "%d".writef(pointset.tvalue());
}
void point_set_filters(T)(T pointset, string[] args)
{
    pointset.toString().write();
    auto sep = ",";
    foreach (arg; args)
    {
        sep.write();
        arg.RapidDick(pointset);
        arg.PrecDick(pointset);
        arg.RMSDick(pointset);
        arg.DickWEP(pointset);
        arg.NRT(pointset);
        arg.RMSNRT(pointset);
        arg.tValue(pointset);
    }
    writeln();
}

void point_sets_filters(string[] args)
{
    foreach (line; stdin.byLine)
    {
        auto precision = line.strip().split()[0].to!size_t();
        if (precision <= 16)
        {
            if (precision <= 8)
                line.fromString!(ubyte).point_set_filters(args);
            else
                line.fromString!(ushort).point_set_filters(args);
        }
        else
        {
            if (precision <= 32)
                line.fromString!(uint).point_set_filters(args);
            else
            {
                assert (precision <= 64);
                line.fromString!(ulong).point_set_filters(args);
            }
        }
    }
}

void display_usage(string arg)
{
    stderr.writeln("usage:");
    stderr.writeln(arg ~ " generate (gen): generate point sets");
    stderr.writeln(arg ~ " filter (filt): compute some function for each point set");
    stderr.writeln(arg ~ " cbc: perform component by component construction of digital net");
}
}
else void main()
{
    version (wep)
    {
        foreach (line; stdin.byLine())
        {
            auto x = line ~ (line.fromString!uint().dick_weight_enumerator_polynomial_csv());
            x.writeln();
            stdout.flush();
            //stderr.writeln(x);
        }
    }
    version (random_search_distribution)
    {
        foreach (line; stdin.byLine())
            line.fromString!uint().write_wafoms(16384);
    }
    version (hamukazu)
    {
        version (small) int c;
        foreach (line; stdin.byLine())
        {
            version (small) if (10 <= c++) break;
            line.fromString!uint().write_performance();
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
                basis[i][j] = buf[i + dimensionF2 * j].fromHex!uint() >> 2;
                assert (basis[i][j] < 1U << 30);
            }
        }
        foreach (dimF2; 1..28)
        {
            stderr.writefln("dimF2 = %d", dimF2);
            auto P = ShiftedBasisPoints!uint(basis[0..dimF2], precision);
            P.write_performance!(Q => bintegral!(default_integrand)(Q))();
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
        immutable size_t precision = 32, dimensionR = 4, start_dimensionF2 = 8, end_dimensionF2 = 16;
        alias ShiftedBasisPoints!uint PST;

        foreach (i; 0..10)
        {
            stderr.writefln("%d%% complete", i * 10);
            auto P = minimum!
                (biwafom, nonshiftedRandomBasisPoints!(PST.ComponentType), PST)
                (10000, precision, dimensionR, start_dimensionF2);
            void rec (PST Q)
            {
                auto R = minimum!
                    (biwafom, (PST S) => S * randomVector!(PST.ComponentType)(S.precision, S.dimensionR), PST)
                    (10000, Q);
                R.write_performance();
                if (R.dimensionF2 < end_dimensionF2)
                    rec(R);
            }
            P.write_performance();
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
                return;
        }
    }
    version (small_wafom)
    {
        foreach (i; 0..100)
        {
            auto P = nonshiftedRandomBasisPoints!ubyte(6, 2, 6);
            auto w = P.biwafom();
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
        DigitalNet!uint[] x;
        double[] wafoms;
        ulong[] tvalues;
        foreach (line; stdin.byLine())
        {
            auto ps = line.to!string().lineToDN!uint();
            x ~= ps;
            wafoms ~= ps.wafom;
            tvalues ~= ps.tvalue;
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
        DigitalNet!uint[][] ts, ws;
        ts.length = thresholdt.length;
        ws.length = thresholdw.length;
        foreach (ps; x)
        {
            foreach (i, t; thresholdt)
            {
                if (ps.tvalue < t)
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
            manytest!uint("t-ordered", ts[i]);
            manytest!uint("wafom-ordered", ws[i]);
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
version(none)
{
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
}
void write_wafoms(R)(R P, size_t count)
{
    auto w = P.biwafom();
    foreach (i; 0..count)
    {
        auto v = randomVector!(R.ComponentType)(P.precision, P.dimensionR);
        "%.15e,%s".writefln((P + v).biwafom() + w, v.toSSV());
    }
    writeln();
}

version (sharase)
void write_performance(alias f, R)(R P)
{
    "%s,%.15e,%.15e,%d".writefln(P.toString(), P.biwafom(), P.bintegral!default_integrand(), P.tvalue());
}
else version (random_search)
void write_performance(R)(R P)
{
    import std.datetime;
    "%s,%s,%.15e,%.15e".writefln(P.toString(), Clock.currTime().toSimpleString(), P.biwafom(), P.bimswafom());
}
else version (hamukazu)
void write_performance(R)(R P)
{
    import std.datetime;
    import testfunction : Hamukazu;
    alias Hamukazu!3.f ham;
    "%s,%s,%.15e,%.15e".writefln(P.toString(), Clock.currTime().toSimpleString(), P.prwafom(), (x => x < 0 ? -x : x)(1.0 - P.bintegral!ham()));
}
else version (large_sobol)
{
void write_performance(R)(R P)
{
    import std.typecons : Tuple;
    import pointset : InfoPointSet;
    alias Tuple!(double, "wafom", ulong, "tvalue", double, "aoPrice") InfoType;
    alias InfoPointSet!(R, InfoType) IPS;
    IPS(P, InfoType(P.biwafom(), P.tvalue(), P.bintegral!default_integrand())).toString().writeln();
}
}
else
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

auto toSSV(T)(const T[] xs)
{
    string ret;
    string sep = "";
    foreach (x; xs)
    {
        ret ~= sep ~ x.to!string();
        sep = " ";
    }
    return ret;
}

auto tocsv(T)(T xs) if (isInputRange!T && !(isInputRangeOfInputRange!T))
{
    string ret;
    string sep = ",";
    foreach (x; xs)
    {
        ret ~= x.to!string() ~ sep;
    }
    return ret;
}

auto tocsv(T)(T xss) if (isInputRangeOfInputRange!T)
{
    string ret;
    string sep = ",";
    foreach (xs; xss)
    {
        ret ~= xs.tocsv() ~ sep;
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
    return -tf.I + P.bintegral!(tf.f, PointSetType)();
}
/** ditto
Params:
Ps = point sets.
*/
auto integrationErrors(alias tf, PointSetTypeRange)(PointSetTypeRange Ps)
{
    return Ps.map!(integrationError!(tf, ElementType!PointSetTypeRange))();
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
    import std.typecons : Tuple;
    alias Tuple!(double, "wafom", ulong, "tvalue") InfoType;

    import std.traits : isUnsigned;
    template DigitalNet(T) if (isUnsigned!T)
    {
        import pointset : InfoPointSet;
        alias InfoPointSet!(ShiftedBasisPoints!T, InfoType) DigitalNet;
    }
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
    void manytest(T)(string header, DigitalNet!T[] dns)
    {
        alias tocsv trf;
        //alias squareRootMeanSquare trf;
        auto pss = dns.map!(x => x.pointSet);
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
    version (none)
    {
    void shifttest(DigitalNet!ulong[] pss)
    {
        import pointset : randomVectors;
        auto shift = randomVectors!ulong(32, 4, 10);
        stderr.writeln(shift.length);
        assert (pss.length);
        foreach (ps; pss)
        {
            auto P = ps.pointSet;
            // old: [wafom, squarerootmeansquarewafom, squarerootmeansquareerror, t]
            // new: [wafom, srmsw, t] [srmse]
            // rationale: wafom, srmsw, t is independent of test funx.
            ps.wafom.write(",");
            P.bimswafom().write(",");
            ps.tvalue.write(",,");
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
    auto srmse(alias tf, PointSetType)(PointSetType P, ulong[][] shifts)
    {
        return P.shifteds(shifts).integrationErrors!tf().squareRootMeanSquare();
    }
    /** Shift a point set by each element of shifts.

    Params:
    P = point set.
    shifts = shifts
    */
    auto shifteds(PointSetType, ShifterRange)(PointSetType P, ShifterRange shifts)
    {
        import std.array : empty, front, popFront;
        struct R
        {
            @property bool empty()
            {
                return shifts.empty;
            }
            @property PointSetType front()
            {
                return P.shifted(shifts.front);
            }
            void popFront()
            {
                shifts.popFront();
            }
        }
        return R();
        //return shifts.map!(x => P.shifted(x));
    }
    }
    auto lineToDN(T)(string line, size_t precision = size_t.max) if (isUnsigned!T)
    {
        InfoType info;
        T[][] basis;
        foreach (i, bufs; line.strip().split(",,"))
        {
            auto buf = bufs.split(",");
            if (!i)
            {
                info.tvalue = buf[0].to!ulong();
                info.wafom = buf[1].to!double();
                continue;
            }
            basis.length += 1;
            foreach (s; bufs.split(","))
                basis[$-1] ~= s.to!T();
        }
        import pointset : guess_precision;
        return DigitalNet!T(ShiftedBasisPoints!T(basis.transpose(), precision = size_t.max ? basis.guess_precision() : precision), info);
    }
}

T[][] transpose(T)(T[][] xs)
{
    import std.exception : enforce;
    if (xs.length == 0)
        return [];
    immutable n = xs[0].length;
    auto ret = new T[][n];
    foreach (x; xs)
    {
        enforce(x.length == n);
        foreach (i; 0..n)
            ret[i] ~= x[i];
    }
    return ret;
}

double[][m] transpose(size_t m)(double[m][] xs)
{
    double[][m] ret;
    foreach (x; xs)
        foreach (i; 0..m)
            ret[i] ~= x[i];
    return ret;
}

template isInputRangeOfInputRange(R)
{
    enum bool isInputRangeOfInputRange = is(typeof((void)
    {
        R r = void;
        if (r.empty){}
        r.popFront();
        auto e = r.front;
        if (e.empty){}
        e.popFront();
        auto ee = e.front;
    }));
}
