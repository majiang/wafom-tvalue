module lib.randomsearch;

import std.stdio;

/** Random search for an object with minimum score.

Execute a random search function repeatedly with parameters specified and yield the best result according to scoring function.

Params:
score = the scoring function, smaller score indicates better.
generate = the random search function.
count = number of search.
params = arguments passed to generate.
*/
R minimum(alias score, alias generate, R, Param...)(size_t count, Param params)
{
    R[] best = [generate(params)];
    auto best_score = score(best[$-1]);
    debug stderr.writefln("0, %.15f", best_score);
    foreach (i; 1..count)
    {
        auto current = generate(params);
        auto current_score = score(current);
        if (best_score <= current_score)
            continue;
        debug stderr.writefln("%s, %.15f", i, current_score);
        best ~= current;
        best_score = current_score;
    }
    debug stderr.writeln();
    return best[$-1];
}

debug (verbose) unittest
{
    immutable size_t precision = 32, dimensionR = 4, dimensionF2 = 10;
    import lib.wafom : biwafom;
    import lib.pointset : ShiftedBasisPoints, nonshiftedRandomBasisPoints, randomVector;
    // count, precision, dimR, dimF2
    alias ShiftedBasisPoints!uint PST;
    PST P = minimum!
        (biwafom, nonshiftedRandomBasisPoints!(PST.ComponentType), PST)
        (1000, precision, dimensionR, dimensionF2);
    PST Q = minimum!
        (biwafom,
         (PST PS) => PS * randomVector!(PST.ComponentType)(PS.precision, PS.dimensionR),
         PST)(1000, P);
}

/// Search randomly for increment and yield the one which attains minimum score.
R increment(alias score, alias incrementor, R, Param...)(R P, size_t count, Param params)
{
    R[] best = [];
    auto best_score = score(P);
    foreach (i; 0..count)
    {
        auto current = P * incrementor(params);
        auto current_score = score(current);
        if (best_score <= current_score)
            continue;
        debug stderr.writefln("%s-th try generated better score: %.15f", i, current_score);
        best ~= current;
        best_score = current_score;
    }
    return best[$-1];
}

debug (verbose) unittest
{
    immutable size_t precision = 32, dimensionR = 4, dimensionF2 = 10;
    import lib.wafom : biwafom;
    import lib.pointset : ShiftedBasisPoints, nonshiftedRandomBasisPoints, randomVector;
    auto P = nonshiftedRandomBasisPoints!uint(precision, dimensionR, dimensionF2);
    auto Q = increment!
        (biwafom, randomVector!uint, ShiftedBasisPoints!uint)
        (P, 1000, precision, dimensionR);
}

