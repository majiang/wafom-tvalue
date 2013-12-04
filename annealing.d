import lib.wafom : bipwafom, bipmswafom;
import lib.pointset : ShiftedBasisPoints, randomVectors, randomVector, randomBits;

import std.exception : enforce;

alias uint U;
alias ShiftedBasisPoints!U PointSet;

import std.stdio : stderr;

version (stand_alone)
{
    import ui.input : getDigitalNet, getDigitalNets;
    alias getDigitalNet!U getPointSet;
    alias getDigitalNets!U getPointSets;
    import ui.output : writeWith, writeWithMany;
    import std.stdio;

    void main()
    {
        auto P = getPointSet();
        auto Q = P.yakinamashi_stay!bipmswafom(10000);
        Q.writeln();
        auto R = Q.yakinamashi_stay!bipmswafom(10000);
        R.writeln();
        auto S = R.yakinamashi_stay!bipmswafom(10000);
        S.writeln();

        P.yakinamashi_proceed!bipmswafom(10000).writeln();
        Q.yakinamashi_proceed!bipmswafom(10000).writeln();
        R.yakinamashi_proceed!bipmswafom(10000).writeln();
        S.yakinamashi_proceed!bipmswafom(10000).writeln();
    }
}

auto ddup(T)(in T[][] arr)
{
    T[][] ret;
    foreach (line; arr)
        ret ~= line.dup;
    return ret;
}

auto yakinamashi_stay(alias criterion)(PointSet P, size_t search_width)
{
    return P.yakinamashi_proceed!criterion(search_width, 0);
}

auto yakinamashi_proceed(alias criterion)(PointSet P, size_t search_width, size_t proceed_width = 1)
{
    auto best_PS = [P];
    auto best_value = criterion(P);
    stderr.writeln(best_value);
    foreach (count; 0..search_width)
    {
        auto basis = P.basis.ddup() ~ P.precision.randomVectors!U(P.dimensionR, proceed_width);
        immutable vector = P.precision.randomVector!U(P.dimensionR);
        immutable xor_at = P.dimensionF2.randomBits!uint();
        foreach (i; 0..P.dimensionF2)
            if (xor_at >> i & 1)
                foreach (j, v; vector)
                    basis[i][j] ^= v;
        auto Q = PointSet(basis, P.precision);
        immutable current_value = criterion(Q);
        if (best_value < current_value)
            continue;
        stderr.writefln("%d: %e -> %e", count, best_value, current_value);
        best_PS = [Q];
        best_value = current_value;
    }
    enforce(best_PS[0].dimensionF2 == P.dimensionF2 + proceed_width);
    return best_PS[0];
}
