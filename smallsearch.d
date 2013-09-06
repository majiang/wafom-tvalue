module smallsearch;

import std.traits : isUnsigned;
import std.exception : enforce;

auto ReducedBasis(U)(in size_t precision, in size_t dimensionR, in size_t dimensionF2) if (isUnsigned!U)
{
    struct Result
    {
        int opApply(int delegate(ref U[][]) dg)
        {
            foreach (m; (precision * dimensionR).ReducedMatrices(dimensionF2))
            {
                U[][] basis;
                basis.length = m.length;
                foreach (i, l; m)
                    foreach (j; 0..dimensionR)
                        basis[i] ~= cast(U)((l >> (j * precision)) & ((1UL << precision) - 1));
                if (auto result = dg(basis))
                    return result;
            }
            return 0;
        }
    }
    return Result();
}

auto ReducedMatrices(in size_t numColumns, in size_t numRows)
{
    import std.algorithm : filter;
    struct Result
    {
        int opApply(int delegate(ref ulong[]) dg)
        {
            foreach (smallest; numColumns.DimensionalVectorsWithZerosAt([]).filter!(v => 0 < v))
                if (auto result = opPartialApply([smallest], dg))
                    return result;
            return 0;
        }
        int opPartialApply(ulong[] partial, int delegate(ref ulong[]) dg)
        {
            if (partial.length == numRows)
                return dg(partial);
            size_t[] zeroPositions;
            foreach (j; partial)
                zeroPositions ~= msb(j);
            immutable m = 1UL << zeroPositions[$-1];
                foreach (next; numColumns.DimensionalVectorsWithZerosAt(zeroPositions).filter!(v => m < v))
                    if (auto result = opPartialApply(partial ~ next, dg))
                        return result;
            return 0;
        }
        size_t msb(ulong j)
        {
            if (j == 1)
                return 0;
            return msb(j >> 1) + 1;
        }
    }
    return Result();
}

auto DimensionalVectorsWithZerosAt(in size_t dimension, in size_t[] zeroPositions)
{
    struct Result
    {
        ulong[2][] zero;
        ulong i;
        ulong max;
        @property bool empty(){return max <= i;}
        @property auto front()
        {
            auto current = i;
            foreach (u; zero)
            {
                if (current & u[0])
                {
                    current &= ~u[0];
                    current |= u[1];
                }
            }
            return current;
        }
        void popFront(){i += 1;}
    }
    enforce(dimension > zeroPositions.length);
    enforce(dimension - zeroPositions.length < 64);
    return Result(zeroPositions.toBits(dimension), 0, 1UL << (dimension - zeroPositions.length));
}

ulong[2][] toBits(in size_t[] zeroPositions, in size_t dimension)
{
    ulong[2][] ret;
    immutable n = zeroPositions.length;
    foreach (i, p; zeroPositions)
        ret ~= [1UL << p, 1UL << (dimension - n + i)];
    return ret;
}

void not_main()
{
    import std.stdio;
    import std.algorithm : filter;
    foreach (basis; 4.ReducedBasis!ubyte(2, 2))
    {
        foreach (v; basis)
        {
            foreach_reverse (l; v)
                "%04b".writef(l);
            writeln();
        }
        writeln();
    }
    version (none) foreach (m; ReducedMatrices(4, 2))
    {
        foreach (j; m)
            "%04b".writefln(j);
        writeln();
    }
    version (none) foreach (j; 6.DimensionalVectorsWithZerosAt([2,3]).filter!(j => j > (1UL << 3)))
        "%06b".writefln(j);
}
