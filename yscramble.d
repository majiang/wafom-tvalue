import lib.scramble : deepcopy, largeRandomIdentity;
import lib.wafom;
import lib.pointsettype;
import std.stdio;
import std.math : lg = log2;
alias U = uint;

void main(string[] args)
{
    args = args[1..$];
    if (args.length != 5)
        return stderr.writeln("yscramble dimR dimB dScramble nInitial nScramble");
    stderr.writefln("%('%s' %)", args);

    immutable size_t
        precision = 32,
        dimR = args[0].to!size_t(),
        dimB = args[1].to!size_t(),
        degree_of_scrambling = args[2].to!size_t(),
        numInitial = args[3].to!size_t(),
        numScramble = args[4].to!size_t();

    assert (numScramble <= 16384);
    foreach (i; 0..numInitial)
    {
        auto initial = randomPointSet!U(Precision(32), DimensionR(dimR), DimensionF2(dimB));
        "%s,%.15f".writef(initial, initial.bimsnrtwafom().lg());
        auto basis = initial.basis.deepcopy();
        foreach (j; 0..numScramble)
        {
            auto current = ShiftedBasisPoints!U(largeRandomIdentity!U(32, dimR, degree_of_scrambling) * basis, Precision(32)).bimsnrtwafom().lg();
            ",%.15f".writef(current);
        }
        writeln();
    }
}
