import std.random : uniform;
import std.algorithm : swap;

auto deepcopy(U)(in U[][] basis)
{
    auto ret = new U[][](basis.length, basis[0].length);
    foreach (i, ref row; ret)
        row[] = basis[i][];
    return ret;
}

/// Given a space by a basis, take another basis of the space so that the first vector of the output is randomly selected from the whole space.
auto randomize(U)(U[][] basis)
{
    immutable m = basis.length;
    immutable k = 1.uniform(1 << m);
    bool first = true;
    foreach (i; 0..m)
        if (k >> i & 1)
            if (first)
            {
                swap(basis[0], basis[i]);
                first = false;
            }
            else
                basis[0][] ^= basis[i][];
    return basis;
}

auto modify(U)(U[][] basis, size_t index_s, U b_idx_r)
{
    basis[0][index_s] ^= 1 << b_idx_r;
    return basis;
}

version (stand_alone):
void main()
{
    auto b = [[8, 12], [4, 6], [3, 1]];
    b.output();
    b.randomize().output();
    foreach_reverse (b_idx_r; 0..4)
        foreach (index_s; 0..2)
            modify(modify(b, index_s, b_idx_r).output(), index_s, b_idx_r);
}

auto output(U)(U[][] basis)
{
    import std.stdio;
    "%(%(%04b %)\n%)\n".writefln(basis[0..1]);
    return basis;
}
