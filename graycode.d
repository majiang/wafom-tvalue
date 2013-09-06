module lib.graycode;
debug import std.stdio;

T graycode(T)(T x)
{
    return x ^ (x >> 1);
}

T edocyarg(T)(T x)
{
    size_t i = 1;
    while (x >> i)
    {
        x = x ^ (x >> i);
        i += 1;
    }
    return x;
}

unittest
{
    foreach (x, y; [0, 1, 3, 2, 6, 7, 5, 4])
    {
        assert (x.graycode == y);
        assert (y.edocyarg == x);
    }
    debug (verbose) "graycode & edocyarg: unittest passed!".writeln();
}

size_t bottom_zeros(ulong x)
{
    assert (x);
    size_t ret;
    while ((x & 1) == 0)
    {
        x >>= 1;
        ret += 1;
    }
    return ret;
}

unittest
{
    auto p = [1, 2, 3, 4, 5, 6, 7, 8, 16, 32];
    auto q = [0, 1, 0, 2, 0, 1, 0, 3, 4, 5];
    foreach (i, x; p)
    {
        assert (x.bottom_zeros() == q[i]);
    }
    debug (verbose) "bottom_zeros: unittest passed!".writeln();
}
