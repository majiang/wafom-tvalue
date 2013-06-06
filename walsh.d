module walsh;

/// walsh coefficients of hamukazu(3) : x => 2 * (3x - floor 3x)
double walsh_coefficient_3(size_t k)
{
    double ret = 0;
    if (k == 0) return 1;
    auto v = k.greater_pow_of_two();
    auto V = 1 << v;
    foreach (i; 0..(1 << v))
    {
        double cur = 0;
        if (3 * i < V && V < 3 * (i + 1))
        {
            if (V % 3 == 1)
                cur = 2 * V + 3; // (V - 1) + V + 0 + 1 + 1 + 2
            else if (V % 3 == 2)
                cur = 4 * V - 3; // (V - 2) + (V - 1) + (V - 1) + V + 0 + 1
            cur /= 3;
        }
        else if (3 * i < 2 * V && 2 * V < 3 * (i + 1))
        {
            if (2 * V % 3 == 1)
                cur = 2 * V + 3; // (V - 1) + V + 0 + 1 + 1 + 2
            else if (2 * V % 3 == 2)
                cur = 4 * V - 3; // (V - 2) + (V - 1) + (V - 1) + V + 0 + 1
            cur /= 3;
        }
        else
        {
            cur = 3 * i % V + 3 * (i + 1) % V;
        }
        ret += cur * k.wal(v, i);
    }
    return ret / (V * V);
}

int wal(size_t k, size_t v, size_t x)
in
{
    assert (k < (1 << v));
    assert (x < (1 << v));
}
body
{
    int ret = 1;
    foreach (i; 0..v)
    {
        if (k >> i & x >> (v-i-1) & 1)
            ret *= -1;
    }
    return ret;
}

unittest
{
    assert (wal(0, 0, 0) == 1); assert (wal(1, 1, 0) == 1); assert (wal(1, 1, 1) == -1);
    assert (wal(2, 2, 0) == 1); assert (wal(2, 2, 1) == -1); assert (wal(2, 2, 2) == 1); assert (wal(2, 2, 3) == -1);
    assert (wal(3, 2, 0) == 1); assert (wal(3, 2, 1) == -1); assert (wal(3, 2, 2) == -1); assert (wal(3, 2, 3) == 1);
    assert (wal(4, 3, 0) == 1);
    assert (wal(5, 3, 0) == 1);
}

size_t greater_pow_of_two(T)(T k)
{
    size_t ret = 0;
    T tmp = 1;
    while (tmp <= k)
    {
        tmp <<= 1;
        ret += 1;
    }
    return ret;
}

unittest
{
    assert (greater_pow_of_two(1) == 1);
    assert (greater_pow_of_two(2) == 2);
    assert (greater_pow_of_two(3) == 2);
    assert (greater_pow_of_two(4) == 3);
    assert (greater_pow_of_two(7) == 3);
    assert (greater_pow_of_two(8) == 4);
    assert (greater_pow_of_two(9) == 4);
}