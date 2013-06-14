module walsh;
import std.stdio;

version (developing){
Oresen walsh(Oresen f, size_t k)
{
    immutable n = f.points.length;
    immutable v = k.greater_pow_of_two();
    immutable V = 1 << v;
    auto d = 1.0 / V;
    double[2][] ret;
    int i, j;
    while (i < V && j < f.points.length)
    {
        auto x = i * d;
        if (x == f.points[j][0])
        {
            
        }
        if (x < f.points[j][0])
        {

            i += 1; continue;
        }
        if (f.points[j][0] < x)
        {

            j += 1; continue;
        }
    }
    return Oresen(ret);
}

struct Oresen
{
    double[2][] points;
    double integral()
    {
        double ret = 0;
        foreach (i; 0..points.length-1)
        {
            ret += (points[i+1][0]-points[i][0]) * (points[i][1] + points[i+1][1]);
        }
        return ret * 0.5;
    }
    double opCall(double x)
    {
        foreach (i; 0..points.length)
        {
            if (points[i][0] < x) continue;
            if (points[i][0] == x)
            {
                if (i + 1 < points.length && points[i+1][0] == x)
                    return (points[i][1] + points[i+1][1]) * 0.5;
                return points[i][1];
            }
            assert (0 < i);
            return (points[i][1] * (x - points[i-1][0]) + points[i-1][1] * (points[i][0] - x)) / (points[i][0] - points[i-1][0]);
        }
        assert (false);
    }
}

auto Hamukazu(size_t n)
{
    double[2][] points;
    foreach (i; 0..(n+1))
    {
        double x = i; x /= n;
        if (0 < i) points ~= [x, 2];
        if (i < n) points ~= [x, 0];
    }
    return Oresen(points);
}
}

/** walsh coefficient of [0..1/3].
*/
double walsh_coefficient_left(size_t k)
{
    double t = 1.0 / 3.0;
    double ret = 0;
    if (k == 0) return t;
    immutable v = k.greater_pow_of_two();
    immutable V = 1 << v;
    immutable d = 1.0 / V;
    foreach (i; 0..V)
    {
        double cur;
        if (3 * i < V && V < 3 * (i + 1))
        {
            version (none) stderr.writefln("%d %s - break", i, k.wal(v, i) > 0 ? "positive" : "negative");
            return ret + (t - i * d) * k.wal(v, i);
        }
        version (none) stderr.writefln("%d %s", i, k.wal(v, i) > 0 ? "positive" : "negative");
        ret += d * k.wal(v, i);
    }
    assert (false);
}

/** walsh coefficient of [2/3..1].
*/
double walsh_coefficient_right(size_t k)
{
    double t = 2.0 / 3.0;
    double ret;
    if (k == 0) return 1.0 / 3.0;
    immutable v = k.greater_pow_of_two();
    immutable V = 1 << v;
    immutable d = 1.0 / V;
    foreach (i; 0..V)
    {
        double cur;
        if (3 * i < 2 * V && 2 * V < 3 * (i + 1))
            ret = ((i+1) * d - t) * k.wal(v, i);
        else if (2 * V < 3 * i)
            ret += d * k.wal(v, i);
    }
    return ret;
}

/** walsh coefficients of hamukazu(3) : x => 2 * (3x - floor 3x)

Bugs:
wrong return value
*/
deprecated double walsh_coefficient_3(size_t k)
{
    double ret = 0;
    if (k == 0) return 1;
    immutable v = k.greater_pow_of_two();
    immutable V = 1 << v;
    foreach (i; 0..V)
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

double[2][] walsh_function(size_t k, size_t v)
{
    double[2][] points;
    foreach (i; 0..((1<<v)+1))
    {
        double x = i * (1.0 / (1 << v));
        if (0 < i) points ~= [x, k.wal(v, i - 1)];
        if (i < (1 << v)) points ~= [x, k.wal(v, i)];
    }
    return points;
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
    size_t ret;
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