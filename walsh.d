module walsh;
import std.stdio;

/// Compte product of f and k-th Walsh function.
Oresen product_walsh(Oresen f, size_t k)
{
    double[2][] points;
    auto i_max = f.points.length;
    auto v = k.greater_pow_of_two();
    size_t j_max = 1U << v;
    auto d = 1.0 / j_max;
    points ~= f.points[0];
    auto i = 1, j = 1;
    double y_jj = 1;
    debug "i < %d, j <= %d".writefln(i_max, j_max);
    while (i < i_max && j <= j_max)
    {
        auto x_i = f.points[i][0];
        auto y_i = f.points[i][1];
        auto x_j = j * d;
        auto y_j = ((j == j_max) ? 0 : k.wal(v, j));
        if (x_i < x_j) // (x_i, x_{j-1})
        {
            points ~= [x_i, y_i * y_jj]; debug "appended %s for x_i < x_j".writefln(points[$-1]);
            i += 1;
        }
        else if (x_i == x_j)
        {
            points ~= [x_i, y_i * y_jj]; debug "appended %s for x_i = x_j".writefln(points[$-1]);
            if (i + 1 < i_max && x_i == f.points[i + 1][0])
            {
                i += 1;
                y_i = f.points[i][1];
            }
            if (j < j_max) {points ~= [x_i, y_i * y_j]; debug "appended %s for x_i = x_j".writefln(points[$-1]);}
            i += 1;
            j += 1;
        }
        else //if (x_j < x_i)
        {
            debug (none)
            {
                "x_i = ".writeln(x_i);
                "y_i = ".writeln(y_i);
                "x_j = ".writeln(x_j);
                "y_j = ".writeln(y_j);
            }
            auto y = x_j.interpolate(f.points[i - 1], f.points[i]);
            points ~= [x_j, y * y_jj]; debug "appended %s for x_i > x_j".writefln(points[$-1]);
            if (j < j_max) {points ~= [x_j, y * y_j]; debug "appended %s for x_i > x_j".writefln(points[$-1]);}
            j += 1;
        }
        if (2 <= points.length && points[$-1] == points[$-2])
            points.length -= 1;
        debug "i = %d, j = %d.".writefln(i, j);
        y_jj = y_j;
    }
    return Oresen(points);
}

unittest
{
    import std.stdio;
    "Oresen test".writeln();
    double[2][] P = [
        [0.0, 0.0],
        [0.3, 1.0],
        [0.5, 0.0],
        [0.8, 1.0],
        [1.0, 0.0],
    ];
    foreach (p; Oresen(P).product_walsh(2).points)
    {
        p.writeln();
    }
}

auto interpolate(double x, double[2] left, double[2] right)
in
{
    assert (left[0] <= x && x <= right[0]);
}
body
{
    return (left[1] * (right[0] - x) + right[1] * (x - left[0])) / (right[0] - left[0]);
}

/// The type for piecewise affine (linear) functions.
struct Oresen
{
    double[2][] points;
    version (none)
    double integral()
    {
        double ret = 0;
        foreach (i; 0..points.length-1)
        {
            ret += (points[i+1][0]-points[i][0]) * (points[i][1] + points[i+1][1]);
        }
        return ret * 0.5;
    }
    version (none)
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

/// Hamukazu function as piecewise affine.
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

/// k-th Walsh function as oresen
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

/// k-th Walsh function value at x
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

unittest ///
{
    assert (wal(0, 0, 0) == 1); assert (wal(1, 1, 0) == 1); assert (wal(1, 1, 1) == -1);
    assert (wal(2, 2, 0) == 1); assert (wal(2, 2, 1) == -1); assert (wal(2, 2, 2) == 1); assert (wal(2, 2, 3) == -1);
    assert (wal(3, 2, 0) == 1); assert (wal(3, 2, 1) == -1); assert (wal(3, 2, 2) == -1); assert (wal(3, 2, 3) == 1);
    assert (wal(4, 3, 0) == 1);
    assert (wal(5, 3, 0) == 1);
}

/// for k, return the smallest value i which satisfy 1 << i > k
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

unittest ///
{
    assert (greater_pow_of_two(1) == 1);
    assert (greater_pow_of_two(2) == 2);
    assert (greater_pow_of_two(3) == 2);
    assert (greater_pow_of_two(4) == 3);
    assert (greater_pow_of_two(7) == 3);
    assert (greater_pow_of_two(8) == 4);
    assert (greater_pow_of_two(9) == 4);
}
