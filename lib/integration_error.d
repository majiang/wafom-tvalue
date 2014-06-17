module lib.integration_error;

import lib.integral;

import std.range : ElementType, hasLength, array;
import std.traits : ReturnType;
import std.algorithm : map, reduce;
import std.math : sqrt;

/** Calculate integration error of a function by a point set.

Params:
tf = an alias. the integrand is given as tf.f, with its true integrated value tf.I.
P = a point set.
*/
real integrationError(alias tf, PointSetType)(PointSetType P)
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

/// sums a NumericRange carefully.
auto sumFromSmallerPositive(NumericRange)(NumericRange r)
    if (is (ElementType!NumericRange : real))
{
    ElementType!NumericRange ret = 0;
    import std.container;
    auto h = r.heapify!"a > b"();
    while (!h.empty())
    {
        ret += h.front;
        h.removeFront();
        if (h.empty())
            break;
        if (ret <= h.front)
            continue;
        h.insert(ret);
        ret = h.front;
        h.removeFront();
    }
    return ret;
}

auto sumSquareDiff(F)(F[] arr)
{
    F[] sd;
    foreach (i, x; arr)
        foreach (j, y; arr)
        {
            if (i == j)
                break;
            sd ~= (x-y) * (x-y);
        }
    return sd.sumFromSmallerPositive();
}

auto stdevPreciseSlow(NumericRange)(NumericRange r)
{
    auto ss = sumSquareDiff(r.array());
    immutable n = r.length;
    return (ss / (n * (n-1))).sqrt();
}
auto stdevSuperlinear(NumericRange)(NumericRange r)
{
    real sum = r.reduce!((a, b) => a + b)();
    immutable ulong n = r.length;
    immutable average = sum / n;
    sum = r.map!(a => (a - average) * (a - average))().sumFromSmallerPositive();
    return (sum / (n - 1)).sqrt();
}
/// Integral f by each P in Ps and return the unbiased variance.
auto integrationStdevPreciseSlow(alias f, PointSetTypeRange)(PointSetTypeRange Ps)
{
    return Ps.map!(bintegral!(f, ElementType!PointSetTypeRange))().stdevPreciseSlow();
}
/// ditto
auto integrationStdevSuperlinear(alias f, PointSetTypeRange)(PointSetTypeRange Ps)
{
    return Ps.map!(bintegral!(f, ElementType!PointSetTypeRange))().stdevLinear();
}
/// Given a range of numbers, return the square root of the mean square of its elements.
auto squareRootMeanSquare(NumericRange)(NumericRange r)
{
    return (r.map!(a => a * a).sumFromSmallerPositive() / r.length).sqrt();
}
