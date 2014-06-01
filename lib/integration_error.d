module lib.integration_error;

import lib.integral;
//import testfunction;

debug import std.stdio;
import std.range : ElementType, hasLength;
import std.traits : ReturnType;
import std.algorithm : map, reduce;

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
auto sumSquareDiff(F)(F[] arr)
{
    import std.container;
    F[] sd;
    foreach (i, x; arr)
        foreach (j, y; arr)
        {
            if (i == j)
                break;
            sd ~= (x-y) * (x-y);
        }
    F ret = 0;
    auto h = sd.heapify!"a > b"();
    while (!h.empty())
    {
        debug "%.5f + %.5f = %.5f".writefln(ret, h.front, h.front+ret);
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
/// Integral f by each P in Ps and return the unbiased variance.
auto integrationStdevPreciseSlow(alias f, PointSetTypeRange)(PointSetTypeRange Ps)
{
    auto I = Ps.map!(bintegral!(f, ElementType!PointSetTypeRange))().array();
    auto ss = sumSquareDiff(I);
    immutable n = I.length;
    return (ss / (n * (n-1))).sqrt();
}
/// ditto
auto integrationStdev(alias f, PointSetTypeRange)(PointSetTypeRange Ps)
{
    auto I = Ps.map!(bintegral!(f, ElementType!PointSetTypeRange))();
    auto sum = I.reduce!((a, b) => a + b)();
    static if (hasLength!(typeof (I)))
        immutable ulong count = I.length;
    else
        immutable ulong count = I.walkLength();
    immutable average = sum / count;
    sum = (cast(ElementType!I)0).reduce((a, b) => a + (b - average) * (b - average))();
     return (sum / (count - 1)).sqrt();
}
/// Given a range of numbers, return the square root of the mean square of its elements.
auto squareRootMeanSquare(NumericRange)(NumericRange r)
{
    import std.math : sqrt;
    static if (is (ElementType!NumericRange == float))
        float sum = 0;
    else static if (is (ElementType!NumericRange == real))
        real sum = 0;
    else
        real sum = 0;
    static if (hasLength!NumericRange)
        ulong count = r.length;
    else
        ulong count;
    foreach (e; r)
    {
        sum += e * e;
        static if (hasLength!NumericRange)
            ++count;
    }
    return (sum / count).sqrt();
}

auto shifts(PointSetType)(in size_t precision, in size_t dimensionR, in size_t count)
{
    import lib.pointset : randomVectors;
    return randomVectors!(PointSetType.ComponentType)(precision, dimensionR, count);
}

/** Shift a point set by each element of shifts.

Params:
P = point set.
shifts = shifts
*/
auto shifteds(PointSetType, ShiftsType)(PointSetType P, ShiftsType shifts)
{
    return shifts.map!(x => P.shifted(x));
}
