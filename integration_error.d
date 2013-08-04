module integration_error;

import integral;
import testfunction;

import std.range : ElementType, hasLength;
import std.traits : ReturnType;
import std.algorithm : map;

/** Calculate integration error of a function by a point set.

Params:
tf = an alias. the integrand is given as tf.f, with its true integrated value tf.I.
P = a point set.
*/
double integrationError(alias tf, PointSetType)(PointSetType P)
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

auto squareRootMeanSquare(NumericRange)(NumericRange r)
{
    import std.math : sqrt;
    static if (is (ElementType!NumericRange == float))
        float sum = 0;
    else static if (is (ElementType!NumericRange == real))
        real sum = 0;
    else
        double sum = 0;
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
    import pointset : randomVectors;
    return randomVectors!(PointSetType.ComponentType)(precision, dimensionR, count);
}

/** Shift a point set by each element of shifts.

Params:
P = point set.
shifts = shifts
*/
auto shifteds(PointSetType)(PointSetType P, ReturnType!(shifts!PointSetType) shifts)
{
    return shifts.map!(x => P.shifted(x));
}
