import std.stdio : stderr;
import std.exception : enforce;
import std.typecons : Tuple, tuple;
import std.functional : pipe;

import lib.wafom : bipmswafom;
import lib.pointset : ShiftedBasisPoints, randomVectors, randomVector, randomBits, nonshiftedRandomBasisPoints;

alias uint U;
alias ShiftedBasisPoints!U PointSet;
alias nonshiftedRandomBasisPoints!U randomDigitalNet;


enum precision = 30;
enum default_start_temperature = 1.0;
enum default_end_temperature = 0.01;
enum min_iteration = 0x10;
enum max_iteration = 0x1000;

/** Generate good point sets.

Read the range of dimB, the set of dimR from arguments.
Do simulated annealing for each (dimB, dimR) and yield the best one.
Precision and temperature condition are hard-coded.
*/
void main(string[] args)
{
import std.conv : to;
import std.algorithm : map;
import std.stdio;

	if (args.length < 4)
	{
		"usage: annealing m_min m_max s...".writeln();
		return;
	}
	auto
		dimensionF2min = args[1].to!size_t(),
		dimensionF2max = args[2].to!size_t(),
		dimensionRs = args[3..$].map!(to!size_t)();
	foreach (dimensionR; dimensionRs)
		foreach (dimensionF2; dimensionF2min..dimensionF2max+1)
			"%s,%.15e".writefln(dimensionR.do_anneal(dimensionF2).expand);
}

/** bisect recursive version of std.algorithm.reduce.

Only accept one alias for reducing function and a nonempty slicable (typically an array).
The purpose of rewriting is to reduce an array whose element is a struct with immutable member.
*/
template reduce(alias f)
{
import std.array : empty;

	auto reduce(T)(T arg)
	{
	import std.exception : enforce;
		enforce(!arg.empty, "Cannot reduce an empty range");
		if (arg.length == 1)
			return arg[0];
		return f(reduce(arg[0 .. $/2]), reduce(arg[$/2 .. $]));
	}
}

/** Repeat simulated annealing, varying the number of iteration.

The temperature condition is "default_start_temperature -> default_end_temperature".
Vary the number of iteration from min_iteration to max_iteration.

Returns:
The best state found in any of the iteration.
*/
auto do_anneal(in size_t dimensionR, in size_t dimensionF2)
{
	auto iteration = min_iteration;
	Tuple!(PointSet, real)[] results;
	while (iteration <= max_iteration)
	{
		results ~= iteration.anneal!
			(bipmswafom,
			(() => iteration.simple_rs!
				(bipmswafom,
					(() => randomDigitalNet(precision, dimensionR, dimensionF2)))[0]
				),
			neighborPointSet!(1, EuclideanRandomVectors), KirkpatrickAcceptance)
			(default_start_temperature, default_end_temperature);
		iteration <<= 1;
	}
	return results.reduce!((P, Q) => P[1] < Q[1] ? P : Q)();
}

version (none)
auto anneal(
	alias criterion, alias neighbor, alias acceptance, T)(
	T start_state, double start_energy,
	size_t iteration,
	double start_temperature, double end_temperature)
{
import std.math;

	auto temperature = start_temperature;
	immutable cooling = (end_temperature / start_temperature) ^^ (1.0 / iteration);
	auto state = [tuple(start_state, start_energy)];
	foreach (i; 0..iteration)
	{
		state ~= neighbor(state[$-1][0]).S!(tuple, criterion);
		if (acceptance(temperature, state[$-2][1], state[$-1][1])) // accept
			if (state[$-1][1] <= state[0][1]) // best
				state = state[$-1 .. $];
			else
				state = state[0 .. 1] ~ state[$-1 .. $];
		else
			state = state[0 .. $-1];
		temperature *= cooling;
	}
	return state[0];
}

double get_cooling(size_t iteration, double ratio)
{
import std.math;
	return ratio ^^ (1.0 ^^ (iteration - 1));
}

/** Simple random search.

Params:
criterion = criterion for random search (smaller is better).
generate = uniform random choice function in the search space.
iteration = the number of iteration (call of generate and criterion)

Returns:
the two-element tuple of best object and its quality.
*/
auto simple_rs(alias criterion, alias generate)(size_t iteration)
{
	auto state = [generate().S!(tuple, criterion)];
	foreach (i; 1..iteration)
	{
		state ~= generate().S!(tuple, criterion);
		if (state[$-1][1] < state[$-2][1])
			state = state[$-1 .. $];
		else
			state = state[0 .. $-1];
	}
	return state[0];
}

/** Simulated annealing.

Params:
criterion = criterion (smaller is better).
initial = initial state (only an object in search space) generator.
neighbor = neighborhood generator.
acceptance = acceptance probability.
iteration = the number of basic iterations in the process.
start_temperature = the temperature at the beginning of the process.
end_temperature = the temperature at the end of the process.

Returns:
the two-element tupleof best object and its quality.

Algorithm:
Cooling schedule is exponential, i.e., the cooling rate c is calculated from the arguments so that after each basic iteration the temperature is multiplied by c and at the end of the process it holds that temperature = end_temperature.
Strictly speaking this algorithm is not the pure SA, but remembers the best state.
*/
auto anneal(
	alias criterion, alias initial, alias neighbor, alias acceptance)(
	size_t iteration, double start_temperature, double end_temperature)
{
	auto temperature = start_temperature;
	immutable cooling = iteration.get_cooling(end_temperature / start_temperature);

	auto state = [initial().S!(tuple, criterion)];
	foreach (i; 0..iteration)
	{
		state ~= neighbor(state[$-1][0]).S!(tuple, criterion);
		if (acceptance(temperature, state[$-2][1], state[$-1][1])) // accept
			if (state[$-1][1] <= state[0][1]) // best
				state = state[$-1 .. $];
			else
				state = state[0 .. 1] ~ state[$-1 .. $];
		else
			state = state[0 .. $-1];
		temperature *= cooling;
	}
	return state[0];
}

///
auto neighborPointSet(size_t distance, alias rvectors = randomVectors)(PointSet P)
{
	auto vectors = rvectors!U(P.precision, P.dimensionR, distance);
	auto xor_at = P.dimensionF2.randomVector!U(distance);
	auto basis = P.basis.ddup();
	foreach (i, vector; vectors)
		foreach (j, ref b; basis)
			if (xor_at[i] >> j & 1)
				foreach (k, v; vector)
					b[k] ^= v;
	return PointSet(basis, P.precision);
}

auto EuclideanRandomVectors(U)(size_t precision, size_t dimensionR, size_t distance)
{
import std.random : uniform;
	U[][] ret;
	foreach (i; 0..distance)
	{
		ret ~= new U[dimensionR];
		U b = 1;
		foreach (j; 0..precision)
		{
			ret[$-1][0.uniform(dimensionR)] |= j;
			j <<= 1;
		}
	}
	return ret;
}

import std.math;
///
bool KirkpatrickAcceptance(double temperature, double current, double next)
{
	//if (next < current)
	//	return true;
import std.random : uniform;
	return uniform(0.0, 1.0) < (current / next) ^^ (1 / temperature);
}


/// Substitution combinator S := \xyz. xz(yz)
auto S(alias x, alias y, T)(T z)
{
	return x(z, y(z));
}

///
T[][] ddup(T)(in T[][] arr)
{
    T[][] ret;
    foreach (line; arr)
        ret ~= line.dup;
    return ret;
}
