import bf;

auto bipmswafom(PointSet P)
{
	return P.standardDickYoshikiWafom(-2);
}

import lib.pointset : ShiftedBasisPoints, randomVectors, randomVector, randomBits, defaultSobols, nonshiftedRandomBasisPoints;

// lib.pointset : defaultSobols(U)(dimR, precision, dimF);

import std.exception : enforce;

alias uint U;
alias ShiftedBasisPoints!U PointSet;
alias nonshiftedRandomBasisPoints!U randomDigitalNet;

import std.stdio : stderr;
import std.typecons : Tuple, tuple;

enum precision = 32;
enum default_start_temperature = 1.0;
enum default_end_temperature = 0.01;
enum min_iteration = 0x10;
enum max_iteration = 0x1000;

void main(string[] args)
{
	import std.conv : to;
	import std.algorithm : map;
	import std.stdio;

	if (args.length < 4)
	{
		"usage: beat_sobol m_min m_max s...".writeln();
		return;
	}
	auto
		dimensionF2min = args[1].to!size_t(),
		dimensionF2max = args[2].to!size_t(),
		dimensionRs = args[3..$].map!(to!size_t)();
	foreach (dimensionR; dimensionRs)
		foreach (dimensionF2; dimensionF2min..dimensionF2max)
		{
			auto result = dimensionR.do_anneal(dimensionF2);
			"%s,%.15e".writefln(result[0], result[1]);
		}
}

import std.array : empty;

template reduce(alias f)
{
	auto reduce(T)(T arg)
	{
		import std.exception : enforce;
		enforce(!arg.empty, "Cannot reduce an empty range");
		if (arg.length == 1)
			return arg[0];
		return f(reduce(arg[0 .. $/2]), reduce(arg[$/2 .. $]));
	}
}

auto do_anneal(in size_t dimensionR, in size_t dimensionF2)
{
	auto iteration = min_iteration;
	Tuple!(PointSet, double)[] results;
	while (iteration <= max_iteration)
	{
		auto start = randomDigitalNet(precision, dimensionR, dimensionF2);
		results ~= start.anneal!
			(bipmswafom, neighborPointSet!1, KirkpatrickAcceptance)
			(start.bipmswafom(), iteration, default_start_temperature, default_end_temperature);
		iteration <<= 1;
	}
	return results.reduce!((P, Q) => P[1] < Q[1] ? P : Q)();
}


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

auto neighborPointSet(size_t distance)(PointSet P)
{
	auto vectors = P.precision.randomVectors!U(P.dimensionR, distance);
	auto xor_at = P.dimensionF2.randomVector!U(distance);
	auto basis = P.basis.ddup();
	foreach (i, vector; vectors)
		foreach (j, ref b; basis)
			if (xor_at[i] >> j & 1)
				foreach (k, v; vector)
					b[k] ^= v;
	return PointSet(basis, P.precision);
}

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

T[][] ddup(T)(in T[][] arr)
{
    T[][] ret;
    foreach (line; arr)
        ret ~= line.dup;
    return ret;
}

