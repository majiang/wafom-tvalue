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

import std.stdio;

void main(string[] args)
{
	import std.conv : to;
	import std.algorithm : map;
	if (args.length < 4)
		return usage();
	auto
		dimensionF2min = args[1].to!size_t(),
		dimensionF2max = args[2].to!size_t(),
		dimensionRs = args[3..$].map!(to!size_t)();
	foreach (dimensionR; dimensionRs)
		foreach (dimensionF2; dimensionF2min .. dimensionF2max+1)
		{
			auto P = defaultSobols!U(dimensionR, precision, dimensionF2);
			auto w = P.bipmswafom();
			"%s,%.15e".writefln(P, w);
//			auto Q = defaultSobols!U(dimensionR, dimensionF2, dimensionF2);
//			"%s,%.15e".writefln(Q, w);
		}
}

void usage()
{
	"usage: refer_sobol m_min m_max s...".writeln();
}

auto get_sobol_w2(in size_t dimensionR, in size_t dimensionF2)
{
	return defaultSobols!U(dimensionR, precision, dimensionF2).bipmswafom();
}
