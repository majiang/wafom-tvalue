module lib.inttest;

import lib.testfunction : genzFactory;
import lib.integral : signedIntegrationError;
import lib.pointsettype;
import std.functional : memoize;

alias memoizedFactory(F) = memoize!(genzFactory!F);

auto genzIntegrationError(F, R)(R P, size_t index, F difficulty)
{
	return P.signedIntegrationError(memoizedFactory!F(index, P.dimensionR, difficulty));
}

auto genzIntegrationErrorsDigitalShift(F, R)(R P, size_t index, size_t numDS, F difficulty = 1)
{
	F[] ret;
	foreach (Q; P.randomlyShifted(numDS))
		ret ~= Q.genzIntegrationError(index, difficulty);
	return ret;
}

void main(string[] args)
{
	import std.algorithm : map;
	import std.functional : pipe;
	import std.stdio;
	import std.conv : to;
	import std.math : abs, lg = log2;
	if (args.length != 4)
		return "inttest genzIndex numDigitalShift difficulty".writeln();
	immutable index = args[1].to!size_t();
	immutable numDS = args[2].to!size_t();
	immutable difficulty = args[3].to!real();
	import ui.input;
	if (numDS)
		foreach (P; getDigitalNets!uint())
			"%s,%(%.15f,%)".writefln(P, P.genzIntegrationErrorsDigitalShift!real(index, numDS, difficulty).map!(pipe!(abs, lg))());
	else
		foreach (P; getDigitalNets!uint())
			"%s,%.15f".writefln(P, P.genzIntegrationError!real(index, difficulty).abs().lg());
}
