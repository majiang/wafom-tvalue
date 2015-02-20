module mvnormtest;

import lib.testfunction, lib.integral, lib.pointsettype;
import ui.input;
import std.stdio, std.range, std.array, std.math;

enum size_t s = 2;

void main(string[] args)
{
	size_t c;
	if (1 < args.length)
		c = args[1].to!size_t;
	real[][] A;
	foreach (i; 0..s)
	{
		A.length += 1;
		foreach (j; 0..i)
			A[$-1] ~= 0.05;
		A[$-1] ~= 1.0;
	}
	real[] b;
	foreach (i; 0..s)
		b ~= 0;
	MultivariateNormalProbability!real[] testfuncs;
	foreach (dd; 0..20)
		testfuncs ~= new MultivariateNormalProbability!real(A, b, dd);
	//auto naive = new MultivariateNormalProbability!real(A, b, 0); // naive method
	//auto pddd = new MultivariateNormalProbability!real(A, b); // use positive definite & diagonally dominant
	import lib.approximate : WAFOM;
	enum real wc = 1;
	immutable real U = -(2 ^^ -wc) * expm1(- (2 ^^ wc));

	foreach (Pp; getDigitalNets!uint())
	{
		auto P = Pp.projectionTo(s.iota.array);
		writefln("%s,%(%.10e,%),%.10e", P, testfuncs.map!(f => P.integral(f)), P.WAFOM(wc) / U);
		foreach (i; 0..c)
		{
			auto Q = P.shiftRandomly;
			writefln("%s,%(%.10e,%),%.10e", Q, testfuncs.map!(f => Q.integral(f)), Q.WAFOM(wc) / U);
		}
	}
}
