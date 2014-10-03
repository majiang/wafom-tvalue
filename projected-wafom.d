module projected_wafom;

import ui.input, lib.pointsettype;
import std.stdio;
import std.math : lg = log2;
import std.range : iota;
import std.algorithm : map;

version (Yoshiki)
{
	auto weight_description = "Yoshiki weight";
	version (RMS)
		import lib.wafom : criterion = bipmswafom;
	else
		import lib.wafom : criterion = bipwafom;
}
else
{
	auto weight_description = "Matsumoto-Saito-Matoba original (Dick weight)";
	version (RMS)
		import lib.wafom : criterion = bimswafom;
	else
		import lib.wafom : criterion = biwafom;
}
version (RMS)
	auto exponent_description = "RMS";
else
	auto exponent_description = "non-RMS";

auto projection(Size)(Size p)
{
	size_t[] ret;
	size_t i;
	while (p)
	{
		if (p & 1)
			ret ~= i;
		p >>= 1;
		i += 1;
	}
	return ret;
}

void main()
{
	stderr.writefln(
"Computes Walsh Figure of Merit of all projection:
    %s
    %s
input:
    http://www.ms.u-tokyo.ac.jp/~ohori/format.html
output:
    lg(W(P)) on each line", weight_description, exponent_description);
	foreach (P; getDigitalNets!uint())
	{
    	auto mask = (1 << P.dimensionR) - 1;
		"%(%.15f,%|%)".writef(
			iota(1, 1 << P.dimensionR).
			map!(p => P.projectionTo(p.projection()).criterion().lg())());
		"%(%.15f,%)".writefln(
			iota(1, mask).
			map!(p => (P.projectionTo(p.projection()) & P.projectionTo((mask - p).projection())).criterion().lg())());
	}
}
