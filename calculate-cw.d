module wafom.calculate;

import lib.approximate, lib.pointsettype;
import ui.input;
alias U = uint;

void main(string[] args)
{
	import std.stdio;
	import std.conv : to;
	import std.math : lg = log2;
	if (args.length < 2)
		return stderr.writeln("specify c-value.");
	auto c = args[1].to!real();
	foreach (P; getDigitalNets!U)
		"%s,%.15f,%.15f".writefln(P, c, P.WAFOM(c).lg());
}
