module hellekalek;

import func.testfunction;
import std.stdio;
import lib.integral;
import lib.pointsettype;
import std.math : lg = log2, abs;
import std.string;

void main()
{
	alias h = Hellekalek!(1.7, 1.9, 2.1, 2.3);
	enum s = 4;
	foreach (line; stdin.byLine())
	{
		auto P = line.fromString!uint();
		if (P.dimensionR != s)
			continue;
		"%s,%.15f".writefln(line.strip(), P.bintegral!(h.f).abs().lg());
	}
}
