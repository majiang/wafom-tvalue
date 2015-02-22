module main;

import lib.wafom, lib.pointsettype;
import std.stdio, std.string, std.array;
import std.math : lg = log2;

void main()
{
	foreach (line; stdin.byLine)
		"%s,%.15e".writefln(line.strip.split(',')[0], line.fromString!uint.biwafom.lg);
}
