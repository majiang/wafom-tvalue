module visualize;

import std.stdio;
import lib.pointset;

void main()
{
	foreach (line; stdin.byLine())
	{
		auto P = line.fromString!ubyte();
		foreach (X; P)
		{
			foreach (x; X)
				x.write(",");
			writeln();
		}
		writeln();
	}
}
