import std.stdio;
import std.functional : pipe;
import std.string : strip;
import std.array : split;
import std.algorithm : map;
import std.conv : to;
import std.array : array;
import core.bitop : bitswap;

void main()
{
	size_t precision = 32, dimensionF2, dimensionR;
	uint[][] vectors;
	foreach (line; stdin.byLine().map!(pipe!(strip, split, map!(pipe!(to!uint, bitswap)), array)))
	{
		if (dimensionR == 0)
		{
			dimensionF2 = line.length;
			vectors = line.map!(a => [a])().array();
		}
		else
		{
			assert (dimensionF2 == line.length);
			foreach (i, ref vector; vectors)
				vector ~= line[i];
		}
		dimensionR += 1;
	}
	"%d %d %d%(%( %d%)%)".writefln(precision, dimensionF2, dimensionR, vectors);
}
