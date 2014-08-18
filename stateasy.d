import std.string, std.array, std.conv, std.algorithm, std.math, std.stdio;
alias lg = log2;

auto manipulateInformation(in char[] line)
{
	auto buf = line.split(',');
	immutable size_t m = buf[0].split()[1].to!size_t();
	real[] lgerr;
	foreach (f; buf[1..$])
		lgerr ~= f.to!real();
	struct Ret
	{
		immutable size_t m;
		immutable real rms, sup;
		string toString()
		{
			return "%d,%.15f,%.15f".format(m, rms, sup);
		}
	}
	immutable sup = lgerr.reduce!max();
	real ss = 0;
	foreach (e; lgerr)
		ss += e.exp2().pow(2);
	immutable rms = (ss / lgerr.length).lg() / 2;
	return Ret(m, rms, sup);
}

void main()
{
	foreach (line; stdin.byLine())
		line.strip().manipulateInformation().writeln();
}
