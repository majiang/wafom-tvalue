import std.array, std.range, std.algorithm, std.string;
import std.stdio;

auto read(T)(T x, string col)
{
	return "'%s' using 1:%s title '%s'".format(x[0], col, x[1]);
}

void main(string[] args)
{
	args.popFront(); // remove the name of this script
	"set terminal png".writefln();
	"set output '%s'".writefln(args[0]);
	"set datafile separator ','".writefln();
	args.popFront(); // remove the output file name
	("plot " ~ args[1..$].chunks(2).map!(x => read(x, args[0]))().join(", ")).writeln();
}
