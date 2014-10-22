module genzerr;
import lib.pointsettype, lib.integral, lib.testfunction;
import std.stdio, std.path;
import std.functional : pipe;
import std.string : strip;

auto idTuple(alias f, T)(T x)
{
	import std.typecons;
	return tuple(x, f(x));
}

void main(string[] args)
{
	if (args.length < 3)
		return stderr.writeln("genz-err pointsetfile testfuncfile [outputfile]");
	auto testfunctions =
		File(args[2]).byLine().map!(pipe!(to!string, stringToGenz!real))().array();

	string outfilename = args[1].absolutePath().stripExtension() ~ '-' ~ args[2].baseName(".ssv") ~ ".csv";
	if (args.length < 4)
		stderr.writeln("default output location is used.");
	else
		outfilename = args[3];
	stderr.writefln("output written to %s", outfilename);
	auto outfile = File(outfilename, "w");
	foreach (p; File(args[1]).byLine().map!(pipe!(to!string, strip, idTuple!(fromString!uint, string)))())
	{
		outfile.write(p[0]);
		foreach (f; testfunctions)
			outfile.writef(",%.15e", p[1].signedIntegrationError(f));
		outfile.writeln();
	}
}
