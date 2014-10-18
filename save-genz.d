module savegenz;
import lib.testfunction;
import std.stdio;
import std.conv : to;
import std.string : format;

version (dynamic)
{}
else
	static assert (false);

void main(string[] args)
{
	if (args.length < 4)
		return stderr.writeln("save-genz index dimR numSamples [coefficient]");
	if (5 < args.length)
		stderr.writeln("warning: ignoring some parameters");
	immutable
		index = args[1].to!size_t(),
		dimR = args[2].to!size_t(),
		numSamples = args[3].to!size_t();
	real coefficient = 1;
	if (args.length == 5)
		coefficient = args[4].to!real();
	auto output = File("genz%d-s%02d-c%.0f-%d.ssv".format(index, dimR, coefficient, numSamples), "w");
	foreach (i; 0..numSamples)
		output.writeln(genzFactory!real(index, dimR, coefficient).toString());
}
