module optimize.filter;

import std.container : heapify;
import std.typecons : Tuple, tuple;
alias Tuple!(real, string) D;

void main(string[] args)
{
	import std.stdio;
	import std.string : strip;
	import std.array : split;
	import std.conv : to;
	stderr.writeln("read from stdin and write the best [capacity] to stdout.");
	stderr.writeln("The smaller third column is the better.");
	if (args.length < 2)
		return stderr.writeln("filter-good capacity");
	auto capacity = args[1].to!size_t();
	auto heap = (new D[capacity]).heapify(0);
	foreach (line; stdin.byLine())
	{
		string value = line.dup;
		real key = line.strip().split(",")[2].to!real();
		heap.conditionalInsert(D(key, value));
	}
	foreach (t; heap.release())
		t[1].writeln();
}
