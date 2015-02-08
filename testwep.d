import std.stdio, std.conv, std.datetime, std.bigint;
import lib.wafom, lib.pointsettype;

alias generateP = randomPointSet!uint;

void main(string[] args)
{
	if (args.length < 4)
		return stderr.writeln("testwep s m count");
	auto
		s = DimensionR(args[1].to!size_t),
		m = DimensionF2(args[2].to!size_t),
		count = args[3].to!size_t;
	StopWatch sw;
	sw.start();
	foreach (c; 0..count)
		auto P = generateP(Precision(32), s, m);
	sw.stop();
	stdout.writefln("empty loop: %d msec", sw.peek().msecs);
	sw.reset();
	sw.start();
	foreach (c; 0..count)
	{
		auto P = generateP(Precision(32), s, m);
		P.yoshiki_wep!BigInt();
	}
	sw.stop();
	stdout.writefln("Use BigInt: %d msec", sw.peek().msecs);
	sw.reset();
	sw.start();
	foreach (c; 0..count)
	{
		auto P = generateP(Precision(32), s, m);
		P.yoshiki_wep!long();
	}
	sw.stop();
	stdout.writefln("Use long: %d msec", sw.peek().msecs);
}
