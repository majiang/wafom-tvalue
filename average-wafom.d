module main;

import std.stdio;
import std.conv : to;
import std.math : lg = log2;
import lib.approximate : WAFOM;

struct EmptyDigitalNet(U)
{
	U[] front;
	static immutable size_t dimensionF2 = 0;
	immutable size_t dimensionR, precision;
	@property bool empty() const
	{
		return _empty;
	}
	void popFront()
	{
		_empty = true;
	}
	this (size_t s, size_t n)
	{
		dimensionR = s;
		front.length = s;
		precision = n;
	}
	static immutable bool bisectable = false;
	typeof (this)[2] bisect()
	{
		assert (false);
	}
private:
	bool _empty;
}

void main(string[] args)
{
	EmptyDigitalNet!uint
		(args[1].to!size_t(), args[2].to!size_t()).
		WAFOM(args[3].to!real()).lg().writeln();
}
