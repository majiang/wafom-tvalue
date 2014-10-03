import ui.input, lib.pointsettype;
import std.stdio;
import std.math : lg = log2;

version (Yoshiki)
{
	auto weight_description = "Yoshiki weight";
	version (RMS)
		import lib.wafom : criterion = bipmswafom;
	else
		import lib.wafom : criterion = bipwafom;
}
else
{
	auto weight_description = "Matsumoto-Saito-Matoba original (Dick weight)";
	version (RMS)
		import lib.wafom : criterion = bimswafom;
	else
		import lib.wafom : criterion = biwafom;
}
version (RMS)
	auto exponent_description = "RMS";
else
	auto exponent_description = "non-RMS";

void main()
{
	stderr.writefln(
"Computes Walsh Figure of Merit:
    %s
    %s
input:
    http://www.ms.u-tokyo.ac.jp/~ohori/format.html
output:
    lg(W(P)) on each line", weight_description, exponent_description);
	foreach (P; getDigitalNets!uint())
		"%.15f".writefln(P.criterion().lg());
}
