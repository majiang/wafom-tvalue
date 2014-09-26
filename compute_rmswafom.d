import ui.input, lib.wafom, lib.pointsettype;
import std.stdio;
import std.math : lg = log2;

void main()
{
	stderr.writeln(
"Computes Walsh Figure of Merit:
    Matsumoto-Saito-Matoba original (Dick weight)
    Non-RMS
input:
    http://www.ms.u-tokyo.ac.jp/~ohori/format.html
output:
    lg(W(P)) on each line");
	foreach (P; getDigitalNets!uint())
		"%.15f".writefln(P.bimswafom().lg());
}
