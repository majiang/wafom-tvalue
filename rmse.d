import lib.wafom, lib.pointsettype, lib.integration_error, ui.input, tf;

alias uint U;

void WRITEERRORS(R, F...)(R P)
{
    import std.stdio;
    "%s,%.15e".writef(P, P.bipmswafom()); // unnecessary recalculate
    auto Qs = P.randomShiftsFor(1024).map!(sigma => P + sigma)();
    foreach (f; F)
        writef(",%.15e", Qs.integrationStdevPreciseSlowNoncentering!(f.f)());
    writeln();
}

void main()
{
    foreach (P; getDigitalNets!U())
        P.WRITEERRORS!(typeof (P),
            x6,
            exponential!(2.0/3.0, S),
            exponential!(3.0/2.0, S),
            coss100,
            gauss100,
            pp100,
            conti,
            disco
        );
}
