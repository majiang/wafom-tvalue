import std.stdio;
alias uint U;
import lib.wafom : bipwafom;
import genz : 
    oscillatoryC,
    prodpeakC,
    cornpeakC,
    gaussianC,
    continuousC,
    discontiC;



string integrationErrors(R)(R P)
{
    import lib.integration_error : integrationError;
    import std.array : appender, join;
    import std.format : formattedWrite;
    string[] ret;
    auto writer = appender!string();
    writer.formattedWrite("%.20e", P.integrationError!(oscillatoryC!(
    0.9,  // u; random
    [0.6, 0.8, 1.0, 1.2] // a; sum = 0.9 * dimR
    ))());
    ret ~= writer.data;
    writer = appender!string();
    writer.formattedWrite("%.20e", P.bipwafom());
    ret ~= writer.data;
    return ret.join(",");
}

void main()
{
    import ui.output : writeWith;
    import ui.input : getDigitalNets;
    alias getDigitalNets!U getDN;
    foreach (P; getDN())
        P.writeWith!integrationErrors();
}
