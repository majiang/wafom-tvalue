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
    writer.formattedWrite("%.20e", P.integrationError!(
        oscillatoryC!(
        0.9,  // u; random
        [0.6, 0.8, 1.0, 1.2] // a; sum = 0.9 * dimR
    ))());
    // 10-dimensional: a.sum = 9.0, 7.25, 1.85, 7.03, 2.04, 4.30 for each function
    ret ~= writer.data;
    writer = appender!string();
    writer.formattedWrite("%.20e", P.integrationError!(
        prodpeakC!(
        [0.2, 0.4, 0.6, 0.8],  // u; random
        [0.9, 0.8, 0.7, 0.5] // a; sum = 7.25 * dimR
    ))());
    // 10-dimensional: a.sum = 9.0, 7.25, 1.85, 7.03, 2.04, 4.30 for each function
    ret ~= writer.data;
    writer = appender!string();
    writer.formattedWrite("%.20e", P.integrationError!(
        cornpeakC!(
        [0.4, 0.2, 0.1, 0.04] // a; sum = 1.85 * dimR
    ))());
    // 10-dimensional: a.sum = 9.0, 7.25, 1.85, 7.03, 2.04, 4.30 for each function
    ret ~= writer.data;
    writer = appender!string();
    writer.formattedWrite("%.20e", P.integrationError!(
        gaussianC!(
        [0.2, 0.1, 0.06, 0.052],  // a; sum = 7.03 * dimR
        [0.2, 0.4, 0.6, 0.8] // u; random
    ))());
    // 10-dimensional: a.sum = 9.0, 7.25, 1.85, 7.03, 2.04, 4.30 for each function
    ret ~= writer.data;
    writer = appender!string();
    writer.formattedWrite("%.20e", P.integrationError!(
        continuousC!(
        [0.4, 0.2, 0.1, 0.116], // a; sum = 2.04 * dimR
        [0.2, 0.4, 0.6, 0.8] // u; random
    ))());
    // 10-dimensional: a.sum = 9.0, 7.25, 1.85, 7.03, 2.04, 4.30 for each function
    ret ~= writer.data;
    writer = appender!string();
    writer.formattedWrite("%.20e", P.integrationError!(
        discontiC!(
        [0.7, 0.4, 0.2, 0.42],  // a; sum = 4.30 * dimR
        [0.2, 0.4, 0.6, 0.8] // u; random
    ))());
    // 10-dimensional: a.sum = 9.0, 7.25, 1.85, 7.03, 2.04, 4.30 for each function
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
