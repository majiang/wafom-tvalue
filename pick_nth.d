import std.stdio, std.conv, std.algorithm, std.array;

void main(string[] args)
{
    auto p = args[1..$].map!(to!size_t)();
    size_t i;
    foreach (line; stdin.byLine(KeepTerminator.yes))
    {
        if (++i == p.front)
        {
            line.write();
            p.popFront();
        }
    }
}
