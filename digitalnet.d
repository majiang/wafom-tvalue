module digitalnet;

debug import std.stdio;

F[m][s] matrices_times_vectors(F, ulong s, ulong m)(F[m][m][s] C, F[m] pn)
{
    F[m][s] ret;
    foreach (i; 0..s)
    {
        ret[i] = C[i].matrix_times_vector!(F, m)(pn);
    }
    return ret;
}

F[m] matrix_times_vector(F, ulong m)(F[m][m] C, F[m] pn)
{
    F[m] ret;
    foreach (i; 0..m)
    {
        ret[i] = F(0);
        foreach (j; 0..m)
        {
            ret[i] += C[i][j] * pn[j];
        }
    }
    return ret;
}

ulong[p ^^ (e * m)] DigitalNet
(
    ulong p, ulong e, ulong m, ulong s
)(
    F!(p, e)[m][m][s] C,
    F!(p, e) delegate (ulong) phi,
    ulong delegate (F!(p, e)) ihp
)
{
    auto q = p ^^ e;
    ulong[m] n;
    F!(p, e)[m] pn;
    foreach (i; 0..m)
    {
        pn[i] = phi(n[i]);
    }
    ulong[q ^^ m] ret;
    foreach (i; 0..q ^^ m)
    {
        auto y = C.matrices_times_vectors!(F!(p, e), s, m)(pn);
        foreach (j; 0..m)
        {
            ret[i] += ihp(y[j]) * q ^^ (m-1-j);
        }

        // count up in q-ary
        n[0] += 1;
        int j = 0;
        while (n[j] == q)
        {
            n[j] = 0;
            j += 1;
            n[j] += 1;
        }
        while (j)
        {
            pn[j] = phi(n[j]);
            j -= 1;
        }
        pn[0] = phi[n[0]];
    }
    return ret;
}
