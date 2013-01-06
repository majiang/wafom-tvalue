module finitefield;

alias F!(2, 1) B;

struct F(ulong p, ulong e) if (supported!(p, e)())
{
    ulong v;
    F opBinary(string op)(F other)
    {
        static if (op == "+")
        {
            auto a = this.v, b = other.v;
            ulong c;
            foreach (i; 0..e)
            {
                c += ((a + b) % p) * (p ^^ i);
                a /= p;
                b /= p;
            }
            return F(c);
        }
        static if (op == "-")
        {
            auto a = this.v, b = other.v;
            ulong c;
            foreach (i; 0..e)
            {
                c += (((a - b) % p + p) % p) * (p ^^ i);
                a /= p;
                b /= p;
            }
            return F(c);
        }
        static if (op == "*")
        {
            static assert (e == 1);
            return F((this.v * other.v) % p);
        }
        static if (op == "/")
        {
            assert (other.v);
            return this * other ^^ (p ^^ e - 2);
        }
    }
    F opBinary(string op)(ulong exponent)
    {
        static if (op == "^^")
        {
            auto ret = F(1);
            auto square = this;
            while (exponent)
            {
                if (exponent & 1)
                {
                    ret *= square;
                }
                exponent >>= 1;
                square *= square;
            }
            return ret;
        }
    }
}

bool supported(ulong p, ulong e)()
{
    if (!(primality(p) && small_enough(p, e)))
    {
        return false;
    }
    return p == 2 && e == 1;
}

bool small_enough(ulong p, ulong e)
{
    if (p == 2)
        return e < 32;
    ulong m = ulong.max;
    foreach (i; 0..e)
    {
        m /= p;
    }
    return cast(bool)m;
}

bool primality(ulong p)
{
    if (p < 4) return 2 <= p;
    if ((p & 1) == 0) return false;
    ulong i = 3;
    while (i * i <= p)
    {
        if (p % i == 0)
            return false;
        i += 2;
    }
    return true;
}

unittest
{
    foreach (p; [2, 3, 5, 7, 11]) assert (primality(p));
    foreach (c; [0, 1, 4, 6, 8, 9, 10]) assert (!primality(c));
}
