main functionality of the files:

pointset.d
    struct ShiftedBasisPoints(T) if (isUnsigned!T) : InputRange
        this(in T[][] basis, in size_t precision, in T[] shift) // dimensionF2 = basis.length; dimensionR = shift.length;
        this(in T[][] basis, in size_t precision) // dimensionF2 = basis.length; dimensionR = basis[0].length;
    ShiftedBasisPoints randomBasisPoints(T)(size_t precision, size_t dimensionR, size_t dimensionF2, Flag!"shift" shift) if (isUnsigned!T)
integral.d
    double integral(alias f, R)(R P)
    double bintegral(alias f, R)(R P) if (isBisectable!R)
wafom.d
    template biwafom(R) // Dick WAFOM
    template prwafom(R) // precise Dick WAFOM
    template bimswafom(R) // root mean square Dick WAFOM
    template binrtwafom(R) // NRT WAFOM
    template bimsnrtwafom(R) // root mean square NRT WAFOM
    string dick_weight_enumerator_polynomial(R)(R P) // Dick weight enumerator polynomial
tvalue.d
    ulong tvalue1(R)(R P)
    ulong tvalue2(R)(R P)
asianoption.d
    double integrand(double r, double T, int S, double P0, double sigma, double K, double[] x)
    double ctintegrand(double r, double T, int S, double P0, double sigma, double K)(double[] x)
    double default_integrand(double[] x)
randomsearch.d
    R minimum(alias score, alias generate, R, Param...)(size_t count, Param params)
walsh.d
    struct Oresen
        double integral()
        double opCall(double x)
    double[2][] walsh_function(size_t k, size_t v)
    int wal(size_t k, size_t v, size_t x)
    Oresen product_walsh(Oresen f, size_t k)
