main functionality of the files:

pointset.d
	struct BasisPoints : InputRange
		this(ulong[][] basis, size_t precision) // dimension = basis.length; length = 1 << basis[0].length
		this(BasisPoints other)
	BasisPoints randomPoints(size_t dimension, size_t precision, size_t lg_length)
sobol.d
	BasisPoints defaultSobols(size_t dimension, size_t precision, size_t lg_length)
	BasisPoints sobols(ulong[][] direction_numbers)
integral.d
	double integral(alias f)(BasisPoints P)
	void integral_sobol(BasisPoints P) // 各パラメータに対して、Asian Option の価格を、与えられた BasisPoints P 上の QMC 積分によって計算する
	double integral_ao(double r, double T, int S, double P0, double sigma, double K, BasisPoints P) // integrand にパラメータと P の各点を渡して積分する
wafom.d
	double wafom(R)(R P)
tvalue.d
	ulong tvalue1(R)(R P)
	ulong tvalue2(R)(R P)
asianoption.d
	double integrand(double r, double T, int S, double P0, double sigma, double K, double[] x)
	double ctintegrand(double r, double T, int S, double P0, double sigma, double K)(double[] x)
