#include <random>					// Required for Mersenne twister
#include "erf.h"					// Required for inverf
// Random variate generation
// Creates an instance of Mersenne Twister with initial seed 123
std::mt19937 generator(123);
double RanMT()
{//Generates random number on (0, 1) using Mersenne Twister
	std::uniform_real_distribution<double> dis(0.0, 1.0);
	return dis(generator);
}
bool IniciaMT()
{//Initializes seed for Mersenne Twister
    generator.seed(123);
    return true;
}
double UnifMT(const double& a, const double& b)
{//Genera uniforme entre a y b
    return (a + (b - a) * RanMT());
}
// Return A random variate from a normal distribution
double NormMT(const double& media, const double& des)
{// Marsaglia-Bray Algorithm for Normal distribution
    double uni, v, w, x, sum = 0.0;
    x = 0;
    uni = RanMT();
    if (uni <= 0.8638)
    {
        v = UnifMT(-1.0, 1.0);
        w = UnifMT(-1.0, 1.0);
        x = 2.3153508 * uni - 1 + v + w;
        return (media + des * x);
    }
    else
    {
        if (uni <= 0.9745)
        {
            v = RanMT();
            x = (3 / 2) * (v - 1 + 9.0334237 * (uni - 0.8638));
            return (media + des * x);
        }
        else
        {
            if (0.9973002 < uni)
            {
                do
                {
                    v = RanMT();
                    w = RanMT();
                    x = 9 / 2 - log(w);
                } while (x * pow(v, 2) <= 9 / 2);
                if (uni >= 0.9986501)
                {
                    x = sqrt(2 * x);
                    return (media + des * x);
                }
                else
                {
                    x = -sqrt(2 * x);
                    return (media + des * x);
                }
            }
            else //0.9745<uni<=0.9973002
            {
                do
                {
                    x = UnifMT(-3, 3);
                    uni = RanMT();
                    v = fabs(x);
                    w = 6.6313339 * (3 - v) * (3 - v);
                    sum = 0;
                    if (v < 3 / 2)
                    {
                        sum = 6.0432809 * (3 / 2 - v);
                    }
                    if (v < 1)
                    {
                        sum = sum + 13.2626678 * (3 - pow(v, 2)) - w;
                    }

                } while (uni <= 49.0024445 * exp(-pow(v, 2) / 2) - sum - w);
                return (media + des * x);
            }

        }
    }
}
// Return a random variate from a lognormal distribution
double LogNMT(const double& media, const double& des)
{// LogNormal con parámetros media y des
	return exp(NormMT(media, des));
}
// Return a random variate from a Johnson-SB distribution
double JohBMT(const double& mu, const double& sigma, const double& a, const double& b)
{// Generador de Johnson B, a es mínimo y b es multiplicador
	double y;
	y = NormMT(mu, sigma);
	y = 1 / (1 + exp(-y));
	return a + b * y;
}
// Return a random variate from a Johnson-SU distribution
double JohUMT(const double& mu, const double& sigma, const double& a, const double& b)
{// Generador de Johnson U, a es mínimo y b es multiplicador
	double y = exp(NormMT(mu, sigma));
	y = (y - (1.0 / y)) / 2;
	return a + b * y;
}
//Funciones de distribución, de densidad y verosimilitudes
//DNor			regresa la densidad de una normal
//CumNHart		regresa la FDA de la distribución normal
//JohUAc		regresa la FDA de la distribución Johnson-SU
//INormAc		regresa el inverso de la distribución acumulada normal
#define PI  3.14159265358979323846
#define RA2 1.41421356237309504880
#define UNDERSCALE (1e-140)				// For KolC function
#define OVERSCALE  (1e+140)				// For KolC function
#define BIGNUMBER  (1e+14)				// For log likelihood
double DNor(const double &x, const double &Mu, const double &StDev)
{//Calcula la densidad de una distribución normal con parámetros Mu, StDev
 //StDev debe ser positivo
	double z = (x - Mu) / StDev;
	return exp(-z * z / 2.0) / sqrt(2.0 * PI);
}
double  CumNHart(const double &y, const double &Mu, const double &Sigma) 
{//Calcula P[X <= y] para X normal con media Mu y Desv. Est. Sigma
	double x, XAbs, Temp, Temp1, Temp2;
	if (Sigma <= 0.0)
	{
		if (y < Mu) return 0.0;
		else return 1.0;
	}
	x = (y - Mu) / Sigma;
	XAbs = abs(x);
	if (XAbs > 37.0) Temp2 = 0.0;
	else
	{
		Temp1 = exp(-XAbs * XAbs / 2.0);
		if (XAbs < 7.07106781186547)
		{
			Temp = 0.0352624965998911 * XAbs + 0.700383064443688;
			Temp = Temp * XAbs + 6.37396220353165;
			Temp = Temp * XAbs + 33.912866078383;
			Temp = Temp * XAbs + 112.079291497871;
			Temp = Temp * XAbs + 221.213596169931;
			Temp = Temp * XAbs + 220.206867912376;
			Temp2 = Temp1 * Temp;
			Temp = 0.0883883476483184 * XAbs + 1.75566716318264;
			Temp = Temp * XAbs + 16.064177579207;
			Temp = Temp * XAbs + 86.7807322029461;
			Temp = Temp * XAbs + 296.564248779674;
			Temp = Temp * XAbs + 637.333633378831;
			Temp = Temp * XAbs + 793.826512519948;
			Temp = Temp * XAbs + 440.413735824752;
			Temp2 = Temp2 / Temp;
		}
		else
		{
			Temp = XAbs + 0.65;
			Temp = XAbs + 4.0 / Temp;
			Temp = XAbs + 3.0 / Temp;
			Temp = XAbs + 2.0 / Temp;
			Temp = XAbs + 1.0 / Temp;
			Temp2 = Temp1 / Temp / 2.506628274631;
		}
	}
	if (x > 0) return 1 - Temp2;
	else return Temp2;
}
double  JohUAc(const double &y, const double &MuN, const double &SigmaN) 
{//Calcula P[X <= y] para X Johnson-SU con parámetros MuN y SigmaN
	double x = log(y + sqrt(1.0 + y*y));
	return CumNHart(x, MuN, SigmaN);
}
/* --------------------------------------------------------------------------
   Funciones para calcular cdf
-----------------------------------------------------------------------------*/
double GammaAc(const double& x, const double& Alpha, const double& Beta)
{//Calcula la probabilidad acumulada de una distribución gamma
	if ((x <= 0.0) || (Alpha <= 0.0) || (Beta <= 0.0)) return 0.0;
	else {
		Gamma gama;
		return gama.gammp(Alpha, x / Beta);
	}
}
double ChiAc(const double& x, const double& GLib)
{//Calcula la probabilidad acumulada de una distribución gamma
	if ((x <= 0.0) || (GLib <= 0.0)) return 0.0;
	else return GammaAc(x, GLib / 2.0, 2.0);
}
double  KolC(double& d, int& n)
{// Algoritmo de Carvalho, regresa P[Kolmogorov-Smirnov_n <= d]
	if ((d <= 0.0) || (n <= 0)) return 0;
	int ns = 0; /* ns is number of exp shifts */
	double u, s;
	int k = (int)(n * d) + 1; /* ceiling(n * d) */
	int m = 2 * k - 1, m2 = 2 * m;
	double h = k - n * d;
	// Initialize v, w, and q 
	VecDoub v(3 * m - 2);					// `v` is at least max{2, m + m + (m - 2)} long 
	//double *q = v + m;				// q[i] = v[i + m]
	//double *w = q + m;				// 'w' is m - 2 long w[i] = v[i + m2]
	for (int i = 0; i < m; i++) v[i + m] = 0;
	v[k + m - 1] = u = s = 1.0;			// q = e_k 
	for (int j = 0; j < m - 1; j++) {
		s /= j + 1; u *= h;
		if (j < m - 2) v[j + m2] = s;
		v[j] = (1 - u) * s;
	}
	v[m - 1] = (1 - 2 * u * h + ((h > 0.5) ? pow(2 * h - 1, m) : 0)) * s / m;
	// End initialize
	// Now iterate
	for (int i = 1; i <= n; i++) {
		s = ((double)i) / n; u = v[m];
		// q[0] = F77_NAME(ddot)(&m, v, &one, q, &one) * s; // No shift 
		v[m] = NRdotprod(m, v, 0, v, m) * s;
		for (int j = 1; j < m - 1; j++) {
			double a = u;				// q[j - 1] 
			int m1 = m - 1 - j;
			u = v[j + m];
			// a += F77_NAME(ddot)(&m1, w, &one, q + j, &one) + v[m1] * q[m - 1];
			int mj = m + j;
			a += NRdotprod(m1, v, m2, v, mj) + v[m1] * v[m2 - 1];
			v[m + j] = a * s;
		}
		if (m > 1)								// shift? 
			v[m2 - 1] = (v[0] * v[m2 - 1] + u) * s;
		// check for under/overflows
		if (v[m + k - 1] > OVERSCALE) {
			double alpha = UNDERSCALE;
			// F77_NAME(dscal)(&m, &alpha, q, &one);
			NRmultconst(m, v, alpha, m);
			ns++;
		}
		if (v[m + k - 1] < UNDERSCALE) {
			double alpha = OVERSCALE;
			// F77_NAME(dscal)(&m, &alpha, q, &one);
			NRmultconst(m, v, alpha, m);
			ns--;
		}
	}
	// Return result
	s = v[m + k - 1];
	if (ns != 0) s *= pow(OVERSCALE, ns);		// Rescale if necessary
	return s;

}
double  Kolc(double& d, long& n1)
{// For external calls from Windows
	int n = static_cast<int>(n1);
	return KolC(d, n);
}
// Kolmogorov-Smirnov Statistics 
double KolStat(const int &N1, const double Value[], const double Cum[])
{// Value must be ordered and cdf Cum is required
	double KS = abs(Cum[0]), XN = static_cast<double>(N1);
	for (int i = 1; i < N1 - 1;i++)
	{
		double Temp = static_cast<double>(i + 1) / XN, Temp1 = abs(Temp - Cum[i]), Temp2 = abs(Temp - Cum[i+1]);
		if (Temp1 < Temp2) Temp1 = Temp2;
		if (Temp1 > KS) KS = Temp1;
	}
	double Temp = abs(1.0 - Cum[N1-1]);
	if (Temp > KS) KS = Temp;
	return KS;
}
double  KolSmir(const long &N1, const double Value[], const double Cum[])
{// For external calls from Windows
	int n = static_cast<int>(N1);
	return KolStat(n, Value, Cum);
}

/* ----------------------------------------------------------------
   Inversa de la probabilidad acumulada de la distribución normal
   ----------------------------------------------------------------*/
double  INormAc(const double &mu, const double &sigma, const double &p)
{
	double p1 = 2*p - 1.0; Erf zeta;
	return mu + sigma * RA2 * zeta.inverf(p1);
}
void sortdll(double it[], const long& L, const long& r, long mark[])
// Ordena el arreglo it entre l y r, en mark están los índices de los elementos
// Se usa para llamadas a librería dll
{
	double x, y;
	long i, j, k;
	long temp;
	i = L;
	j = r;
	temp = (L + r) / 2;
	x = it[temp];
	do
	{
		while (it[i] < x)
		{
			i = i + 1;
		}
		while (x < it[j])
		{
			j = j - 1;
		}
		if (i <= j)
		{
			y = it[i];
			it[i] = it[j];
			it[j] = y;
			k = mark[i];
			mark[i] = mark[j];
			mark[j] = k;
			i = i + 1;
			j = j - 1;
		}
	} while (i <= j);
	if (L < j)
	{
		sortdll(it, L, j, mark);
	}
	if (i < r)
	{
		sortdll(it, i, r, mark);
	}
}
/* ----------------------------------------------------------------
   Logaritmo de la Verosimilitud de Johnson SU
   ----------------------------------------------------------------*/
double LogVerJoU(const int &N1, const double Value[], const double &MuN, const double &StDevN, const double &Shift, const double &Mult)
{// N1 is sample size, Value is sample
	if ((Mult == 0.0) || (StDevN <= 0.0)) return -BIGNUMBER;
	double Mu = MuN, b = Mult;
	if (b < 0.0) { b = -b, Mu = -Mu; }
	double XN = static_cast<double>(N1), Score = -XN *((log(2.0 * PI) / 2.0) + log(b) + log (StDevN));
	double T0 = 0.0, T1 = 0.0, T2 = 0.0, T3 = 0.0;
	for (int i = 0; i < N1; i++)
	{
		double a = Value[i];
		T0 = (Value[i] - Shift) / b;		
		T1 = sqrt(1.0 + T0 * T0);		Score -= log(T1);
		T2 = asinh(T0) - Mu;			T3 += T2 * T2;
	}
	Score = Score - T3 / (2 * StDevN * StDevN);
	return Score;
}