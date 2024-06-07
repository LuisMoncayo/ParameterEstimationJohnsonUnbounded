double polev1(double &x, double coef[], int &N1)
{// Used in function psi
	double ans; long i; double *p;
	p = coef; ans = *p++; i = N1;
	do
		ans = ans * x + *p++;
	while (--i);
	return(ans);
}
#define MAXNUM   9999999999999999					// Used in function psi
#define PI  3.14159265358979323846					// Se usa en varios procedimientos
#define EUL 0.57721566490153286061					// Constante de Euler
double psi(double x)     // NOTE: do not call for x negative
{// Digamma function
	static double A[] = {
	 8.33333333333333333333E-2,
	-2.10927960927960927961E-2,
	 7.57575757575757575758E-3,
	-4.16666666666666666667E-3,
	 3.96825396825396825397E-3,
	-8.33333333333333333333E-3,
	 8.33333333333333333333E-2
	};
	/* double floor(), log(), tan(), polevl();
	extern double PI, MAXNUM; */

	double p, q, nz, s, w, y, z;
	int n, negative, seis;
	negative = 0; nz = 0.0; seis = 6;
	if (x <= 0.0)
	{
		negative = 1; q = x; p = floor(q);
		if (p == q)
		{ //mtherr( "psi", SING );
			return((double)MAXNUM);
		}
		/* Remove the zeros of tan(PI x)
		 * by subtracting the nearest integer from x
		*/
		nz = q - p;
		if (nz != 0.5)
		{
			if (nz > 0.5)
			{
				p += 1.0;
				nz = q - p;
			}
			nz = PI / tan(PI*nz);
		}
		else
		{
			nz = 0.0;
		}
		x = 1.0 - x;
	}

	/* check for positive integer up to 10 */

	if ((x <= 10.0) && (x == floor(x)))
	{
		y = 0.0; n = static_cast<int>(x);
		for (int i = 1; i < n; i++)
		{
			w = i; y += 1.0 / w;
		}
		y -= EUL;
		goto done;
	}
	s = x; w = 0.0;
	while (s < 10.0)
	{
		w += 1.0 / s;
		s += 1.0;
	}

	if (s < 1.0e17)
	{
		z = 1.0 / (s * s);
		y = z * polev1(z, A, seis);
	}
	else
		y = 0.0;

	y = log(s) - (0.5 / s) - y - w;

done:

	if (negative)
	{
		y -= nz;
	}

	return(y);
}
double gammln(const double &xx) {
	int j;
	double x, tmp, y, ser;
	static const double cof[14] = { 57.1562356658629235,-59.5979603554754912,
		14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
		.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
		-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
		.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5 };
	if (xx <= 0) throw("bad arg in gammln");
	y = x = xx;
	tmp = x + 5.24218750000000000;
	tmp = (x + 0.5)*log(tmp) - tmp;
	ser = 0.999999999999997092;
	for (j = 0;j<14;j++) ser += cof[j] / ++y;
	return tmp + log(2.5066282746310005*ser / x);
}
struct Gauleg18 {
	static const int ngau = 18;
	static const double y[18];
	static const double w[18];
};
const double Gauleg18::y[18] = {0.0021695375159141994,
0.011413521097787704,0.027972308950302116,0.051727015600492421,
0.082502225484340941, 0.12007019910960293,0.16415283300752470,
0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
0.87126389619061517, 0.95698180152629142};
const double Gauleg18::w[18] = {0.0055657196642445571,
0.012915947284065419,0.020181515297735382,0.027298621498568734,
0.034213810770299537,0.040875750923643261,0.047235083490265582,
0.053244713977759692,0.058860144245324798,0.064039797355015485,
0.068745323835736408,0.072941885005653087,0.076598410645870640,
0.079687828912071670,0.082187266704339706,0.084078218979661945,
0.085346685739338721,0.085983275670394821};
struct Gamma : Gauleg18 {
	static const int ASWITCH=100;
	static const double EPS;
	static const double FPMIN;
	double gln;

	double gammp(const double &a, const double &x) {
		if (x < 0.0 || a <= 0.0) throw("bad args in gammp");
		if (x == 0.0) return 0.0;
		else if ((int)a >= ASWITCH) return gammpapprox(a,x,1);
		else if (x < a+1.0) return gser(a,x);
		else return 1.0-gcf(a,x);
	}

	double gammq(const double &a, const double &x) {
		if (x < 0.0 || a <= 0.0) throw("bad args in gammq");
		if (x == 0.0) return 1.0;
		else if ((int)a >= ASWITCH) return gammpapprox(a,x,0);
		else if (x < a+1.0) return 1.0-gser(a,x);
		else return gcf(a,x);
	}

	double gser(const double &a, const double &x) {
		double sum,del,ap;
		gln=gammln(a);
		ap=a;
		del=sum=1.0/a;
		for (;;) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				return sum*exp(-x+a*log(x)-gln);
			}
		}
	}

	double gcf(const double &a, const double &x) {
		int i;
		double an,b,c,d,del,h;
		gln=gammln(a);
		b=x+1.0-a;
		c=1.0/FPMIN;
		d=1.0/b;
		h=d;
		for (i=1;;i++) {
			an = -i*(i-a);
			b += 2.0;
			d=an*d+b;
			if (fabs(d) < FPMIN) d=FPMIN;
			c=b+an/c;
			if (fabs(c) < FPMIN) c=FPMIN;
			d=1.0/d;
			del=d*c;
			h *= del;
			if (fabs(del-1.0) <= EPS) break;
		}
		return exp(-x+a*log(x)-gln)*h;
	}

	double gammpapprox(const double &a, const double &x, const int &psig) {
		int j;
		double xu,t,sum,ans;
		double a1 = a-1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
		gln = gammln(a);
		if (x > a1) xu = MAX(a1 + 11.5*sqrta1, x + 6.0*sqrta1);
		else xu = MAX(0.,MIN(a1 - 7.5*sqrta1, x - 5.0*sqrta1));
		sum = 0;
		for (j=0;j<ngau;j++) {
			t = x + (xu-x)*y[j];
			sum += w[j]*exp(-(t-a1)+a1*(log(t)-lna1));
		}
		ans = sum*(xu-x)*exp(a1*(lna1-1.)-gln);
		return (psig?(ans>0.0? 1.0-ans:-ans):(ans>=0.0? ans:1.0+ans));
	}

	double invgammp(const double &p, const double &a);

};
const double Gamma::EPS = std::numeric_limits<double>::epsilon();
const double Gamma::FPMIN = DBL_MIN/EPS; //const double Gamma::FPMIN = numeric_limits<double>::min()/EPS;
double Gamma::invgammp(const double &p, const double &a) {
	int j;
	double x,err,t,u,pp,lna1,afac,a1=a-1;
	const double EPS=1.e-8;
	gln=gammln(a);
	if (a <= 0.) throw("a must be pos in invgammap");
	if (p >= 1.) return MAX(100.,a + 100.*sqrt(a));
	if (p <= 0.) return 0.0;
	if (a > 1.) {
		lna1=log(a1);
		afac = exp(a1*(lna1-1.)-gln);
		pp = (p < 0.5)? p : 1. - p;
		t = sqrt(-2.*log(pp));
		x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
		if (p < 0.5) x = -x;
		x = MAX(1.e-3,a*pow(1.-1./(9.*a)-x/(3.*sqrt(a)),3));
	} else {
		t = 1.0 - a*(0.253+a*0.12);
		if (p < t) x = pow(p/t,1./a);
		else x = 1.-log(1.-(p-t)/(1.-t));
	}
	for (j=0;j<12;j++) {
		if (x <= 0.0) return 0.0;
		err = gammp(a,x) - p;
		if (a > 1.) t = afac*exp(-(x-a1)+a1*(log(x)-lna1));
		else t = exp(-x+a1*log(x)-gln);
		u = err/t;
		x -= (t = u/(1.-0.5*MIN(1.,u*((a-1.)/x - 1))));
		if (x <= 0.) x = 0.5*(x + t);
		if (fabs(t) < EPS*x ) break;
	}
	return x;
}
struct Beta : Gauleg18 {
	static const int SWITCH=3000;
	static const double EPS, FPMIN;

	double betai(const double &a, const double &b, const double &x) {
		double bt;
		if (a <= 0.0 || b <= 0.0) throw("Bad a or b in routine betai");
		if (x < 0.0 || x > 1.0) throw("Bad x in routine betai");
		if (x == 0.0 || x == 1.0) return x;
		if (a > SWITCH && b > SWITCH) return betaiapprox(a,b,x);
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
		if (x < (a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
		else return 1.0-bt*betacf(b,a,1.0-x)/b;
	}

	double betacf(const double &a, const double &b, const double &x) {
		int m,m2;
		double aa,c,d,del,h,qab,qam,qap;
		qab=a+b;
		qap=a+1.0;
		qam=a-1.0;
		c=1.0;
		d=1.0-qab*x/qap;
		if (fabs(d) < FPMIN) d=FPMIN;
		d=1.0/d;
		h=d;
		for (m=1;m<10000;m++) {
			m2=2*m;
			aa=m*(b-m)*x/((qam+m2)*(a+m2));
			d=1.0+aa*d;
			if (fabs(d) < FPMIN) d=FPMIN;
			c=1.0+aa/c;
			if (fabs(c) < FPMIN) c=FPMIN;
			d=1.0/d;
			h *= d*c;
			aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
			d=1.0+aa*d;
			if (fabs(d) < FPMIN) d=FPMIN;
			c=1.0+aa/c;
			if (fabs(c) < FPMIN) c=FPMIN;
			d=1.0/d;
			del=d*c;
			h *= del;
			if (fabs(del-1.0) <= EPS) break;
		}
		return h;
	}

	double betaiapprox(const double &a, const double &b, const double &x) {
		int j;
		double xu,t,sum,ans;
		double a1 = a-1.0, b1 = b-1.0, mu = a/(a+b);
		double lnmu=log(mu),lnmuc=log(1.-mu);
		t = sqrt(a*b/(SQR(a+b)*(a+b+1.0)));
		if (x > a/(a+b)) {
			if (x >= 1.0) return 1.0;
			xu = MIN(1.,MAX(mu + 10.*t, x + 5.0*t));
		} else {
			if (x <= 0.0) return 0.0;
			xu = MAX(0.,MIN(mu - 10.*t, x - 5.0*t));
		}
		sum = 0;
		for (j=0;j<18;j++) {
			t = x + (xu-x)*y[j];
			sum += w[j]*exp(a1*(log(t)-lnmu)+b1*(log(1-t)-lnmuc));
		}
		ans = sum*(xu-x)*exp(a1*lnmu-gammln(a)+b1*lnmuc-gammln(b)+gammln(a+b));
		return ans>0.0? 1.0-ans : -ans;
	}

	double invbetai(const double &p, const double &a, const double &b) {
		const double EPS = 1.e-8;
		double pp,t,u,err,x,al,h,w,afac,a1=a-1.,b1=b-1.;
		int j;
		if (p <= 0.) return 0.;
		else if (p >= 1.) return 1.;
		else if (a >= 1. && b >= 1.) {
			pp = (p < 0.5)? p : 1. - p;
			t = sqrt(-2.*log(pp));
			x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
			if (p < 0.5) x = -x;
			al = (SQR(x)-3.)/6.;
			h = 2./(1./(2.*a-1.)+1./(2.*b-1.));
			w = (x*sqrt(al+h)/h)-(1./(2.*b-1)-1./(2.*a-1.))*(al+5./6.-2./(3.*h));
			x = a/(a+b*exp(2.*w));
		} else {
			double lna = log(a/(a+b)), lnb = log(b/(a+b));
			t = exp(a*lna)/a;
			u = exp(b*lnb)/b;
			w = t + u;
			if (p < t/w) x = pow(a*w*p,1./a);
			else x = 1. - pow(b*w*(1.-p),1./b);
		}
		afac = -gammln(a)-gammln(b)+gammln(a+b);
		for (j=0;j<10;j++) {
			if (x == 0. || x == 1.) return x;
			err = betai(a,b,x) - p;
			t = exp(a1*log(x)+b1*log(1.-x) + afac);
			u = err/t;
			x -= (t = u/(1.-0.5*MIN(1.,u*(a1/x - b1/(1.-x)))));
			if (x <= 0.) x = 0.5*(x + t);
			if (x >= 1.) x = 0.5*(x + t + 1.);
			if (fabs(t) < EPS*x && j > 0) break;
		}
		return x;
	}

};
const double Beta::EPS = std::numeric_limits<double>::epsilon();
const double Beta::FPMIN = DBL_MIN/EPS; //const double Beta::FPMIN = numeric_limits<double>::min()/EPS;
struct Gammadist : Gamma {
	double alph, bet, fac;
	Gammadist(double aalph, double bbet = 1.) : alph(aalph), bet(bbet) {
		if (alph <= 0. || bet <= 0.) throw("bad alph,bet in Gammadist");
		fac = alph*log(bet)-gammln(alph);
	}
	double p(double x) {
		if (x <= 0.) throw("bad x in Gammadist");
		return exp(-bet*x+(alph-1.)*log(x)+fac);
	}
	double cdf(double x) {
		if (x < 0.) throw("bad x in Gammadist");
		return gammp(alph,bet*x);
	}
	double invcdf(double p) {
		if (p < 0. || p >= 1.) throw("bad p in Gammadist");
		return invgammp(p,alph)/bet;
	}
};
struct Betadist : Beta {
	double alph, bet, fac;
	Betadist(double aalph, double bbet) : alph(aalph), bet(bbet) {
		if (alph <= 0. || bet <= 0.) throw("bad alph,bet in Betadist");
		fac = gammln(alph+bet)-gammln(alph)-gammln(bet);
	}
	double p(double x) {
		if (x <= 0. || x >= 1.) throw("bad x in Betadist");
		return exp((alph-1.)*log(x)+(bet-1.)*log(1.-x)+fac);
	}
	double cdf(double x) {
		if (x < 0. || x > 1.) throw("bad x in Betadist");
		return betai(alph,bet,x);
	}
	double invcdf(double p) {
		if (p < 0. || p > 1.) throw("bad p in Betadist");
		return invbetai(p,alph,bet);
	}
};
struct Studenttdist : Beta {
	double nu, mu, sig, np, fac;
	Studenttdist(const double &nnu, double mmu = 0., double ssig = 1.)
	: nu(nnu), mu(mmu), sig(ssig) {
		if (sig <= 0. || nu <= 0.) throw("bad sig,nu in Studentdist");
		np = 0.5*(nu + 1.);
		fac = gammln(np)-gammln(0.5*nu);
	}
	double p(const double &t) {
		return exp(-np*log(1.+SQR((t-mu)/sig)/nu)+fac)
			/(sqrt(3.14159265358979324*nu)*sig);
	}
	double cdf(const double &t) {
		double p = 0.5*betai(0.5*nu, 0.5, nu/(nu+SQR((t-mu)/sig)));
		if (t >= mu) return 1. - p;
		else return p;
	}
	double invcdf(const double &p) {
		if (p <= 0. || p >= 1.) throw("bad p in Studentdist");
		double x = invbetai(2.*MIN(p,1.-p), 0.5*nu, 0.5);
		x = sig*sqrt(nu*(1.-x)/x);
		return (p >= 0.5? mu+x : mu-x);
	}
	double aa(const double &t) {
		if (t < 0.) throw("bad t in Studentdist");
		return 1.-betai(0.5*nu, 0.5, nu/(nu+SQR(t)));
	}
	double invaa(const double &p) {
		if (p < 0. || p >= 1.) throw("bad p in Studentdist");
		double x = invbetai(1.-p, 0.5*nu, 0.5);
		return sqrt(nu*(1.-x)/x);
	}
};
struct Poissondist : Gamma {
	double lam;
	Poissondist(const double &llam) : lam(llam) {
		if (lam <= 0.) throw("bad lam in Poissondist");	
	}
	double p(int n) {
		if (n < 0) throw("bad n in Poissondist");
		return exp(-lam + n*log(lam) - gammln(n+1.));
	}
	double cdf(int n) {
		if (n < 0) throw("bad n in Poissondist");
		if (n == 0) return 0.;
		return gammq((double)n,lam);
	}
	int invcdf(double p) {
		int n,nl,nu,inc=1;
		if (p <= 0. || p >= 1.) throw("bad p in Poissondist");
		if (p < exp(-lam)) return 0;
		n = (int)MAX(sqrt(lam),5.);
		if (p < cdf(n)) {
			do {
				n = MAX(n-inc,0);
				inc *= 2;
			} while (p < cdf(n));
			nl = n; nu = n + inc/2;
		} else {
			do {
				n += inc;
				inc *= 2;
			} while (p > cdf(n));
			nu = n; nl = n - inc/2;
		}
		while (nu-nl>1) {
			n = (nl+nu)/2;
			if (p < cdf(n)) nu = n;
			else nl = n;
		}
		return nl;
	}
};
struct Binomialdist : Beta {
	int n;
	double pe, fac;
	Binomialdist(int nn, double ppe) : n(nn), pe(ppe) {
		if (n <= 0 || pe <= 0. || pe >= 1.) throw("bad args in Binomialdist");
		fac = gammln(n+1.);
	}
	double p(int k) {
		if (k < 0) throw("bad k in Binomialdist");
		if (k > n) return 0.;
		return exp(k*log(pe)+(n-k)*log(1.-pe)
			+fac-gammln(k+1.)-gammln(n-k+1.));
	}
	double cdf(int k) {
		if (k < 0) throw("bad k in Binomialdist");
		if (k == 0) return 0.;
		if (k > n) return 1.;
		return 1. - betai((double)k,n-k+1.,pe);
	}
	int invcdf(double p) {
		int k,kl,ku,inc=1;
		if (p <= 0. || p >= 1.) throw("bad p in Binomialdist");
		k = MAX(0,MIN(n,(int)(n*pe)));
		if (p < cdf(k)) {
			do {
				k = MAX(k-inc,0);
				inc *= 2;
			} while (p < cdf(k));
			kl = k; ku = k + inc/2;
		} else {
			do {
				k = MIN(k+inc,n+1);
				inc *= 2;
			} while (p > cdf(k));
			ku = k; kl = k - inc/2;
		}
		while (ku-kl>1) {
			k = (kl+ku)/2;
			if (p < cdf(k)) ku = k;
			else kl = k;
		}
		return kl;
	}
};
struct Chisqdist : Gamma {
	double nu,fac;
	Chisqdist(double nnu) : nu(nnu) {
		if (nu <= 0.) throw("bad nu in Chisqdist");
		fac = 0.693147180559945309*(0.5*nu)+gammln(0.5*nu);
	}
	double p(double x2) {
		if (x2 <= 0.) throw("bad x2 in Chisqdist");
		return exp(-0.5*(x2-(nu-2.)*log(x2))-fac);
	}
	double cdf(double x2) {
		if (x2 < 0.) throw("bad x2 in Chisqdist");
		return gammp(0.5*nu,0.5*x2);
	}
	double invcdf(double p) {
		if (p < 0. || p >= 1.) throw("bad p in Chisqdist");
		return 2.*invgammp(p,0.5*nu);
	}
};
struct Fdist : Beta {
	double nu1,nu2;
	double fac;
	Fdist(double nnu1, double nnu2) : nu1(nnu1), nu2(nnu2) {
		if (nu1 <= 0. || nu2 <= 0.) throw("bad nu1,nu2 in Fdist");
		fac = 0.5*(nu1*log(nu1)+nu2*log(nu2))+gammln(0.5*(nu1+nu2))
			-gammln(0.5*nu1)-gammln(0.5*nu2);
	}
	double p(double f) {
		if (f <= 0.) throw("bad f in Fdist");
		return exp((0.5*nu1-1.)*log(f)-0.5*(nu1+nu2)*log(nu2+nu1*f)+fac);
	}
	double cdf(double f) {
		if (f < 0.) throw("bad f in Fdist");
		return betai(0.5*nu1,0.5*nu2,nu1*f/(nu2+nu1*f));
	}
	double invcdf(double p) {
		if (p <= 0. || p >= 1.) throw("bad p in Fdist");
		double x = invbetai(p,0.5*nu1,0.5*nu2);
		return nu2*x/(nu1*(1.-x));
	}
};
