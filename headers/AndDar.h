/* This code is to compute the cdf of the Anderson-Darling statistic 
 It is constructed from the code of 
		Evaluating the Anderson-Darling Distribution
		George Marsaglia and John C. W. Marsaglia
*/
#include <stdio.h>
#include <math.h>

double  ADinf(const double &z);
/*
A procedure for evaluating the limiting distribution of the Anderson-Darling statistic
 A_n= -n-(1/n)[ln(x_1(1-x_n)+3ln(x_2(1-x_{n-1})+5ln(x_3(1-x_{n-2})+...+(2n-1)ln(x_n(1-x_1))]
where x_1<x_2<...<x_n is an ordered set of purported uniform [0,1) variates.
The function is ADinf(z)=lim_{n->infty} Pr[A_n<z]. About 15 digit accuracy.
If you don't need that much accuracy, use the quick-and-easy adinf(z).
ADinf uses a two-term recursion for coefficients in series for which initial values
require the complementary normal integral, included as cPhi(z).
Otherwise, use erfc() if your C compiler has one with adequate accuracy.
*/
double adinf(const double &z);                  // Short version of ADinf but less accuracy (4 decimals)
double errfix(const int &n, const double &x);
double AD(const int &n,const double &z);        // To perform Anderson-Darling test

double cPhi(const double &z);                   // prototype; listing follows main 

double ADf(const double &z,const int &j)
{                     // called by ADinf(); see article. 
double t,f,fnew,a,b,c,r, temp;
t=(4*j+1)*(4*j+1)*1.23370055013617/z;
if(t>150.) return 0.;
a=2.22144146907918*exp(-t)/sqrt(t);
temp = sqrt(2 * t);
b=3.93740248643060*2.*cPhi(temp);// initialization requires cPhi 
// if you have erfc(), replace 2*cPhi(sqrt(2*t)) with erfc(sqrt(t))
r=z*.125; f=a+b*r;
for(int i=1; i < 200; i++)
{
    c=((i-.5-t)*b+t*a)/i;
    a=b; b=c; r*=z/(8*i+8);
    if(fabs(r)<1e-40 || fabs(c)<1.e-40) return f;
    fnew=f+c*r;
    if(f==fnew) return f;
    f=fnew;            
}
return f;				
}

double ADinf(const double &z)
{
int j;
double ad,adnew,r;
if(z<.01) return 0.; /* avoids exponent limits; ADinf(.01)=.528e-52 */
r=1./z;
ad=r*ADf(z,0);
for(j=1;j<100;j++)
{
    r*=(.5-j)/j;
    adnew=ad+(4*j+1)*r*ADf(z,j);
    if(ad==adnew) return ad;
    ad=adnew;         
}
return ad;
}
/*
int main()
{
    double z;
    while(1) // try a few values for z: 
    {
    printf("Enter value for z:");
    scanf("%lf",&z);
    printf("     lim n->infty Prob(A_n < %f) = %18.16f\n",z,ADinf(z));
    }
}
*/

/*
Complementary normal distribution function
   cPhi(x) = integral from x to infinity of phi(x)=exp(-.5*t^2)/sqrt(2*pi)
 13-15 digit accuracy for abs(x)<16.
Stores R(0),R(2),R(4),...,R(16), with cPhi(x)=R(x)*phi(x), phi normal density,
then uses Taylor series for R(z+h)=R(z)+hR'(z)+(1/2)h^2R''(z)+...
with -1<h<1,and R(z) one of R(0),R(2),R(4),...,R(16) stored as v[0],v[1],...v[8].
Examples: cPhi(2.75) needs R(2.75) and 2.75=2+.75 so use h=.75 and R(2)=v[1],
          cPhi(3.3)  needs R(3.3) and 3.3=4-.7, so use h=-.7 and R(4)=v[2].
*/
   double cPhi(const double &x)
   {
     long double v[]={1.25331413731550025,.421369229288054473,.236652382913560671,
     .162377660896867462,.123131963257932296,.0990285964717319214,
     .0827662865013691773,.0710695805388521071,.0622586659950261958};
     double h,a,b,z,t,s,pwr;
     int i,j;
      j=(fabs(x)+1.)/2.;
      a=v[j];    z=2*j;  h=fabs(x)-z;
      b=z*a-1;   pwr=1;  s=a+h*b;
      for(i=2;i<100;i+=2)
      {/* begin i loop */
           a=(a+z*b)/i;
           b=(b+z*a)/(i+1);
           pwr=pwr*h*h;
           t=s;
           s=s+pwr*(a+h*b);
           if(s==t){ s*=exp(-.5*x*x-.91893853320467274178);
                     return (x>0) ? s: 1-s;
                   }
      }/* end i loop */
      return (x > 0) ? s : 1 - s;
   }

/* Short, practical version of full ADinf(z), z>0.   */
double adinf(const double &z)
{ if(z<2.) return exp(-1.2337141/z)/sqrt(z)*(2.00012+(.247105-  \
          (.0649821-(.0347962-(.011672-.00168691*z)*z)*z)*z)*z);
      /* max |error| < .000002 for z<2, (p=.90816...) */
 return
 exp(-exp(1.0776-(2.30695-(.43424-(.082433-(.008056 -.0003146*z)*z)*z)*z)*z));
      /* max |error|<.0000008 for 4<z<infinity */
}

/*
The procedure  errfix(n,x)  corrects the error caused
by using the asymptotic approximation, x=adinf(z).
Thus x+errfix(n,x) is uniform in [0,1) for practical purposes;
accuracy may be off at the 5th, rarely at the 4th, digit.
*/
double errfix(const int &n, const double &x)
{
double c,t;
	if(x>.8) return
(-130.2137+(745.2337-(1705.091-(1950.646-(1116.360-255.7844*x)*x)*x)*x)*x)/n;
    c=.01265+.1757/n;
       if(x<c){ t=x/c;
         t=sqrt(t)*(1.-t)*(49*t-102);
       return t*(.0037/(n*n)+.00078/n+.00006)/n;
	       }
t=(x-c)/(.8-c);
t=-.00022633+(6.54034-(14.6538-(14.458-(8.259-1.91864*t)*t)*t)*t)*t;
return t*(.04213+.01365/n)/n;
}

/* The function AD(n,z) returns Prob(A_n<z) where
    A_n = -n-(1/n)*[ln(x_1*z_1)+3*ln(x_2*z_2+...+(2*n-1)*ln(x_n*z_n)]
          z_1=1-x_n, z_2=1-x_(n-1)...z_n=1-x_1, and
    x_1<x_2<...<x_n is an ordered set of iid uniform [0,1) variates.
*/

double AD(const int &n,const double &z)
{
 double c,v,x, xn;
 xn = static_cast<double>(n);
 x = adinf(z);
   /* now x=adinf(z). Next, get v=errfix(n,x) and return x+v; */
 if(x>.8)
    {
     v=(-130.2137+(745.2337-(1705.091-(1950.646-(1116.360-255.7844*x)*x)*x)*x)*x)/xn;
     return x+v;
    }
 c=.01265+.1757/xn;
 if(x<c)
    { 
     v=x/c;
     v=sqrt(v)*(1.-v)*(49*v-102);
     return x+v*(.0037/(n*n)+.00078/n+.00006)/xn;
    }
 v=(x-c)/(.8-c);
 v=-.00022633+(6.54034-(14.6538-(14.458-(8.259-1.91864*v)*v)*v)*v)*v;
 return x+v*(.04213+.01365/n)/xn;
}

/* You must give the ADtest(int n, double *x) routine a sorted array
       Fx[0]<=Fx[1]<=..<=Fx[n-1]
    that you are testing for uniformity.
   It will return the p-value associated
   with the Anderson-Darling test, using
    the above adinf() and errfix( ,   )
         Not well-suited for n<7,
     (accuracy could drop to 3 digits).
*/
double ADtest(const int &n1, const double Fx[], double &z)
{// Computes Anderson-Darling statistics and returns p-value
    double n2 = static_cast<double>(n1);
    z = 0.0;
    for(int i=0; i < n1; i++)   
    {
      double i1 = static_cast<double>(i), t;
      t = Fx[i]*(1.- Fx[n1-1-i]);
      z -= (i1 + i1 +1.0)*log(t);
    }
  z = -n2 + z / n2;
  return AD(n1,z);
}
double AndDarAc(const long &n, const double &z)
{// Returns CDF of Anserson-Darling, z is the Anderson-Darling statistic
    int n1 = static_cast<int>(n);
    return AD(n1, z);
}
double CMtest(const int& n1, const double Fx[])
{// Computes Cramer-von Mises statistic, Fx-values must be ordered
    double n2 = static_cast<double>(n1), z = 0.0;
    for (int i = 0; i < n1; i++)
    {
        double i1 = static_cast<double>(i+1), t = i1 - 0.5 - n2 * Fx[i];
        z += t * t;
    }
    z /= (n2 * n2);
    z += (1.0 / (12.0 * n2));
    return z;
}
double CM(const long& n, const double Fx[])
{
    int n1 = static_cast<int>(n);
    return CMtest(n1, Fx);
}
int SkewKurt(const int& n1, const double DataVal[], double &Mean, double &StDev, double &Skew, double &Kurt)
{// Computes Mean, StDev, Skewness and Kurtosis of DataVal
    Mean = 0.0; StDev = 0.0; Skew = 0.0; Kurt = 0.0;
    if (n1 < 2) return -1;
    double xn = static_cast<double>(n1);
    for (int i = 0; i < n1; i++)
    {
        Mean += DataVal[i];     StDev += DataVal[i] * DataVal[i];
    }
    StDev = sqrt((StDev - Mean * Mean / xn) / (xn - 1.0));
    Mean /= xn;
    if (n1 < 3) return -2;
    for (int i = 0; i < n1; i++)
    {
        double Temp = DataVal[i] - Mean, Temp2 = Temp * Temp;
        Skew += Temp * Temp2;                                    // To compute sum of third deviations
        Kurt += Temp2 * Temp2;                                   // To compute sum of fourth deviations
    }
    Skew = sqrt((xn - 1.0) / xn) * Skew / ((xn - 2.0) * StDev * StDev * StDev);
    if (n1 < 4) return -3;
    double Temp = (xn - 1.0) * StDev * StDev;
    Kurt = xn * (xn * xn - 1) * Kurt / (Temp * Temp * (xn - 2.0) * (xn - 3.0));
    return 0;
}


