#include <Eigen/Core>
#include <iostream>
#include <nr3.h>
#include <LBFGSB1.h>
#include <fungammabeta.h>
#include <Functions.cpp>
#include <AndDar.h>
#include <roots.h>
using namespace LBFGSpp;
typedef double Scalar;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
const Scalar SMALLML = 1e-7;
class FunJoUCvM
{
private:
    int NDataJoU;
    Vector XJoU, GYJoU, W12JoU, ZJoU, F1JoU, UJoU, ValJoU; // yi, g(yi), 1/g'(yi), zi, F'(yi), 1 + yi * yi
    Scalar XN = 0.0;
    bool FaultJoU = false;
public:
    Scalar a, b, c, d;
    FunJoUCvM(const int N1, Vector Value) : NDataJoU(N1), ValJoU(Value) {
        XJoU.resize(N1); GYJoU.resize(N1); W12JoU.resize(N1);
        UJoU.resize(N1); ZJoU.resize(N1); F1JoU.resize(N1);
        XN = static_cast<double>(NDataJoU);
    }
    double operator()(const Vector& x, Vector& grad)
    {
        a = x[0]; b = x[1]; c = x[2];  d = x[3];
        if ((x[1] == 0.0) || (NDataJoU < 2) || (x[3] == 0.0)) { FaultJoU = true; return 0.0; }
       // if (b <= SMALLML) b = SMALLML;      if (d <= SMALLML) d = SMALLML;
        Scalar b2 = SQR(b);
        double Score = 0.0;     XN = static_cast<Scalar>(NDataJoU);
        for (int i = 0; i < NDataJoU; i++)
        {
            XJoU[i] = (ValJoU[i] - a);
            W12JoU[i] = sqrt(b2 + XJoU[i] * XJoU[i]);
            UJoU[i] = log(XJoU[i] + W12JoU[i]) - log(b) - c;
            ZJoU[i] = UJoU[i] / d;
            double xi = static_cast<double>(i);
            GYJoU[i] = (xi + 0.5) / XN - CumNHart(ZJoU[i], 0.0, 1.0);			// (i - 1/2) / N - F(Xi)
            Score += GYJoU[i] * GYJoU[i];
            F1JoU[i] = DNor(ZJoU[i], 0.0, 1.0);									// F'(Zi) 
        }
        grad[0] = 0.0;		grad[1] = 0.0;		grad[2] = 0.0;		grad[3] = 0.0;
        for (int i = 0; i < NDataJoU; i++)
        {
            Scalar Temp = GYJoU[i] * F1JoU[i], Temp1 = Temp / W12JoU[i];
            grad[0] += Temp1;
            grad[1] += XJoU[i] * Temp1;
            grad[2] += Temp;
            grad[3] += Temp * UJoU[i];
        }
        Scalar Divisor = 0.5 * d, Direction = Divisor * b;
        grad[0] /= Divisor;
        grad[1] /= Direction;
        grad[2] /= Direction;
        grad[3] /= (Direction * d);
        a = grad[0]; b = grad[1]; c = grad[2]; d = grad[3];
        return Score;
    }
    bool getFault() { return FaultJoU; }
};
inline double InSinh(const double& A, const double& B)
{// Returns function for MultJo
    double Value = A + sqrt(SQR(A) + SQR(B));
    if (Value <= numeric_limits<double>::lowest()) return numeric_limits<double>::lowest();
    else return Value;
}
/* ---------------------------------------------------------------------
   Calcula estimadores de una Johnson SU por máxima verosimilitud
   dados parámetros de traslado y escala
   ---------------------------------------------------------------------*/
bool FitJoUMuSt(const int& N1, Vector Value, const double& Shift, const double& Mult, double& MuN, double& StDevN)
{//Estima los parámetros de una distribución Johnson no acotada.
 //Entrada: N1, Value, Shift, Mult. Salida: MuN, StDevN, Shift debe ser menor que Mínimo Value
    if ((Mult <= 0.0) || (N1 < 2)) return false;
    double LogMult = log(Mult);
    MuN = 0.0; StDevN = 0.0;
    double XN1 = static_cast<double>(N1);
    for (int i = 0; i < N1; i++)
    {
        double Temp = (Value[i] - Shift);
        Temp = log(InSinh(Temp, Mult));
        MuN += Temp; StDevN += SQR(Temp);
        double aa = MuN, bb = StDevN;
    }
    double Temp = (StDevN - SQR(MuN) / XN1) / (XN1 - 1.0);
    if (Temp <= numeric_limits<double>::lowest()) StDevN = sqrt(numeric_limits<double>::lowest());
    else StDevN = sqrt(Temp);
    MuN = (MuN / XN1) - log(Mult);
    return true;
}
double FunJo1, FunJo2, FunJo3;
double MultJo(const double& X)
{// left hand to solve MultJou(X) = 0 using zbrent
    return  InSinh(FunJo1, X) * InSinh(FunJo3, X) - SQR(InSinh(FunJo2, X));
}
Scalar TryInvalC(const int& N1, const int& I1, Vector Value, const double& Quan1, const double& Quan2,
    const double& Quan3, double& MuY, double& SigmaY, double& Shift, double& Mult)
{// To try initial value for Shift = Value[i], and Cramér-von Mises
    Shift = Value[I1];
    FunJo1 = Quan1 - Shift;				FunJo2 = Quan2 - Shift;				FunJo3 = Quan3 - Shift;
    Doub Tol = 3.0e-8, MaxA = abs(FunJo3), B1 = abs(FunJo1);
    if (B1 > MaxA)	MaxA = B1;
    int NIter = 0;
    if (MultJo(0.0) == 0.0)		Mult = 1.0;
    else Mult = zbrent(MultJo, 0.0, MaxA, Tol, NIter);
    if (NIter < 0)	Mult = 1.0;
    if (!FitJoUMuSt(N1, Value, Shift, Mult, MuY, SigmaY))		return 9.0e16;
    FunJoUCvM funJoU(N1, Value);
    Vector p(4), grad(4);		p[0] = Shift;		p[1] = Mult;		p[2] = MuY;		p[3] = SigmaY;
    Scalar CvM = funJoU(p, grad);
    return CvM;
}
int InValC(const int& N1, Vector Value, double& MuY, double& SigmaY, double& ShiftY, double& MultY)
{// Initial Johnson parameters from quantiles trying to minimize Cramér-von Mises
// Computation of quantiles
    if (N1 < 5) return -1;
    double XN = static_cast<double>(N1);
    int N0 = static_cast<int>(XN / 2.0);
    double XN0 = static_cast<double>(N0), Quan1, Quan2, Quan3, Zn = 0.524, Pn = 0.699860730068;		// Pn = CuMNHart(Zn)
    if (XN0 == XN / 2) Quan2 = (Value[N0] + Value[N0 - 1]) / 2.0;
    else Quan2 = Value[N0];
    if (N1 > 100) {
        Zn = 0.85; Pn = 0.802337456877;
    }
    double XM = XN * Pn - 0.5;
    int M = static_cast<int>(XM);
    Quan3 = Value[M];		Quan1 = Value[N1 - M];
    // Searching for initial values 
    Scalar ScoreY = TryInvalC(N1, 0, Value, Quan1, Quan2, Quan3, MuY, SigmaY, ShiftY, MultY);
    for (int i = 1; i < N1; i++)
    {
        double Mu, Sigma, Shift, Mult;
        Scalar Score = TryInvalC(N1, i, Value, Quan1, Quan2, Quan3, Mu, Sigma, Shift, Mult);
        if (Score < ScoreY) {
            MuY = Mu;	SigmaY = Sigma;		ShiftY = Shift;		MultY = Mult;
            ScoreY = Score;
        }
    }
    return 0;
}
/* ---------------------------------------------------------------------
   Calcula estimadores de forma de una Johnson SU
   ---------------------------------------------------------------------*/
bool FitParWheeler(const int& N1, const Vector Value, Scalar& MuN, Scalar& StDevN)
{//Estima los parámetros de forma de una distribución Johnson no acotada.
 //Entrada: N1, Value. Salida: MuN, StDevN
    if (N1 < 5) return false;
    MuN = 0.0; StDevN = 1.0;
    double XN = static_cast<double>(N1), Pn = (XN - 0.5) / XN, Zn = INormAc(0.0, 1.0, Pn), Temp = Pn * XN + 0.5;
    int TempI = static_cast<int>(Temp);
    if (TempI > N1 - 1) TempI = N1 - 1;			double Quan4 = Value[TempI];
    Temp = (1 - Pn) * XN + 0.5;					TempI = static_cast<int>(Temp);			double Quan0 = Value[TempI];
    Pn = CumNHart(Zn / 2.0, 0.0, 1.0);
    Temp = Pn * XN + 0.5;						TempI = static_cast<int>(Temp);
    if (TempI > N1 - 1) TempI = N1 - 1;			double Quan3 = Value[TempI];
    Temp = (1 - Pn) * XN + 0.5;					TempI = static_cast<int>(Temp);			double Quan1 = Value[TempI];
    Temp = 0.5 * (XN + 1.0);					TempI = static_cast<int>(Temp);			double Quan2 = Value[TempI];
    double Tu = Quan4 - Quan0;
    if (Quan3 != Quan1)		Tu /= (Quan3 - Quan1);
    Temp = Tu / 2;
    if (Temp > 1.0) {										// StDevN = 1 otherwise
        Temp += sqrt(Temp * Temp - 1);	// Temp es b
        StDevN = 2 * log(Temp) / Zn;
    }
    Tu = Quan4 - Quan2;
    if (Quan2 != Quan0)		Tu /= (Quan2 - Quan0);			// Tu es t
    Temp = (1 - Tu * Temp * Temp) / (Tu - Temp * Temp);		// Temp es a^2
    if (Temp >= 0)	MuN = log(sqrt(Temp));					// MuN = 0 otherwise
    return true;
}
bool InValJoU(const int& N1, const Vector Value, Scalar& MuY, Scalar& SigmaY, Vector& p)
{// Initial values using proposal of Wheeler
    double ND = static_cast<double>(N1), P1 = 0.699860730068, P1neg = 1 - P1;
    //P1 = CumNHart(0.524, 0.0, 1.0);
    MuY = 0.0; SigmaY = 1.0;
    if (!FitParWheeler(N1, Value, MuY, SigmaY)) return false;
    int i1 = static_cast<int>(P1 * ND + 0.5), i1neg = static_cast<int>(P1neg * ND + 0.5);
    double x1 = Value[i1], x1neg = Value[i1neg], z1 = sinh(MuY + 0.524 * SigmaY), z1neg = sinh(MuY - 0.524 * SigmaY);
    if (SigmaY > 0) p[1] = (x1 - x1neg) / (z1 - z1neg);
    p[0] = x1 - p[1] * z1;
    return true;
}
class FunJoUML
{public:
    int NDataJoU;								bool Fault = false;
    double XN, StDevJo, MuJo;
    Vector XJoU, WJoU, W2JoU, GUJoU, ValJoU;	// xi, wi, wi^0.5, g
    FunJoUML(int N1, Vector DataVal) : NDataJoU(N1), ValJoU(DataVal) {
        XJoU.resize(N1);					WJoU.resize(N1);	W2JoU.resize(N1);	GUJoU.resize(N1);
        XN = static_cast<double>(NDataJoU);		MuJo = 0.0;				StDevJo = 0.0;
    }
    Scalar operator()(const Vector& x, Vector& grad)
    {
        Scalar b = x[1];
        if (b <= SMALLML)  b = SMALLML;
        if (NDataJoU < 2) {Fault = true; return 3.0e16; }
        Scalar XN = static_cast<Scalar>(NDataJoU), Score = 0.0, S1 = 0.0, S2 = 0.0, S3 = 0.0, 
            SP1 = 0.0, SP2 = 0.0, b2 = SQR(b);
        MuJo = 0.0, StDevJo = 0.0;
        for (int i = 0; i < NDataJoU; i++)
        {
            XJoU[i] = (ValJoU[i] - x[0]);					WJoU[i] = b2 + SQR(XJoU[i]);	 // SQR is square
            W2JoU[i] = sqrt(WJoU[i]);                   Scalar Temp = XJoU[i] + W2JoU[i],Temp1 = W2JoU[i] * Temp;
            GUJoU[i] = log(Temp);						// Error if Temp <=0
            MuJo += GUJoU[i];								StDevJo += SQR(GUJoU[i]);
            Score += log(WJoU[i]);
            S1 += GUJoU[i];								S2 += (1.0 / W2JoU[i]);					S3 += (1.0 / Temp1);
            SP1 += (GUJoU[i] / W2JoU[i]);				SP2 += GUJoU[i] / Temp1;
        }
        Score *= 0.5;
        StDevJo = (StDevJo - MuJo * MuJo / XN) / XN;
        SP1 -= (S1 * S2) / XN;							    SP2 -= (S1 * S3) / XN;
        SP1 /= StDevJo;										SP2 /= StDevJo;
        StDevJo = sqrt(StDevJo);							MuJo = (MuJo / XN) - log(b);
        Score += XN * log(StDevJo);
        grad[0] = 0.0;									    grad[1] = 0.0;
        for (int i = 0; i < NDataJoU; i++)
        {
            Scalar Temp = XJoU[i] / WJoU[i];
            grad[0] -= Temp;
            grad[1] += 1 / WJoU[i];
        }
        grad[0] -= SP1;		grad[1] -= SP2;		grad[1] *= b;
        double a1 = grad[0], a2 = grad[1];
        return Score;
    }
    Scalar getMuJo() { return MuJo; }
    Scalar getStDevJo() { return StDevJo; }
    bool getFault() { return Fault; }
}; Scalar TryInvalM(const int& N1, const int& I1, Vector Value, const double& Quan1, const double& Quan2,
    const double& Quan3, double& MuY, double& SigmaY, double& Shift, double& Mult)
{// To try initial value for Shift = Value[i], and ML
    Shift = Value[I1];
    FunJo1 = Quan1 - Shift;				FunJo2 = Quan2 - Shift;				FunJo3 = Quan3 - Shift;
    Doub Tol = 3.0e-8, MaxA = abs(FunJo3), B1 = abs(FunJo1);
    if (B1 > MaxA)	MaxA = B1;
    int NIter = 0;
    if (MultJo(0.0) == 0.0)		Mult = 1.0;
    else Mult = zbrent(MultJo, 0.0, MaxA, Tol, NIter);
    if (NIter < 0)	Mult = 1.0;
    FunJoUML funJoU(N1, Value);
    Vector p(2), grad(2);		p[0] = Shift;		p[1] = Mult;
    Doub LogLik = funJoU(p, grad);			    MuY = funJoU.MuJo;	SigmaY = funJoU.StDevJo;
    return LogLik;
}
int InValM(const int& N1, Vector Value, double& MuY, double& SigmaY, double& ShiftY, double& MultY)
{// Initial Johnson parameters from quantiles trying to maximize LogLikelihood
// Computation of quantiles
    if (N1 < 5) return -1;
    double XN = static_cast<double>(N1);
    int N0 = static_cast<int>(XN / 2.0);
    double XN0 = static_cast<double>(N0), Quan1, Quan2, Quan3, Zn = 0.524, Pn = 0.699860730068;		// Pn = CuMNHart(Zn)
    if (XN0 == XN / 2) Quan2 = (Value[N0] + Value[N0 - 1]) / 2.0;
    else Quan2 = Value[N0];
    if (N1 > 100) {
        Zn = 0.85; Pn = 0.802337456877;
    }
    double XM = XN * Pn - 0.5;
    int M = static_cast<int>(XM);
    Quan3 = Value[M];		Quan1 = Value[N1 - M];
    // Searching for initial values 
    Scalar ScoreY = TryInvalM(N1, 0, Value, Quan1, Quan2, Quan3, MuY, SigmaY, ShiftY, MultY);
    for (int i = 1; i < N1; i++)
    {
        Scalar Mu, Sigma, Shift, Mult, Score = TryInvalM(N1, i, Value, Quan1, Quan2, Quan3, Mu, Sigma, Shift, Mult);
        if (Score < ScoreY) {
            MuY = Mu;	SigmaY = Sigma;		ShiftY = Shift;		MultY = Mult;
            ScoreY = Score;
        }
    }
    return 0;
}
bool ReadVectorI(std::string file_name, int& num_data, Vector& data_read)
{//Input vector from file file_name and store it in data_read 
    std::ifstream data_file(file_name);					//Create an object for input file file_name
    if (!data_file.good()) return false;
    std::istream_iterator<int> start(data_file), end;	//Iterators to read input file
    std::vector<int> data(start, end);					//Input integers in file to vector data
    std::vector<int>::iterator p = data.begin();		//Pointer to the first integer in data
    num_data = *p;										//First integer is number of data
    data_read.resize(num_data);						    //Space for data
    for (int i = 0; i < num_data; ++i)
    {
        ++p;	data_read[i] = static_cast<Scalar>(*p);	// Next value
    }
    return true;
}
bool ReportVector(Vector& data_read)
{
    Eigen::EigenBase<Vector>::Index size = data_read.size();
    if (size == 0) std::cout << "Vector has no data" << std::endl;
    else
    {
        std::cout << "Vector has " << size << " values" << std::endl;
        for (int i = 0; i < size; ++i)
        {
            std::cout <<  data_read[i] << " ";
        }
        std::cout << std::endl;
    }
    return true;
}
bool input_data(int& N1, double& Mu, double& Sigma, double& a, double& b)
{
    std::cout << "Enter sample size " << std::endl; std::cin >> N1;
    std::cout << "Enter Mu " << std::endl;          std::cin >> Mu;
    std::cout << "Enter Sigma " << std::endl;       std::cin >> Sigma;
    std::cout << "Enter Shift " << std::endl;       std::cin >> a;
    std::cout << "Enter Multiplier " << std::endl;  std::cin >> b;
    return true;
}
/*
int main()
{// For ML
    std::ifstream indata;					            //Create an object for input file file_name
    int ErrorCode = 0, N1 = 0, N0 = 0;                           double Mu, Sigma, a, b;
    indata.open("SimDataL.txt");
    if (!indata) { // file couldn't be opened
        std::cout << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    indata >> N0;   indata >> N1;        indata >> Mu;        indata >> Sigma;     indata >> a;     indata >> b;
    indata.close();
    std::cout << " N0 = " << N0 << " N1 = " << N1 << " Mu = " << Mu << " Sigma = " << Sigma << " Shift = " << a << " Mult = " << b << std::endl;
    if ((N1 < 5) || (Sigma <= 0.0) || (b <= 0.0)) { std::cout << "Data are not valid " << std::endl; return -2; }
    Vector DataVal;
    DataVal.resize(N1);	    if (!IniciaMT()) return 2;

    for (int j = 0; j <= N0; j++)	for (int i = 0; i < N1; i++)	DataVal[i] = JohUMT(Mu, Sigma, a, b);
    if (!ReportVector(DataVal)) { std::cout << "Error reporting data " << std::endl; return -3; }
    std::sort(DataVal.begin(), DataVal.end());
    LBFGSBParam<Scalar> param;                                                  // Default parameters to control de L-BFGS-B algorithm
    param.max_iterations = 100;         // Default is infinite
    param.max_linesearch = 100;         // Default is 20
    param.epsilon = Scalar(1e-8);       // Default is 1e-6
    param.epsilon_rel = Scalar(1e-8);   // Default is 1e-6
    param.past = 1;                     // Default is 1
    LBFGSBSolver<Scalar> solver(param);
    FunJoUML fun(N1, DataVal);                                                    // For ML
    // Variable bounds
    Vector lb = Vector::Constant(2, 0.0);                                      // For ML
    Vector ub = Vector::Constant(2, std::numeric_limits<double>::infinity());  // For ML
    lb[0] = -std::numeric_limits<double>::infinity();                                   // First variable is unbounded
    Vector x = Vector::Constant(2, 0.0);                                       // For ML
    Scalar MuJo, StDevJo;
    // Initial values 
    if (InValM(N1, DataVal, MuJo, StDevJo, x[0], x[1])) { std::cout << "Error getting initial values " << std::endl; ErrorCode = 4; return 4; };
    std::cout << "Initial Mu = " << MuJo << " Initial StDev = " << StDevJo << " Initial Shift = " << x[0] << " Initial Mult = " << x[1] <<  std::endl; // std::cin.get();
    Scalar fx; 
    Vector grad = Vector::Constant(2, 0.0);     // fx = fun(x, grad);
    int niter = solver.minimize(fun, x, fx, lb, ub);
    std::cout << niter << " iterations" << std::endl;
    MuJo = fun.getMuJo();       StDevJo = fun.getStDevJo();
    std::cout << " Mu = " << MuJo << " StDev = " << StDevJo << " Shift = " << x[0] << " Multiplier = " << x[1] <<  std::endl;
    std::cout << "f(x) = " << fx << std::endl;
    std::cout << "grad = " << solver.final_grad().transpose() << std::endl;
    std::cout << "projected grad norm = " << solver.final_grad_norm() << std::endl;
    if (fun.getFault()) { std::cout << "Error evaluating objective function " << std::endl; ErrorCode = 5; return 5; };
    // std::cin.get();
    return 0;
}
*/
//*
int main()
{// For Cramer von Mises,
    std::ifstream indata;					            //Create an object for input file file_name
    int N0 = 0, N1 = 0;                           double Mu, Sigma, a, b;
    indata.open("SimDataL.txt");
    if (!indata) { // file couldn't be opened
        std::cout << "Error: file could not be opened" << std::endl;
        exit(-1);
    }
    indata >> N0;   indata >> N1;        indata >> Mu;        indata >> Sigma;     indata >> a;     indata >> b;
    indata.close();
    std::cout << " N0 = " << N0 << " N1 = " << N1 << " Mu = " << Mu << " Sigma = " << Sigma << " Shift = " << a << " Mult = " << b << std::endl;
    if ((N1 < 5) || (Sigma <= 0.0) || (b <= 0.0)) { std::cout << "Data are not valid " << std::endl; return -2; }
    Vector DataVal;
    DataVal.resize(N1);	    if (!IniciaMT()) return -3;
    for (int j = 0; j <= N0; j++)	for (int i = 0; i < N1; i++)	DataVal[i] = JohUMT(Mu, Sigma, a, b);
    if (!ReportVector(DataVal)) { std::cout << "Error reporting data " << std::endl; return -3; }
    std::sort(DataVal.begin(), DataVal.end());
    // Parameters for LBFGSBP
    LBFGSBParam<Scalar> param;            // Default parameters to control de L-BFGS-B algorithm
    param.m = 6;                          // Default is 6
    param.max_iterations = 0;             // Default is infinite (0)
    param.max_linesearch = 20;             // Default is 20
    param.epsilon = Scalar(1e-6);          // Default is 1e-6
    param.epsilon_rel = Scalar(1e-6);      // Default is 1e-6
    param.min_step = Scalar(1e-20);         // Default is 1e-20
    param.ftol = Scalar(1e-4);               // Default is 1e-4
    param.delta = Scalar(1e-10);             // Default is 1e-10
    LBFGSBSolver<Scalar> solver(param);
    FunJoUCvM fun(N1, DataVal);                                                         // For CvM
    // Variable bounds
    Vector lb = Vector::Constant(4, 0.0);                                               // For CvM
    Vector ub = Vector::Constant(4, std::numeric_limits<double>::infinity());           // For CvM
    lb[0] = -std::numeric_limits<double>::infinity();                                   // First variable is unbounded
    lb[2] = -std::numeric_limits<double>::infinity();                                   // Third variable is unbounded for CvM
    Vector x = Vector::Constant(4, 0.0);                                                // For CvM
    Scalar MuJo, StDevJo;
    // Initial values 
    if (InValC(N1, DataVal, MuJo, StDevJo, x[0], x[1]) < 0) { std::cout << "Error getting initial values " << std::endl; return -4; };
    std::cout << "Initial Mu = " << MuJo << " Initial StDev = " <<  StDevJo << " Initial Shift = " << x[0] << " Initial Mult = " << x[1] <<  std::endl;
    x[2] = MuJo;        x[3]  = StDevJo;                                     // For CvM
    double fx;
    int niter = solver.minimize(fun, x, fx, lb, ub);
    std::cout << niter << " iterations" << std::endl;
    std::cout << " Mu = " << x[2] << " StDev = " << x[3] << " Shift = " << x[0] << " Multiplier = " << x[1] << std::endl;
    std::cout << "f(x) = " << fx << std::endl;
    if ((x[1] > 0) && (x[3] > 0))
    {
        double *Cum = new double[N1], *Value = new double[N1];
        for (int i = 0; i < N1; ++i)
        {
            Value[i] = DataVal[i];
            Cum[i] = JohUAc((Value[i] - x[0]) / x[1], x[2], x[3]);  // sample cdf 
        }
        double KS = KolStat(N1, Value, Cum), LogVer = LogVerJoU(N1, Value, x[2], x[3], x[0], x[1]), CvM = CMtest(N1, Cum);
        std::cout << " KS = " << KS << " LL = " << LogVer << " CvM = " << CvM;
    
    }
    std::cout << " grad = " << solver.final_grad().transpose() << std::endl;
    std::cout << "Projected grad norm = " << solver.final_grad_norm() << std::endl;
    std::cout << "Error code = " << solver.CodeE << std::endl;
    if (fun.getFault()) { std::cout << "Error evaluating objective function " << std::endl; return -5; };
    //std::cin.get();
    return 0;
}
// */
/*
class Rosenbrock
{
private:
    int n;
public:
    Rosenbrock(int n_) : n(n_) {}
    double operator()(const Vector& x, Vector& grad)
    {
        double fx = (x[0] - 1.0) * (x[0] - 1.0);
        grad[0] = 2 * (x[0] - 1) + 16 * (x[0] * x[0] - x[1]) * x[0];
        for (int i = 1; i < n; i++)
        {
            fx += 4 * std::pow(x[i] - x[i - 1] * x[i - 1], 2);
            if (i == n - 1)
            {
                grad[i] = 8 * (x[i] - x[i - 1] * x[i - 1]);
            }
            else {
                grad[i] = 8 * (x[i] - x[i - 1] * x[i - 1]) + 16 * (x[i] * x[i] - x[i + 1]) * x[i];
            }
        }
        return fx;
    }
};
int main()
{
    const int n = 25;
    LBFGSBParam<Scalar> param;
    LBFGSBSolver<Scalar> solver(param);
    Rosenbrock fun(n);

    // Variable bounds
    Vector lb = Vector::Constant(n, 2.0);
    Vector ub = Vector::Constant(n, 4.0);
    // The third variable is unbounded
    lb[2] = -std::numeric_limits<double>::infinity();
    ub[2] = std::numeric_limits<double>::infinity();
    // Initial values
    Vector x = Vector::Constant(n, 3.0);
    // Make some initial values at the bounds
    x[0] = x[1] = 2.0;
    x[5] = x[7] = 4.0;

    double fx;
    int niter = solver.minimize(fun, x, fx, lb, ub);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;
    std::cout << "grad = " << solver.final_grad().transpose() << std::endl;
    std::cout << "projected grad norm = " << solver.final_grad_norm() << std::endl;
    std::cin.get();
    return 0;
}
*/
