o use LBFGS++, one needs to first define a functor to represent the multivariate function to be minimized. It should return the objective function value on a vector x and overwrite the vector grad with the gradient evaluated on x. For example we could define the Rosenbrock function in the following way:

#include <Eigen/Core>
#include <iostream>
#include <LBFGS.h>

using Eigen::VectorXd;
using namespace LBFGSpp;

class Rosenbrock
{
private:
    int n;
public:
    Rosenbrock(int n_) : n(n_) {}
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        double fx = 0.0;
        for(int i = 0; i < n; i += 2)
        {
            double t1 = 1.0 - x[i];
            double t2 = 10 * (x[i + 1] - x[i] * x[i]);
            grad[i + 1] = 20 * t2;
            grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
            fx += t1 * t1 + t2 * t2;
        }
        return fx;
    }
};
Then we just need to set up parameters, create solver object, provide initial guess, and then run the minimization function.

int main()
{
    const int n = 10;
    // Set up parameters
    LBFGSParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 100;

    // Create solver and function object
    LBFGSSolver<double> solver(param);
    Rosenbrock fun(n);

    // Initial guess
    VectorXd x = VectorXd::Zero(n);
    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(fun, x, fx);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    return 0;
}
The example can then be compiled and run.

$ g++ -I/path/to/eigen -I/path/to/lbfgspp/include -O2 example.cpp
$ ./a.out
23 iterations
x =
1 1 1 1 1 1 1 1 1 1
f(x) = 1.87948e-19
You can also use a different line search algorithm by providing a second template parameter to LBFGSSolver. For example, the code below illustrates the bracketing line search algorithm (contributed by @DirkToewe).

int main()
{
    const int n = 10;
    // Set up parameters
    LBFGSParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 100;

    // Create solver and function object
    LBFGSSolver<double, LineSearchBracketing> solver(param);
    Rosenbrock fun(n);

    // Initial guess
    VectorXd x = VectorXd::Zero(n);
    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(fun, x, fx);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    return 0;
}
Box-constrained Problem
If the parameters to be optimized have simple bounds, then the L-BFGS-B solver class LBFGSBSolver can be used. The code is very similar to that of LBFGSSolver. Below is the same Rosenbrock example, but we require that all variables should be between 2 and 4.

#include <Eigen/Core>
#include <iostream>
#include <LBFGSB.h>  // Note the different header file

using Eigen::VectorXd;
using namespace LBFGSpp;

class Rosenbrock
{
private:
    int n;
public:
    Rosenbrock(int n_) : n(n_) {}
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        double fx = 0.0;
        for(int i = 0; i < n; i += 2)
        {
            double t1 = 1.0 - x[i];
            double t2 = 10 * (x[i + 1] - x[i] * x[i]);
            grad[i + 1] = 20 * t2;
            grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
            fx += t1 * t1 + t2 * t2;
        }
        return fx;
    }
};

int main()
{
    const int n = 10;
    // Set up parameters
    LBFGSBParam<double> param;  // New parameter class
    param.epsilon = 1e-6;
    param.max_iterations = 100;

    // Create solver and function object
    LBFGSBSolver<double> solver(param);  // New solver class
    Rosenbrock fun(n);

    // Bounds
    VectorXd lb = VectorXd::Constant(n, 2.0);
    VectorXd ub = VectorXd::Constant(n, 4.0);

    // Initial guess
    VectorXd x = VectorXd::Constant(n, 3.0);

    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(fun, x, fx, lb, ub);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    return 0;
}