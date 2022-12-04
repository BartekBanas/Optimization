#include "opt_alg.h"
#include <iostream>
#include <cmath>
#include "VectorUtilities.h"

using namespace std;

#define M_PI 3.14159265358979323846

double Function1(double x);
double Function2(vector<double> x);
double Function3(double x1, double x2);

double Fibonacci(double function(double), double a, double b, double precision);
double lagrange(double aInput, double bInput, double eps, double gamma, int Nmax);
vector<double> HookeJeeves(double function(vector<double>), vector<double> x, double step, double alfa, double epsilon,
                           int nMax);
vector<double> Trying(vector<double> x, double step, double function(vector<double>));

int n = 2;
int fcalls = 0;

void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
    try
    {
        //lab1();
        lab2();
    }
    catch (string EX_INFO)
    {
        cerr << "ERROR:\n";
        cerr << EX_INFO << endl
            << endl;
    }

    return 0;
}

int iterations;

void lab1()
{
    double result = Fibonacci(Function1, -10, 10, 0.00001);
    cout << "The result of the Fibonacci method: " << endl;
    cout << "x = " << result << endl;
    cout << "y = " << Function1(result) << endl;
    cout << "Iterations: " << iterations << endl;

    double resultLagrange = lagrange(-100, 100, 0.001, 1e-7, 1000);
    cout << "The result of the lagrange method: " << endl;
    cout << "x = " << resultLagrange << endl;
    cout << "y = " << Function1(resultLagrange) << endl;
    // cout << "Iterations: " << iterations << endl;
}

void lab2()
{
    double step = 0.5;
    vector<double> x = {1, 0.5};
    PrintVector(x);
    double alfa = 0.5;
    double epsilon = 1e-3;
    int Nmax = 1000;

    vector<double> result = (HookeJeeves(Function2, x, step, alfa, epsilon, Nmax));
    
    PrintVector(result);
    cout << "f(x) = " << Function2(result) << endl;

    cout << "Fcalls: " << fcalls << endl;
}

void lab3()
{
}

void lab4()
{
}

void lab5()
{
}

void lab6()
{
}

double Function1(double x)
{
    return -cos(0.1 * x) * exp(-pow(0.1 * x - 2 * M_PI, 2)) +
        0.002 * pow(0.1 * x, 2);
}

double Function2(vector<double> x)
{
    fcalls++;
    return x[0] * x[0] + x[1] * x[1] - cos(2.5 * M_PI * x[0]) - cos(2.5 * M_PI * x[1]) + 2;
}

double Function3(vector<double> x)
{
    fcalls++;
    return sin(M_PI * sqrt(pow(x[0] / M_PI, 2) + pow(x[1] / M_PI, 2))) /
        M_PI * sqrt(pow(x[0] / M_PI, 2) + pow(x[1] / M_PI, 2));
}

double Fibonacci(double function(double), double a, double b, double precision)
{
    double* fibonacci = new double[100]{0, 1};

    vector<double> x1(5, 5);


    vector<double> x2;
    x2.push_back(6.0);
    x2.push_back(2);
    PrintVector(x1);
    PrintVector(x2);

    //vector<double> x3 = vektorek.AddVectors(x1, x2);


    vector<double> x3 = AddVectors(x1, x2);
    PrintVector(x3);

    int k = 1;

    while (fibonacci[k] <= (b - a) / precision)
    {
        fibonacci[k + 1] = fibonacci[k] + fibonacci[k - 1];
        k++;
    }

    double* As = new double[k];
    double* bs = new double[k];
    double* cs = new double[k];
    double* ds = new double[k];

    As[0] = a;
    bs[0] = b;
    cs[0] = b - fibonacci[k - 1] / fibonacci[k] * (bs[0] - As[0]);
    ds[0] = As[0] + bs[0] - cs[0];

    int i = 0;

    for (i = 0; i < k - 3; i++)
    {
        iterations++;
        if (function(cs[i]) < function(ds[i]))
        {
            As[i + 1] = As[i];
            bs[i + 1] = ds[i];
        }
        else
        {
            bs[i + 1] = bs[i];
            As[i + 1] = cs[i];
        }

        cs[i + 1] = bs[i + 1] - fibonacci[k - i - 2] / fibonacci[k - i - 1] * (bs[i + 1] - As[i + 1]);

        ds[i + 1] = As[i + 1] + bs[i + 1] - cs[i + 1];
    }

    return cs[i];
}

double lagrange(double aInput, double bInput, double eps, double gamma, int Nmax)
{
    int i = 0;
    double a = aInput;
    double b = bInput;
    double c = (b + a) / 2;
    double l, m;
    double fcalls = 0;
    double d0, d = 2;
    do
    {
        d0 = d;
        fcalls++;
        l = Function1(a) * (b * b - c * c) + Function1(b) * (c * c - a * a) + Function1(c) * (a * a - b * b);
        m = Function1(a) * (b - c) + Function1(b) * (c - a) + Function1(c) * (a - b);

        if (m <= 0)
            return 404;

        d = 0.5 * l / m;

        if (a < d && c > d)
        {
            if (Function1(d) < Function1(c))
            {
                a = a;
                b = c;
                c = d;
            }
            else
            {
                a = d;
                b = b;
            }
        }
        else
        {
            if (c < d && b > d)
                if (Function1(d) < Function1(c))
                {
                    a = c;
                    c = d;
                }
                else
                {
                    a = a;

                    b = d;
                }
            else
                return 4042;
        }
        i++;
        if (fcalls > Nmax)
            return 40422;
    } // while (b - a < eps || abs(d - d0) < gamma);
    while (b - a >= eps && abs(d - d0) >= gamma);
    return d;
}

vector<double> HookeJeeves(double function(vector<double>), vector<double> x, double step, double alfa, double epsilon,
                           int nMax)
{
    vector<double> error;
    error.push_back(0.00002137);
    vector<double> xB;
    vector<double> xB_;
    do
    {
        xB = x;
        x = Trying(xB, step, function);
        if (function(x) < function(xB))
        {
            while (function(x) < function(xB))
            {
                xB_ = xB;
                xB = x;
                x = SubtractVectors(MultiplyVector(xB, 2), xB_);
                x = Trying(x, step, function);
                if (fcalls > nMax)
                {
                    return error;
                }
            }
            
            x = xB;
        }
        else
        {
            step = alfa * step;
        }
        if (fcalls > nMax)
        {
            return error;
        }
    }
    while (step < epsilon);
    return xB;
}

vector<double> Trying(vector<double> x, double step, double function(vector<double>))
{
    // = {(1, 0), (0, 1), (-1, 0), (0, -1)}
    vector<double> eJ[4];
    eJ[0].push_back(1);
    eJ[0].push_back(0);
    eJ[1].push_back(0);
    eJ[1].push_back(1);
    eJ[2].push_back(-1);
    eJ[2].push_back(0);
    eJ[3].push_back(0);
    eJ[3].push_back(1);

    vector<double> v1;
    vector<double> v2;

    for (int i = 1; i <= n; ++i)
    {
        v1 = MultiplyVector(eJ[i], step);
        v2 = AddVectors(x, v1);

        if (function(v2) < function(x))
        {
            x = AddVectors(x, MultiplyVector(eJ[i], step));
        }
        else
        {
            if (function(SubtractVectors(x, MultiplyVector(eJ[i], step))) < function(x))
            {
                x = SubtractVectors(x, MultiplyVector(eJ[i], step));
            }
        }
    }

    return x;
}

// vector<double> NelderMeadMethod(vector<double> x0, double s, double α, double β, double γ, double δ, double ε, int nMax)
// {
//     auto p = new vector<double>[n];
//     p[0] = x0;
//     for (int i = 1; i < n; ++i)
//     {
//         p[i] = p[0] + s * 
//     }
// }