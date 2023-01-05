#include "opt_alg.h"
#include <iostream>
#include <cmath>
#include "VectorUtilities.h"
#include "Algorithms.h"

using namespace std;

#define M_PI 3.14159265358979323846

double Function1(double x);
double Function2(vector<double> x);
double funtion2_real(vector<double> x);
double Function3(vector<double> x);
double Function4(vector<double> x);

double Fibonacci(double function(double), double a, double b, double precision);
double lagrange(double aInput, double bInput, double eps, double gamma, int Nmax);
vector<double> HookeJeeves(double function(vector<double>), vector<double> x, double step, double alfa, double epsilon,
                           int nMax);
vector<double> Trying(vector<double> x, double step, double function(vector<double>));
vector<double> NelderMeadMethod(vector<double> x0, double s, double alfa, double beta, double gamma, double delta,
                                double epsilon, int nMax);
vector<double> NelderMeadMethodAi(vector<double> x0, double s, double alfa, double beta, double gamma, double delta,
                                  double epsilon, int nMax);

matrix function2realistic(matrix K, matrix alfaT, matrix empty);
matrix df(double t, matrix Y, matrix empty, matrix ud2);


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
        //lab2();
        lab3();
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
    double step = 0.25;
    vector<double> x = {1, 0.75};
    PrintVector(x);
    double alfa = 0.5;
    double epsilon = 1e-3;
    int Nmax = 1000;

    vector<double> result = HookeJeeves(Function2, x, step, alfa, epsilon, Nmax);

    PrintVector(result);
    cout << "f(x) = " << Function2(result) << endl;

    cout << "Fcalls: " << fcalls << endl;

    matrix k(2, 1, 0.5);
    matrix x0(2, 1, 0.5);
    double alphaHJ = 0.5;

    solution resultHJ = HJmethod(function2realistic, x0, step, alphaHJ, epsilon, Nmax);
    cout << resultHJ.x << endl;
    cout << resultHJ.y << endl;
    cout << resultHJ.f_calls << endl;


    matrix* simulation = solve_ode(df, 0, 0.1, 100, matrix(2, 1), matrix(0), resultHJ.x);
    cout << "HJ-a" << "Hj-w" << endl;
    for (int i = 0; i < 1000; ++i)
    {
        cout << simulation[1][0](i) << " ; ";
        cout << simulation[1][1](i) << endl;
    }
}

void lab3()
{
    vector<double> x = {1, 0.75};
    double step = 0.25;
    PrintVector(x);
    double alfa = 0.5;
    double beta = 0;
    double gamma = 0;
    double delta = 0;
    double epsilon = 1e-3;
    int Nmax = 1000;


    vector<double> result = NelderMeadMethodAi(x, step, alfa, beta, gamma, delta, epsilon, Nmax);
    PrintVector(result);

    result = NelderMeadMethod(x, step, alfa, beta, gamma, delta, epsilon, Nmax);
    PrintVector(result);
    
    result = nelderMead(Function3 ,result, alfa, gamma, delta, epsilon, Nmax);
    PrintVector(result);
}

void lab4()
{
    std::vector<double> a = {0, 0};
    std::vector<double> b = {10, 10};
    
    PrintVector(a);
    PrintVector(b);

    PrintVector(GoldenSectionSearch(Function4, a, b, 0.001, 1000));
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

matrix df(double t, matrix Y, matrix empty, matrix ud2)
{
    double mr = 1, mc = 9, l = 0.5, b = 0.5, a_ref = M_PI, w_ref = 0;
    double I = mr * l * l / 3 + mc * l * l;
    double k1 = ud2(0), k2 = ud2(1);
    double M = k1 * (a_ref - Y(0)) + k2 * (w_ref - Y(1));
    matrix dY(2, 1);
    dY(0) = Y(1);
    dY(1) = (M - b * Y(1)) / I;
    return dY;
}

matrix function2realistic(matrix K, matrix alfaT, matrix empty)
{
    matrix y;
    matrix Y0(2, 1);
    matrix* Y = solve_ode(df, 0, 0.1, 100, Y0, alfaT, K);
    int n = get_len(Y[0]);
    double a_ref = M_PI, w_ref = 0;
    y = 0;
    for (int i = 0; i < n; i++)
        y = y + 10 * pow(a_ref - Y[1](i, 0), 2) + pow(w_ref - Y[1](i, 1), 2) +
            pow(K(0) * (a_ref - Y[1](i, 0)) + K(1) * (w_ref - Y[1](i, 1)), 2);
    y = y * 0.1;

    return y;

    // vector<double> result;
    //
    // result.push_back(y[0](0, 0));
    // result.push_back(y[1](50, 1));
    //
    // // result.push_back(0);
    // // result.push_back(0);
    //
    // return result;
}

double Function3(vector<double> x)
{
    fcalls++;
    return sin(M_PI * sqrt(pow(x[0] / M_PI, 2) + pow(x[1] / M_PI, 2))) /
        M_PI * sqrt(pow(x[0] / M_PI, 2) + pow(x[1] / M_PI, 2));
}

double Function4(vector<double> x)
{
    fcalls++;
    return pow(x[0] + 2 * x[1] - 7 ,2) + pow(2 * x[0] + x[1] - 5, 2);
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

vector<double> NelderMeadMethod(vector<double> x0, double s, double alfa, double beta, double gamma, double delta,
                                double epsilon, int nMax)
{
    vector<double> e[4];
    {
        e[0].push_back(1);
        e[0].push_back(0);
        e[1].push_back(0);
        e[1].push_back(1);
        e[2].push_back(-1);
        e[2].push_back(0);
        e[3].push_back(0);
        e[3].push_back(1);
    }

    auto p = new vector<double>[n + 1];
    p[0] = x0;
    for (int i = 1; i <= n; ++i)
    {
        p[i] = AddVectors(p[0], MultiplyVector(e[i], s));
    }

    
    vector<double> pMin(2, 0), pMax(2, 0);
    double min = 0, max = 0;

    while (max < epsilon || Function3(pMin) - Function3(p[n]) < epsilon)
    {
        for (int i = 0; i <= n; ++i)    //Finding maximum and minimum values
        {
            if (Function3(p[i]) > max)
            {
                pMax = p[i];
                max = Function3(p[i]);
                cout << "Point nr 1, fcalls: " << fcalls << endl;
            }

            if (Function3(p[i]) < min)
            {
                pMin = p[i];
                min = Function3(p[i]);
            }
        }

        vector<double> p_;

        for (int i = 0; i < n; ++i)
        {
            p_ = AddVectors(p_, p[i]);
            cout << "Point nr 2, fcalls: " << fcalls << endl;
        }   cout << "Point nr 2.1, fcalls: " << fcalls << endl;
        p_ = MultiplyVector(p_, 1 / static_cast<double>(n));
        

        vector<double> pOdb = AddVectors(p_, MultiplyVector(SubtractVectors(p_, pMax), alfa));
        cout << "Point nr 2.2, fcalls: " << fcalls << endl;

        if (Function3(pOdb) < Function3(pMin))
        {
            cout << "Point nr 2.3, fcalls: " << fcalls << endl;
            vector<double> pe(2, 0);
            pe = AddVectors(p_, MultiplyVector(SubtractVectors(pOdb, p_), gamma));
            cout << "Point nr 2.5, fcalls: " << fcalls << endl;

            if (Function3(pe) < Function3(pOdb))
            {
                pMax = pe;
            }
            else
            {
                pMax = pOdb;
            }   cout << "Point nr 3, fcalls: " << fcalls << endl;
        }
        else
        {
            if (Function3(pMin) <= Function3(pOdb) < Function3(pMax))
            {
                pMax = pOdb;
            }
            else
            {
                vector<double> pz = AddVectors(p_, MultiplyVector(SubtractVectors(pMax, p_), beta));
                if (Function3(pz) >= Function3(pMax))
                {
                    for (int i = 0; i < n; ++i)
                    {
                        p[i] = MultiplyVector(AddVectors(p[i], pMin), delta);
                    }
                }
                else
                {
                    pMax = pz;
                }
            }
        }
    }


    solution solution(2, x0);
    cout << "Solutione: " << solution << endl;

    return pMin;
}

vector<double> NelderMeadMethodAi(vector<double> x0, double s, double alfa, double beta, double gamma, double delta,
                                  double epsilon, int nMax)
{
    vector<double> e[4];
    {
        e[0].push_back(1);
        e[0].push_back(0);
        e[1].push_back(0);
        e[1].push_back(1);
        e[2].push_back(-1);
        e[2].push_back(0);
        e[3].push_back(0);
        e[3].push_back(-1);
    }

    auto p = new vector<double>[n + 1];
    p[0] = x0;
    for (int i = 1; i <= n; ++i)
    {
        p[i] = AddVectors(p[0], MultiplyVector(e[i], s));
    }

    int pMin = 0, pMax = 0;
    double min = 0, max = 0;

    while (max < epsilon || Function3(p[pMin]) - Function3(p[n]) < epsilon) {
        for (int i = 0; i <= n; ++i) //Finding maximum and minimum values
        {
            if (Function3(p[i]) > max) {
                pMax = i;
                max = Function3(p[i]);
            }
            if (Function3(p[i]) < min) {
                pMin = i;
                min = Function3(p[i]);
            }
        }

        vector<double> p_(2, 0);

        for (int i = 0; i < n; ++i)
            p_ = AddVectors(p_, p[i]);
        
        p_ = MultiplyVector(p_, 1 / static_cast<double>(n));

        vector<double> pOdb = AddVectors(p_, MultiplyVector(SubtractVectors(p_, p[pMax]), alfa));

        if (Function3(pOdb) < Function3(p[pMin]))
        {
            vector<double> pe = AddVectors(p_, MultiplyVector(SubtractVectors(pOdb, p_), gamma));
            if (Function3(pe) < Function3(p[pMin])) {
                p[pMax] = pe;
            }
            else {
                p[pMax] = pOdb;
            }
        }
        else if (Function3(pOdb) < Function3(p[pMax])) {
            p[pMax] = pOdb;
        }
        else {
            vector<double> pc = AddVectors(p_, MultiplyVector(SubtractVectors(p_, p[pMax]), beta));
            if (Function3(pc) < Function3(p[pMax])) {
                p[pMax] = pc;
            }
            else {
                vector<double> pd = AddVectors(p[pMin], MultiplyVector(SubtractVectors(p[pMax], p[pMin]), delta));
                p[pMax] = pd;
            }
        }
    }
    
    return p[pMin];
}

double penalty(const vector<double>& x, double r)
{
    double sum = 0;
    for (size_t i = 0; i < x.size(); i++)
    {
        if (x[i] < 0)
        {
            sum += r * std::pow(-x[i], 2);
        }
    }
    return sum;
}

double PunishMeDaddyOut(const std::vector<double>& x, double r)
{
    double sum = 0;
    for (size_t i = 0; i < x.size(); i++)
    {
        if (x[i] < 0)
        {
            sum += r * std::pow(-x[i], 2);
        }
    }
    return sum;
}

double PunishMeDaddyIn(const std::vector<double>& x, double r)
{
    double sum = 0;
    for (size_t i = 0; i < x.size(); i++)
    {
        if (x[i] < 0)
        {
            sum += r * std::pow(-x[i], 2);
        }
    }
    return sum + Function3(x);
}