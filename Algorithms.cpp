#include "Algorithms.h"

#include <algorithm>
#include <vector>
#include "VectorUtilities.h"
#include <cmath>
#include <complex>
#include <iostream>

using namespace std;

vector<double> nelderMead(double objectiveFunction(std::vector<double>), const std::vector<double>& initialPoint,
                          double alpha, double gamma, double rho, double sigma, int maxIterations)
{
    // Inicjalizujemy sympleks jako tablicę wektorów
    std::vector<std::vector<double>> simplex;
    simplex.push_back(initialPoint); // Punkt startowy jest pierwszym wierzchołkiem sympleksu
    int n = initialPoint.size(); // Określamy liczbę zmiennych (rozmiar wektora)

    // Dodajemy pozostałe wierzchołki sympleksu, które są równe punktowi startowemu, ale z jednym wymiarem przesuniętym o wartość delta
    double delta = 0.1; // Możesz zmienić tę wartość na inne
    for (int i = 0; i < n; i++)
    {
        std::vector<double> newPoint = initialPoint;
        newPoint[i] += delta;
        simplex.push_back(newPoint);
    }

    // Główna pętla iteracyjna
    for (int iteration = 0; iteration < maxIterations; iteration++)
    {
        // Posortuj wierzchołki sympleksu według wartości funkcji celu
        std::sort(simplex.begin(), simplex.end(), [&](const std::vector<double>& a, const std::vector<double>& b)
        {
            return objectiveFunction(a) < objectiveFunction(b);
        });

        // Wylicz średnią arytmetyczną wszystkich wierzchołków sympleksu, z wyjątkiem najgorszego
        std::vector<double> mean(n, 0.0);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                mean[j] += simplex[i][j] / n;
            }
        }

        // Oblicz nowy punkt Xr przez odbicie najgorszego wierzchołka sympleksu od średniej arytmetycznej
        std::vector<double> xr(n);
        for (int i = 0; i < n; i++)
        {
            xr[i] = mean[i] + alpha * (mean[i] - simplex[n][i]);
        }

        // Porównaj wartość funkcji celu dla nowego punktu Xr z wartością dla najgorszego wierzchołka sympleksu
        double fxr = objectiveFunction(xr);
        double fworst = objectiveFunction(simplex[n]);
        if (fxr < fworst)
        {
            // Jeśli nowy punkt jest lepszy, to zamień go z najgorszym wierzchołkiem sympleksu
            simplex[n] = xr;

            // Sprawdź, czy nowy punkt jest lepszy od najlepszego wierzchołka sympleksu
            double fbest = objectiveFunction(simplex[0]);
            if (fxr < fbest)
            {
                // Jeśli tak, to oblicz nowy punkt Xe przez dalsze odbicie od średniej arytmetycznej
                std::vector<double> xe(n);
                for (int i = 0; i < n; i++)
                {
                    xe[i] = mean[i] + gamma * (xr[i] - mean[i]);
                }

                // Porównaj wartość funkcji celu dla nowego punktu Xe z wartością dla najlepszego wierzchołka sympleksu
                double fxe = objectiveFunction(xe);
                if (fxe < fbest)
                {
                    // Jeśli Xe jest lepsze, to zamień najgorszy wierzchołek sympleksu na Xe
                    simplex[n] = xe;
                }
            }
        }
        else
        {
            // Jeśli nowy punkt Xr jest gorszy od najgorszego wierzchołka sympleksu, to oblicz punkt Xc przez kontrakcję
            std::vector<double> xc(n);
            for (int i = 0; i < n; i++)
            {
                xc[i] = mean[i] + rho * (simplex[n][i] - mean[i]);
            }

            // Porównaj wartość funkcji celu dla nowego punktu Xc z wartością dla najgorszego wierzchołka sympleksu
            double fxc = objectiveFunction

                (xc);
            if (fxc < fworst)
            {
                // Jeśli Xc jest lepsze, to zamień najgorszy wierzchołek sympleksu na Xc
                simplex[n] = xc;
            }
            else
            {
                // Jeśli Xc jest gorsze, to zmniejsz każdy wierzchołek sympleksu o sigma * (średnia arytmetyczna - wierzchołek)
                for (int i = 1; i < n + 1; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        simplex[i][j] = simplex[i][j] + sigma * (mean[j] - simplex[i][j]);
                    }
                }
            }
        }
    }

    // Po zakończeniu pętli iteracyjnej zwróć najlepszy wierzchołek sympleksu jako wynik
    return simplex[0];
}

vector<double> GoldenRatioMethod(double f(vector<double>), double A, double B, double epsilon, int nMax)
{
    double alfa = (sqrt(5) - 1) / 2;

    auto a = new vector<double>[100];
    auto b = new vector<double>[100];
    auto c = new vector<double>[100];
    auto d = new vector<double>[100];

    a[0] = vector<double>{0, 0};
    b[0] = vector<double>{A, B};
    c[0] = SubtractVectors(b[0], MultiplyVector(SubtractVectors(b[0], a[0]), alfa));
    d[0] = AddVectors(a[0], MultiplyVector(SubtractVectors(b[0], a[0]), alfa));

    int i;
    for (i = 0; VectorLength(SubtractVectors(b[i], a[i])) > epsilon; ++i)
    {
        if (f(c[i]) < f(d[i]))
        {
            a[i + 1] = a[i];
            b[i + 1] = d[i];
            c[i + 1] = SubtractVectors(b[i + 1], MultiplyVector(SubtractVectors(b[i + 1], a[i + 1]), alfa));
            d[i + 1] = c[i];
        }
        else
        {
            a[i + 1] = c[i];
            b[i + 1] = b[i];
            c[i + 1] = d[i];
            d[i + 1] = SubtractVectors(a[i + 1], MultiplyVector(SubtractVectors(b[i + 1], a[i + 1]), alfa));
        }
    }

    return MultiplyVector(AddVectors(a[i], b[i]), 0.5);
}

constexpr double GoldenRatio = 1.61803398874989;
constexpr double InvGoldenRatio = 0.61803398874989;

vector<double> GradientMethod(vector<double> x0, double epsilon, int* fcalls, int nMax)
{
    int i = 0;
    double d;
    vector<double> h;
    auto x = new vector<double>[100];

    for (int i = 0; VectorLength(SubtractVectors(x[i], x[i - 1])) < epsilon || *fcalls > nMax; ++i)
    {
        d = VectorLength(x[i]);
        h = MultiplyVector(vector<double>{1, 1}, x0[0]);

        x[i + 1] = AddVectors(x[i], MultiplyVector(h, d));
    }

    return x[i];
}

std::vector<double> GoldenSectionSearch(double f(std::vector<double>), std::vector<double> a, std::vector<double> b,
                                        double epsilon, int nMax)
{
    std::vector<double> d = MultiplyVector(SubtractVectors(b, a), InvGoldenRatio);
    std::vector<double> c = AddVectors(a, d);
    double fc = f(c);
    std::vector<double> e = MultiplyVector(SubtractVectors(b, a), GoldenRatio);
    double fd = f(AddVectors(a, e));

    while (VectorLength(SubtractVectors(b, a)) > epsilon)
    {
        if (fc < fd)
        {
            b = AddVectors(a, e);
            e = d;
            d = MultiplyVector(SubtractVectors(b, a), InvGoldenRatio);
            fd = fc;
            c = AddVectors(a, d);
            fc = f(c);
        }
        else
        {
            a = c;
            c = AddVectors(a, e);
            e = d;
            d = MultiplyVector(SubtractVectors(b, a), InvGoldenRatio);
            fc = fd;
            fd = f(AddVectors(a, e));
        }
    }

    return MultiplyVector(AddVectors(a, b), 0.5);
}

vector<double> Newton(double f(vector<double>), vector<double> x0, double epsilon, double alpha)
{
    int n = 2;

    while (true)
    {
        vector<double> grad(n);
        for (int i = 0; i < n; i++)
        {
            vector<double> x_perturbed = x0;
            x_perturbed[i] += epsilon;
            grad[i] = (f(x_perturbed) - f(x0)) / epsilon;
        }
        PrintVector(grad);

        vector<vector<double>> hessian(n, vector<double>(n));
        int i = 0;
        for (i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                vector<double> x_perturbed = x0;
                x_perturbed[i] += epsilon;
                x_perturbed[j] += epsilon;
                hessian[i][j] = (f(x_perturbed) - f(x0) - grad[i] * epsilon - grad[j] * epsilon) / (epsilon * epsilon);
            }
        }

        vector<double> x1 = SubtractVectors(x0, MultiplyVector(AddVectors(hessian[i], grad), alpha));

        if (VectorLength(SubtractVectors(x1, x0)) < epsilon)
        {
            return x1;
        }
        x0 = x1;
    }
}

vector<double> PowellMethod(double f(vector<double>, double), double f2(vector<double>, double), vector<double> x0,
                            double a, double epsilon, int Nmax)
{
    int i = 0;
    vector<vector<double>> d(x0.size());
    vector<vector<double>> p(x0.size() + 1);
    vector<double> x(x0);
    vector<double> h(x0.size());

    p[0] = x0;

    for (int i = 0; i < x0.size(); i++)
    {
        d[i].resize(x0.size());
        d[i][i] = 1;
    }

    while (i < Nmax)
    {
        int imax = 0;
        double fmax = f(x, a);
        for (int j = 0; j < x0.size(); j++)
        {
            h = MultiplyVector(d[j], f2(x, a));
            p[j + 1] = AddVectors(x, h);
            double fj = f(p[j + 1], a);
            if (fj > fmax)
            {
                imax = j;
                fmax = fj;
            }
        }
        if (fmax - f(x, a) < epsilon)
            return x;

        x = p[imax + 1];
        d[imax] = SubtractVectors(p[imax + 1], p[0]);
        for (int j = 0; j < imax; j++)
            d[j] = d[j + 1];
        
        h = SubtractVectors(p[imax + 1], p[0]);
        d[x0.size() - 1] = MultiplyVector(h, 1.0 / VectorLength(h));
        i++;
    }
    return x;
}
