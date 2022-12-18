#include "Algorithms.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

// Funkcja celu, którą chcemy maksymalizować/minimalizować
double objectiveFunction(const std::vector<double>& x) 
{
    // Tutaj możesz zdefiniować własną funkcję celu
    // Pamiętaj, że x to wektor argumentów funkcji
    // Możesz go używać jak dowolnej tablicy
    return std::pow(x[0], 2) + std::pow(x[1], 2);
}

// Metoda sympleksu Neldera-Meada
std::vector<double> nelderMead(const std::vector<double>& initialPoint, double alpha, double gamma, double rho,
                               double sigma, int maxIterations) 
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
