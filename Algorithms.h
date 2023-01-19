#pragma once
#include <vector>

std::vector<double> nelderMead(double objectiveFunction(std::vector<double>), const std::vector<double>& initialPoint,
                               double alpha, double gamma, double rho, double sigma, int maxIterations);

std::vector<double> GoldenRatioMethod(double f(std::vector<double>), double A, double B, double epsilon, int nMax);
std::vector<double> GoldenSectionSearch(double f(std::vector<double>), std::vector<double> a, std::vector<double> b,
                                        double epsilon, int nMax);

std::vector<double> GradientMethod(std::vector<double> x0, double epsilon, int* fcalls, int nMax);

std::vector<double> PowellMethod(double f(std::vector<double>, double), double f2(std::vector<double>, double),
                                 std::vector<double> x0, double a, double epsilon, int Nmax);
