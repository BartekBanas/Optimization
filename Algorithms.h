#pragma once
#include <vector>

std::vector<double> nelderMead(double objectiveFunction(std::vector<double>), const std::vector<double>& initialPoint,
    double alpha, double gamma, double rho,double sigma, int maxIterations);
std::vector<double> GoldenSectionSearch(double f(std::vector<double>), std::vector<double> a, std::vector<double> b, double epsilon, int nMax);