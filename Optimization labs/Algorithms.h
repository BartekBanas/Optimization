#pragma once
#include <vector>

double objectiveFunction(const std::vector<double>& x);

std::vector<double> nelderMead(const std::vector<double>& initialPoint, double alpha, double gamma, double rho,
                               double sigma, int maxIterations);
