#pragma once
#include <vector>

std::vector<double> nelderMead(const std::vector<double>& initialPoint, double alpha, double gamma, double rho,
                               double sigma, int maxIterations);
