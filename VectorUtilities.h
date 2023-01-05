#pragma once
#include <vector>

std::vector<double> AddVectors(std::vector<double> vector1, std::vector<double> vector2);
std::vector<double> AddToVector(std::vector<double> vector, double number);
std::vector<double> SubtractVectors(std::vector<double> vector1, std::vector<double> vector2);
std::vector<double> MultiplyVector(std::vector<double> vector, double multiplier);
void PrintVector(std::vector<double> vector);
double VectorLength(std::vector<double> vector)