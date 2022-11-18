#pragma once
#include <vector>

class VectorUtilities
{
public:
    std::vector<double> AddVectors(std::vector<double> vector1, std::vector<double> vector2);
    std::vector<double> SubtractVectors(std::vector<double> vector1, std::vector<double> vector2);
    std::vector<double> MultiplyVector(std::vector<double> vector, double multiplier);
};
