#include "VectorUtilities.h"

#include <vector>

std::vector<double> AddVectors(std::vector<double> vector1, std::vector<double> vector2)
{
    return std::vector<double> {vector1[0]+vector2[0], vector1[1]+vector2[1]};
}

std::vector<double> AddToVector(std::vector<double> vector, double number)
{
    return std::vector<double> {vector[0] + number, vector[1] + number};
}

std::vector<double> SubtractVectors(std::vector<double> vector1, std::vector<double> vector2)
{
    return std::vector<double> {vector1[0]-vector2[0], vector1[1]-vector2[1]};
}

std::vector<double> MultiplyVector(std::vector<double> vector, double multiplier)
{
    return std::vector<double> {vector[0] * multiplier, vector[1] * multiplier};
}