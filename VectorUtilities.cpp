#include "VectorUtilities.h"

#include <iostream>
#include <vector>

using namespace std;

std::vector<double> AddVectors(std::vector<double> vector1, std::vector<double> vector2)
{
    vector1[0] += vector2[0];
    vector1[1] += vector2[1];

    return vector1;
}

std::vector<double> AddToVector(std::vector<double> vector, double number)
{
    vector[0] += number;
    vector[1] += number;

    return vector;
}

std::vector<double> SubtractVectors(std::vector<double> vector1, std::vector<double> vector2)
{
    vector1[0] -= vector2[0];
    vector1[1] -= vector2[1];

    return vector1;
}

std::vector<double> MultiplyVector(std::vector<double> vector, double multiplier)
{
    vector[0] *= multiplier;
    vector[1] *= multiplier;

    return vector;
}

void PrintVector(std::vector<double> vector)
{
    std::cout << "Vector {" << vector[0] << ", " << vector[1] << "}" << endl;
}
