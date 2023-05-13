#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

using namespace std;

double EuclideanDistance(const vector<double>& vector1, const vector<double>& vector2)
{
    if (vector1.size() != vector2.size())
        throw invalid_argument("Vectors must have the same dimensions.");

    double sumOfSquares = 0.0;

#pragma omp parallel for reduction(+:sumOfSquares)
    for (int i = 0; i < vector1.size(); i++)
    {
        double difference = vector1[i] - vector2[i];
        sumOfSquares += difference * difference;
    }

    return sqrt(sumOfSquares);
}

double ManhattanDistance(const vector<double>& vector1, const vector<double>& vector2)
{
    if (vector1.size() != vector2.size())
        throw invalid_argument("Vectors must have the same dimensions.");

    double sumOfDifferences = 0.0;

#pragma omp parallel for reduction(+:sumOfDifferences)
    for (int i = 0; i < vector1.size(); i++)
    {
        double difference = abs(vector1[i] - vector2[i]);
        sumOfDifferences += difference;
    }

    return sumOfDifferences;
}

int main()
{
    const vector<double> vector1 = {1.5, 2.0, 3.2};
    const vector<double> vector2 = {2.0, 3.5, 1.8};

    double euclideanDist = EuclideanDistance(vector1, vector2);
    double manhattanDist = ManhattanDistance(vector1, vector2);

    cout << "Euclidean Distance: " << euclideanDist << endl;
    cout << "Manhattan Distance: " << manhattanDist << endl;

    return 0;
}
