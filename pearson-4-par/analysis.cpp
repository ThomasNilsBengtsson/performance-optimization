/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>

namespace Analysis
{

    std::vector<double> correlation_coefficients(std::vector<Vector> datasets)
    {
    size_t n = datasets.size();
    size_t total_pairs = (n * (n - 1)) / 2;
    std::vector<double> result(total_pairs);

    #pragma omp parallel for schedule(dynamic)
    for (size_t sample1 = 0; sample1 < n - 1; sample1++)
    {
        size_t row_start = sample1 * n - (sample1 * (sample1 + 1)) / 2 - sample1;
        
        for (size_t sample2 = sample1 + 1; sample2 < n; sample2++)
        {
            size_t idx = row_start + sample2 - 1;
            result[idx] = pearson(datasets[sample1], datasets[sample2]);
        }
    }

    return result;
    }

    double pearson(Vector vec1, Vector vec2)
    {

        unsigned size = vec1.get_size();

        double *dataV1 = vec1.get_data();
        double *dataV2 = vec2.get_data();

        double sumX = 0.0;
        double sumY = 0.0;

        for (auto i = 0; i < size; i++)
        {
            sumX += dataV1[i];
            sumY += dataV2[i];
        }

        double xMean = sumX / size;
        double yMean = sumY / size;

        double dotProduct = 0.0;
        double xMagnitudeSquare = 0.0;
        double yMagnitudeSquare = 0.0;

        for (unsigned int i = 0; i < size; i++)
        {
            double xCentered = dataV1[i] - xMean;
            double yCentered = dataV2[i] - yMean;

            dotProduct += xCentered * yCentered;
            xMagnitudeSquare += xCentered * xCentered;
            yMagnitudeSquare += yCentered * yCentered;
        }

        double x_mag = std::sqrt(xMagnitudeSquare);
        double y_mag = std::sqrt(yMagnitudeSquare);
        double r = dotProduct / (x_mag * y_mag);

        return std::max(std::min(r, 1.0), -1.0);
    }
};
