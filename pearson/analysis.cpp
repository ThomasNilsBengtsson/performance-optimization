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
        std::vector<double> result{};

        for (auto sample1{0}; sample1 < datasets.size() - 1; sample1++)
        {
            for (auto sample2{sample1 + 1}; sample2 < datasets.size(); sample2++)
            {
                auto corr{pearson(datasets[sample1], datasets[sample2])};
                result.push_back(corr);
            }
        }

        return result;
    }

    double pearson(Vector vec1, Vector vec2)
    {
        auto x_mean{vec1.mean()};
        auto y_mean{vec2.mean()};

        double numerator = 0.0;
        double x_sum_sq = 0.0;
        double y_sum_sq = 0.0;

        for (unsigned i = 0; i < vec1.get_size(); i++)
        {
            double x_diff = vec1[i] - x_mean;
            double y_diff = vec2[i] - y_mean;

            numerator += x_diff * y_diff;
            x_sum_sq += x_diff * x_diff;
            y_sum_sq += y_diff * y_diff;
        }

        double denominator = std::sqrt(x_sum_sq * y_sum_sq);
        double r = numerator / denominator;

        return std::max(std::min(r, 1.0), -1.0);
    }
};
