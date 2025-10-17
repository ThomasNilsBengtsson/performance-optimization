/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>
#include <omp.h>

namespace Analysis
{
    std::vector<double> correlation_coefficients(std::vector<Vector> datasets)
    {
        size_t n = datasets.size();
        size_t total_pairs = (n * (n - 1)) / 2;
        std::vector<double> result(total_pairs);

        // Pre-compute all means (parallelized)
        std::vector<double> means(n);
        #pragma omp parallel for
        for (size_t i = 0; i < n; i++)
        {
            means[i] = datasets[i].mean();
        }

        // Parallelize the outer loop only (remove collapse(2))
        #pragma omp parallel for schedule(dynamic)
        for (size_t sample1 = 0; sample1 < n - 1; sample1++)
        {
            // Calculate starting index for this row in the result vector
            size_t row_start = sample1 * n - (sample1 * (sample1 + 1)) / 2 - sample1;
            
            for (size_t sample2 = sample1 + 1; sample2 < n; sample2++)
            {
                size_t idx = row_start + sample2 - 1;
                result[idx] = pearson(datasets[sample1], datasets[sample2], 
                                     means[sample1], means[sample2]);
            }
        }

        return result;
    }

    double pearson(Vector vec1, Vector vec2, double x_mean, double y_mean)
    {
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
}