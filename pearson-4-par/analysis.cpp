/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <pthread.h>

namespace Analysis
{
    // Structure to pass data to each thread
    struct ThreadData {
        std::vector<Vector>* datasets;
        std::vector<double>* result;
        size_t start_row;
        size_t end_row;
        size_t n;
    };

    // Thread function that computes correlations for assigned rows
    void* compute_correlations_thread(void* arg)
    {
        ThreadData* data = static_cast<ThreadData*>(arg);
        
        for (size_t sample1 = data->start_row; sample1 < data->end_row; sample1++)
        {
            size_t row_start = sample1 * data->n - (sample1 * (sample1 + 1)) / 2 - sample1;
            
            for (size_t sample2 = sample1 + 1; sample2 < data->n; sample2++)
            {
                size_t idx = row_start + sample2 - 1;
                (*data->result)[idx] = pearson((*data->datasets)[sample1], 
                                               (*data->datasets)[sample2]);
            }
        }
        
        return nullptr;
    }

    std::vector<double> correlation_coefficients(std::vector<Vector> datasets, int num_threads)
    {
        size_t n = datasets.size();
        size_t total_pairs = (n * (n - 1)) / 2;
        std::vector<double> result(total_pairs);

        // Create threads
        pthread_t* threads = new pthread_t[num_threads];
        ThreadData* thread_data = new ThreadData[num_threads];

        // Calculate rows per thread (n-1 total rows to process)
        size_t total_rows = n - 1;
        size_t rows_per_thread = total_rows / num_threads;
        size_t extra_rows = total_rows % num_threads;

        size_t current_row = 0;

        // Launch threads
        for (int i = 0; i < num_threads; i++)
        {
            thread_data[i].datasets = &datasets;
            thread_data[i].result = &result;
            thread_data[i].start_row = current_row;
            thread_data[i].n = n;
            
            // Distribute extra rows to first threads
            size_t rows_for_this_thread = rows_per_thread + (i < extra_rows ? 1 : 0);
            thread_data[i].end_row = current_row + rows_for_this_thread;
            
            current_row = thread_data[i].end_row;

            pthread_create(&threads[i], nullptr, compute_correlations_thread, &thread_data[i]);
        }

        // Wait for all threads to complete
        for (int i = 0; i < num_threads; i++)
        {
            pthread_join(threads[i], nullptr);
        }

        // Clean up
        delete[] threads;
        delete[] thread_data;

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