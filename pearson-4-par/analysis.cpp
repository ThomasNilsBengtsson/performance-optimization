/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <chrono>
#include <atomic>

namespace Analysis
{
    struct ThreadData {
        std::vector<Vector>* datasets;
        std::vector<double>* result;
        std::atomic<size_t>* current_row;
        size_t n;
        size_t chunk_size;
    };

    void* compute_correlations_thread(void* arg)
    {
        ThreadData* data = static_cast<ThreadData*>(arg);
        
        auto start = std::chrono::high_resolution_clock::now();
        size_t pairs_computed = 0;
        
        while (true)
        {
            // Grab a CHUNK of rows at once
            size_t start_row = data->current_row->fetch_add(data->chunk_size);
            
            if (start_row >= data->n - 1)
                break;
            
            size_t end_row = std::min(start_row + data->chunk_size, data->n - 1);
            
            // Process all rows in this chunk
            for (size_t sample1 = start_row; sample1 < end_row; sample1++)
            {
                size_t row_start = sample1 * data->n - (sample1 * (sample1 + 1)) / 2 - sample1;
                
                for (size_t sample2 = sample1 + 1; sample2 < data->n; sample2++)
                {
                    size_t idx = row_start + sample2 - 1;
                    (*data->result)[idx] = pearson((*data->datasets)[sample1], 
                                                   (*data->datasets)[sample2]);
                    pairs_computed++;
                }
            }
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        std::cerr << "Thread computed " << pairs_computed << " pairs in " << duration << "ms" << std::endl;
        
        return nullptr;
    }

    std::vector<double> correlation_coefficients(std::vector<Vector> datasets, int num_threads)
    {
        size_t n = datasets.size();
        size_t total_pairs = (n * (n - 1)) / 2;
        std::vector<double> result(total_pairs);

        std::atomic<size_t> current_row(0);

        pthread_t* threads = new pthread_t[num_threads];
        ThreadData* thread_data = new ThreadData[num_threads];

        // Chunk size: aim for ~20-50 chunks total
        // For 1024 datasets: chunk_size = 1024 / (4 threads × 10) = ~25 rows per chunk
        size_t chunk_size = std::max(size_t(1), (n - 1) / (num_threads * 10));

        for (int i = 0; i < num_threads; i++)
        {
            thread_data[i].datasets = &datasets;
            thread_data[i].result = &result;
            thread_data[i].current_row = &current_row;
            thread_data[i].n = n;
            thread_data[i].chunk_size = chunk_size;

            pthread_create(&threads[i], nullptr, compute_correlations_thread, &thread_data[i]);
        }

        for (int i = 0; i < num_threads; i++)
        {
            pthread_join(threads[i], nullptr);
        }

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