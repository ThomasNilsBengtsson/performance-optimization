/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>

struct Data
{
    int start_row;
    int end_row;
    int y_size;
    int x_size;
    int radius;
    double *weights;
    unsigned char *dstR;
    unsigned char *dstG;
    unsigned char *dstB;
    unsigned char *scratchR;
    unsigned char *scratchG;
    unsigned char *scratchB;
};

namespace Filter
{

    namespace Gauss
    {
        void get_weights(int n, double *weights_out)
        {
            for (auto i = 0; i <= n; i++)
            {
                double x{static_cast<double>(i) * max_x / n};
                weights_out[i] = exp(-x * x * pi);
            }
        }
    }

    Matrix blur(Matrix m, const int radius, const int num_threads)
    {
        Matrix scratch{PPM::max_dimension};
        auto dst{m};

        auto x_size = dst.get_x_size();
        auto y_size = dst.get_y_size();

        double w[Gauss::max_radius]{};
        Gauss::get_weights(radius, w); // Flyttat denna

        auto *dstR = dst.get_R();
        auto *dstG = dst.get_G();
        auto *dstB = dst.get_B();
        auto *scratchR = scratch.get_R();
        auto *scratchG = scratch.get_G();
        auto *scratchB = scratch.get_B();

        int rowsPerThread = y_size / num_threads;
        pthread_t threads[num_threads];
        Data threadData[num_threads];

        for (int i = 0; i < num_threads; i++)
        {
            int startRow = i * rowsPerThread;
            threadData[i].start_row = i * rowsPerThread;
            threadData[i].end_row = (i == num_threads - 1) ? y_size : startRow + rowsPerThread;
            threadData[i].x_size = x_size;
            threadData[i].y_size = y_size;
            threadData[i].radius = radius;
            threadData[i].weights = w;
            threadData[i].dstR = dstR;
            threadData[i].dstG = dstG;
            threadData[i].dstB = dstB;
            threadData[i].scratchR = scratchR;
            threadData[i].scratchG = scratchG;
            threadData[i].scratchB = scratchB;
            pthread_create(&threads[i], NULL, horizontalBlur, &threadData[i]);
        }
        for (int i = 0; i < num_threads; i++)
        {
            pthread_join(threads[i], NULL);
        }
        for (int i = 0; i < num_threads; i++)
        {
            int startRow = i * rowsPerThread;
            threadData[i].start_row = i * rowsPerThread;
            threadData[i].end_row = (i == num_threads - 1) ? y_size : startRow + rowsPerThread;
            threadData[i].x_size = x_size;
            threadData[i].y_size = y_size;
            threadData[i].radius = radius;
            threadData[i].weights = w;
            threadData[i].dstR = dstR;
            threadData[i].dstG = dstG;
            threadData[i].dstB = dstB;
            threadData[i].scratchR = scratchR;
            threadData[i].scratchG = scratchG;
            threadData[i].scratchB = scratchB;
            pthread_create(&threads[i], NULL, verticalBlur, &threadData[i]);
        }
        for (int i = 0; i < num_threads; i++)
        {
            pthread_join(threads[i], NULL);
        }

        return dst;
    }
    void *horizontalBlur(void *arg)
    {

        Data *data = (Data *)arg;

        int radius = data->radius;
        int x_size = data->x_size;
        int start_row = data->start_row;
        int end_row = data->end_row;
        double *w = data->weights;
        unsigned char *dstR = data->dstR;
        unsigned char *dstG = data->dstG;
        unsigned char *dstB = data->dstB;
        unsigned char *scratchR = data->scratchR;
        unsigned char *scratchG = data->scratchG;
        unsigned char *scratchB = data->scratchB;

        for (auto y = start_row; y < end_row; y++)
        {
            for (auto x = 0; x < x_size; x++)
            {

                // unsigned char Matrix::r(unsigned x, unsigned y) const
                // {
                //     return R[y * x_size + x];
                // }
                auto idx = y * x_size + x;
                auto r{w[0] * dstR[idx]};
                auto g{w[0] * dstG[idx]};
                auto b{w[0] * dstB[idx]};
                auto n{w[0]};

                for (auto wi = 1; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto x2{x - wi};
                    auto idx = y * x_size + x2;
                    if (__builtin_expect(x2 >= 0, 1))
                    {

                        r += wc * dstR[idx];
                        g += wc * dstG[idx];
                        b += wc * dstB[idx];
                        n += wc;
                    }
                    x2 = x + wi;
                    auto idx2 = y * x_size + x2;
                    if (__builtin_expect(x2 < x_size, 1))
                    {
                        r += wc * dstR[idx2];
                        g += wc * dstG[idx2];
                        b += wc * dstB[idx2];
                        n += wc;
                    }
                }
                scratchR[idx] = r / n;
                scratchG[idx] = g / n;
                scratchB[idx] = b / n;
            }
        }
        return NULL;
    }

    void *verticalBlur(void *arg)
    {
        Data *data = (Data *)arg;

        int radius = data->radius;
        int x_size = data->x_size;
        int y_size = data->y_size;
        int start_row = data->start_row;
        int end_row = data->end_row;
        double *w = data->weights;
        unsigned char *dstR = data->dstR;
        unsigned char *dstG = data->dstG;
        unsigned char *dstB = data->dstB;
        unsigned char *scratchR = data->scratchR;
        unsigned char *scratchG = data->scratchG;
        unsigned char *scratchB = data->scratchB;
        for (auto x = 0; x < x_size; x++)
        {
            for (auto y = start_row; y < end_row; y++)
            {
                /* Den kallar get weights andra gången här vilket inte behövs
                double w[Gauss::max_radius]{};
                Gauss::get_weights(radius, w);
                */
                auto idx = y * x_size + x;
                auto r{w[0] * scratchR[idx]};
                auto g{w[0] * scratchG[idx]};
                auto b{w[0] * scratchB[idx]};
                auto n{w[0]};

                for (auto wi = 1; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto y2{y - wi};
                    auto idx = y2 * x_size + x;
                    if (__builtin_expect(y2 >= 0, 1))
                    {
                        r += wc * scratchR[idx];
                        g += wc * scratchG[idx];
                        b += wc * scratchB[idx];
                        n += wc;
                    }
                    y2 = y + wi;
                    auto idx2 = y2 * x_size + x;
                    if (__builtin_expect(y2 < y_size, 1))
                    {
                        r += wc * scratchR[idx2];
                        g += wc * scratchG[idx2];
                        b += wc * scratchB[idx2];
                        n += wc;
                    }
                }
                dstR[idx] = r / n;
                dstG[idx] = g / n;
                dstB[idx] = b / n;
            }
        }
        return NULL;
    }
}
