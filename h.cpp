/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "matrix.hpp"
#include "ppm.hpp"
#include "filters.hpp"
#include <cstdlib>
#include <iostream>
#include <omp.h>

int main(int argc, char const* argv[])
{
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " [radius] [infile] [outfile] [threads]" << std::endl;
        std::exit(1);
    }

    PPM::Reader reader {};
    PPM::Writer writer {};

    int num_threads = std::atoi(argv[4]);
    omp_set_num_threads(num_threads);

    auto m { reader(argv[2]) };
    auto radius { static_cast<unsigned>(std::stoul(argv[1])) };

    auto blurred { Filter::blur(m, radius) };
    writer(blurred, argv[3]);

    return 0;
}


-----


/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>

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

    Matrix blur(Matrix m, const int radius)
    {
        Matrix scratch{PPM::max_dimension};
        auto dst{m};
        
        const auto x_size = dst.get_x_size();
        const auto y_size = dst.get_y_size();

        double w[Gauss::max_radius]{};
        Gauss::get_weights(radius, w);
        
        // Horizontal pass - USE = instead of {}
        #pragma omp parallel for schedule(dynamic)
        for (int x = 0; x < x_size; x++)  // Changed: auto x{0} → int x = 0
        {
            for (int y = 0; y < y_size; y++)  // Changed: auto y{0} → int y = 0
            {
                auto r{w[0] * dst.r(x, y)}, g{w[0] * dst.g(x, y)}, b{w[0] * dst.b(x, y)}, n{w[0]};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto x2{x - wi};
                    if (__builtin_expect(x2 >= 0, 1))
                    {
                        r += wc * dst.r(x2, y);
                        g += wc * dst.g(x2, y);
                        b += wc * dst.b(x2, y);
                        n += wc;
                    }
                    x2 = x + wi;
                    if (__builtin_expect(x2 < x_size, 1))
                    {
                        r += wc * dst.r(x2, y);
                        g += wc * dst.g(x2, y);
                        b += wc * dst.b(x2, y);
                        n += wc;
                    }
                }
                auto inv_n = 1.0 / n;
                scratch.r(x, y) = r * inv_n;
                scratch.g(x, y) = g * inv_n;
                scratch.b(x, y) = b * inv_n;
            }
        }

        // Vertical pass - USE = instead of {}
        #pragma omp parallel for schedule(dynamic)
        for (int y = 0; y < y_size; y++)  // Changed: auto y{0} → int y = 0
        {
            for (int x = 0; x < x_size; x++)  // Changed: auto x{0} → int x = 0
            {
                auto r{w[0] * scratch.r(x, y)}, g{w[0] * scratch.g(x, y)}, b{w[0] * scratch.b(x, y)}, n{w[0]};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto y2{y - wi};
                    if (__builtin_expect(y2 >= 0, 1))
                    {
                        r += wc * scratch.r(x, y2);
                        g += wc * scratch.g(x, y2);
                        b += wc * scratch.b(x, y2);
                        n += wc;
                    }
                    y2 = y + wi;
                    if (__builtin_expect(y2 < y_size, 1))
                    {
                        r += wc * scratch.r(x, y2);
                        g += wc * scratch.g(x, y2);
                        b += wc * scratch.b(x, y2);
                        n += wc;
                    }
                }
                auto inv_n = 1.0 / n;
                dst.r(x, y) = r * inv_n;
                dst.g(x, y) = g * inv_n;
                dst.b(x, y) = b * inv_n;
            }
        }

        return dst;
    }

}
