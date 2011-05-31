#include "som.hpp"

#include <boost/timer.hpp>

int main (void)
{
    boost::timer timer;

    unsigned map_size = 200;
    unsigned sample_size = 3; // rgb colors (3 components)

    som::map map(map_size, sample_size, som::norm_2<double>());
    som::point3 p;
    std::cout << "Best matching unit for " << map(2,2,2)->get_weights() << std::endl;
    timer.restart();
    p = map.best_mathing_unit(map(2,2,2)->get_weights());
    std::cout << "position: " << p << std::endl;
    std::cout << "t1: " << timer.elapsed() << std::endl;

    return EXIT_SUCCESS;
}
