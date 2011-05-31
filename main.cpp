#include "som.hpp"

#include <boost/timer.hpp>

int main (void)
{
    boost::timer timer;

    som::map map(200, 3);
    som::point3 p;
    std::cout << "Best matching unit for " << map(2,2,2) << std::endl;
    timer.restart();
    p = map.best_matching_unit(map(2,2,2));
    std::cout << "position: " << p << std::endl;
    std::cout << "t1: " << timer.elapsed() << std::endl;

    return EXIT_SUCCESS;
}
