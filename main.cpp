#include "som.hpp"

#include <boost/timer.hpp>

int main (void)
{
    boost::timer timer;

    const int N=49; // number of elements per dimension
    const int S=3;  // number of dimensions

    som::map<N,S> map(500, N/2, 0.75, 1.0, 0.1);
    map.init();
    map.save_image("map0.png");
    map.load_samples(som::util::load_samples_from_file("samples.txt"));
    std::cout << "start learning" << std::endl;
    timer.restart();
    map.learn();
    std::cout << "t1: " << timer.elapsed() << std::endl;

    map.save_image("map1.png");

    return EXIT_SUCCESS;
}
