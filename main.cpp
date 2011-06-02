#include "som.hpp"

#include <boost/timer.hpp>

int main (void)
{
    boost::timer timer;

    const int N=21; // number of elements per dimension
    const int S=3;  // number of dimensions

    som::map<N,S> map(1000, N/2, 0.75, 1.0, 0.01);
    map.init();
    map.load_samples(som::util::load_samples_from_file("samples.txt"));
    map.save_image("map0.png");
    timer.restart();
    map.learn();
    std::cout << "t1: " << timer.elapsed() << std::endl;
    map.save_image("map1.png");

    return EXIT_SUCCESS;
}
