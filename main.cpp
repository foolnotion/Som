#include "som.hpp"

#include <boost/timer.hpp>

int main (void)
{
    boost::timer timer;

    const int N = 400; // number of elements per dimension
    const int S = 3;  // length of the input samples

    som::map<N,S> map(1, N/2, 0.8, 1.0, 0.01);
//    timer.restart();
    map.init();
    map.load_samples(som::util::load_samples_from_file("samples.txt"));
//    som::util::save_to_image(map,"map0.png");
//    timer.restart();
    map.learn();
//    std::cout << "t1: " << timer.elapsed() << std::endl;
    som::util::save_to_image(map,"map1.png");

    return EXIT_SUCCESS;
}
