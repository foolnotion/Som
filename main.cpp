#include "som.hpp"

#include <boost/timer.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

int main (void)
{
    boost::timer timer;
    som::util::samples samples = som::util::load_samples_from_file("samples.txt");
    som::map<200, 3> map;
    map.load_samples(samples);
    som::point3 p;
    std::cout << "Best matching unit for " << map.samples_[0];
    timer.restart();
    p = map.best_matching_unit(map.samples_[1]);
    std::cout << " is " << map(p.x, p.y, p.z) << " on position: " << p
              << " norm = " << ublas::norm_2(map(p.x,p.y,p.z)-map.samples_[1]) << std::endl;
    std::cout << "t1: " << timer.elapsed() << std::endl;

    return EXIT_SUCCESS;
}
