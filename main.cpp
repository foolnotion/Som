#include "som.hpp"

#include <boost/timer.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <SFML/Graphics.hpp>

template <int N, int S>
void
save_map(const som::map<N,S>& map, const std::string& filename)
{
    sf::Image image(400,300, sf::Color::Black);
    int x = 1, y = 1;
    for (int i = 0; i != N; ++i)
        for (int j = 0; j != N; ++j)
            for (int k = 0; k != N; ++k)
            {
                ublas::c_vector<double,3> v = map.element(i,j,k);
                sf::Color c(v[0]*255, v[1]*255, v[2]*255);
                image.SetPixel(x, y, c);
                ++x;
                if (x == 400)
                {
                    x = 1;
                    ++y;
                }
            }
    image.SaveToFile(filename);
}

int main (void)
{
    boost::timer timer;
    som::util::samples samples = som::util::load_samples_from_file("samples1.txt");
    const int N=49, S=3;
    som::map<N,S> map(100, N/2, 0.8, 1.0, 0.1);
    map.init();

    save_map<N,S>(map,"map0.png");

    map.load_samples(samples);
    std::cout << "start learning" << std::endl;
    timer.restart();
    map.learn();
    std::cout << "t1: " << timer.elapsed() << std::endl;

    save_map<N,S>(map,"map1.png");

    return EXIT_SUCCESS;
}
